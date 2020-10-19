version development

import "./run_association_test_on_subsets_workflow.wdl" as run_association_test 


workflow run_meta_analysis {
    input {
        File batch_tsv_file 
        File covariate_tsv_file
        File variable_info_tsv_file
        File sample_sets_json_file

        Boolean? useBOLT
        Boolean? useSAIGE

        #for bolt
        File? genetic_map_file
        File? ld_scores_file
    
        Float? minMAF=0.001
        Float? minMAC=1

        File? chain_file

        String? id_delim #delim character for vcf file, if not defined, double ID is assumed
    }

    #read in batch information. If there are more than 1 batch we will meta analyze the summary results
    Array[Array[String]] batch_tsv = read_tsv(batch_tsv_file) 
    Int n_batches = length(batch_tsv)

    scatter (idx in range(n_batches-1)) { 
        Array[String] batch_tsv_rows = batch_tsv[(idx)] 
    }
    Map[String, Array[String]] batch_tbl = as_map(zip(batch_tsv[0], transpose(batch_tsv_rows)))

    scatter (idx in range(n_batches-1)) {
        File genotype_bed = batch_tbl["genotype_bed"][idx]
        File genotype_bim = batch_tbl["genotype_bim"][idx]
        File genotype_fam = batch_tbl["genotype_fam"][idx]
        File genotype_samples_to_keep_file = batch_tbl["genotype_samples_to_keep_file"][idx]
        File imputed_samples_to_keep_file = batch_tbl["imputed_samples_to_keep_file"][idx]
        File imputed_list_of_files = batch_tbl["imputed_list_of_files"][idx]
        String impute_file_format = batch_tbl["impute_file_format"][idx]

        String batch_name = batch_tbl["batch_name"][idx] 

        call run_association_test.run_association_test {
            input:
                genotype_bed = genotype_bed,
                genotype_bim = genotype_bim,
                genotype_fam = genotype_fam,
                genotype_samples_to_keep_file = genotype_samples_to_keep_file,
                imputed_samples_to_keep_file = imputed_samples_to_keep_file,
                covariate_tsv_file = covariate_tsv_file,
                variable_info_tsv_file = variable_info_tsv_file,
                sample_sets_json_file = sample_sets_json_file,
                imputed_list_of_files = imputed_list_of_files,
                impute_file_format = impute_file_format,

                useBOLT = useBOLT,
                useSAIGE = useSAIGE,

                #for bolt
                genetic_map_file = genetic_map_file,
                ld_scores_file = ld_scores_file,
    
                minMAF = minMAF,
                minMAC = minMAC,

                chain_file = chain_file,

                id_delim = id_delim

        }
    }    
}

task run_metal {
    input {
    Array[File] association_summary_files
    String prefix 

    String? FREQLABEL = "Freq_Tested_Allele_in_TOPMed"
    String? MARKERLABEL = "SNP"
   

    Int? memory = 32
    Int? disk = 20
    Int? threads = 2
    Int? preemptible_tries = 3
    }

     command <<<
        metal_command="MARKERLABEL SNP\n \
                        ALLELELABELS Tested_Allele Other_Allele \
                        EFFECTLABEL BETA
                        STDERRLABEL SE
                        FREQLABEL Freq_Tested_Allele_in_TOPMed
                        CUSTOMVARIABLE TotalSampleSize
                        LABEL TotalSampleSize as N
                        SCHEME STDERR
                        GENOMICCONTROL ON
                        AVERAGEFREQ ON
                        MINMAXFREQ ON\n"


                        PROCESSFILE /mnt/projects/plco/bq_bmi_curr_co.GSA_batch1.boltlmm.rawids.tsv
                        PROCESSFILE /mnt/projects/plco/bq_bmi_curr_co.GSA_batch2.boltlmm.rawids.tsv
                        PROCESSFILE /mnt/projects/plco/bq_bmi_curr_co.Oncoarray.boltlmm.rawids.tsv
     
        metal_command=${metal_command} +  "OUTFILE " + ~{prefix} +  ".metal.tsv\n \
                     ANALYZE HETEROGENEITY\n \
                     QUIT"


        metal < ${metal_command}
    >>>

    runtime {
		docker: "quay.io/shukwong/metal:2011-03-25"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: "${preemptible_tries}"
	}

    output {
        File metal_output_file =  "${prefix}.metal.tsv"
    }
}