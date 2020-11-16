version development

import "./run_association_test_workflow.wdl" as run_association_test 
import "./tasks/preprocess_workflow.wdl" as gwas_tasks

workflow run_meta_analysis {
    input {
        File batch_tsv_file 
        File covariate_tsv_file
        File variable_info_tsv_file

        String binary_covar_list
        String continuous_covar_list 
        String phenoCol 
        String covar_sampleID_colname
        String phenotype_type
        String setname 

        Boolean? useBOLT
        Boolean? useSAIGE

        #for bolt
        File? genetic_map_file
        File? ld_scores_file
    
        Float? minMAF=0.001
        Float? minMAC=1

        File? chain_file

        String? id_delim #delim character for vcf file, if not defined, double ID is assumed
        String? dosageField
    }

    #read in batch information. If there are more than 1 batch we will meta analyze the summary results
    Array[Array[String]] batch_tsv = read_tsv(batch_tsv_file) 
    Int n_batches = length(batch_tsv)-1

    scatter (idx in range(n_batches)) { 
        Array[String] batch_tsv_rows = batch_tsv[(idx+1)] 
    }
    Map[String, Array[String]] batch_tbl = as_map(zip(batch_tsv[0], transpose(batch_tsv_rows)))

    scatter (idx in range(n_batches)) {
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
                imputed_list_of_files = imputed_list_of_files,
                impute_file_format = impute_file_format,

                batch_name = batch_name,

                useBOLT = useBOLT,
                useSAIGE = useSAIGE,

                #for bolt
                genetic_map_file = genetic_map_file,
                ld_scores_file = ld_scores_file,
    
                minMAF = minMAF,
                minMAC = minMAC,

                chain_file = chain_file,

                id_delim = id_delim,

                binary_covar_list = binary_covar_list,
                continuous_covar_list = continuous_covar_list,
                phenoCol = phenoCol,
                covar_sampleID_colname = covar_sampleID_colname,
                phenotype_type = phenotype_type,
                setname = setname,
                dosageField = dosageField
        }
    }  

    if (defined(useBOLT) && useBOLT ) {
        call run_metal as run_metal_bolt {
            input:
                association_summary_files = run_association_test.merged_bolt_file, 
                prefix = setname,
                FREQLABEL = "A1FREQ"
        }

        call gwas_tasks.make_summary_plots as make_bolt_plots {
            input: 
                association_summary_file = run_metal_bolt.metal_output_file,
                BP_column = "POS",
                CHR_column = "CHR",
                pval_col = "P-value",
                minrep_col = "MarkerName",
        } 
    }

    if (defined(useSAIGE) && useSAIGE ) {
        call run_metal as run_metal_saige {
            input:
                association_summary_files = run_association_test.merged_saige_file, 
                prefix = setname,
                FREQLABEL = ""
        }

        call gwas_tasks.make_summary_plots as make_sage_plots {
            input: 
                association_summary_file = run_metal_saige.metal_output_file,
                BP_column = "POS",
                CHR_column = "CHR",
                pval_col = "P-value",
                minrep_col = "MarkerName"
        } 

    }

    output {
        File? bolt_metal_output_file =  run_metal_bolt.metal_output_file
        File? bolt_metal_info_file = run_metal_bolt.metal_info_file
        File? bolt_metal_manhattan_file = make_bolt_plots.manhattan_file
        File? bolt_metal_qqplot_file = make_bolt_plots.qqplot_file
        File? saige_metal_output_file =  run_metal_saige.metal_output_file
        File? saige_metal_info_file = run_metal_saige.metal_info_file
    }
}

#TODO: add allele frequency for saige
task run_metal {
    input {
    Array[File] association_summary_files
    String prefix 

    String FREQLABEL
    String? MARKERLABEL = "SNP"
   

    Int? memory = 32
    Int? disk = 20
    Int? threads = 2
    Int? preemptible_tries = 3
    }

     command <<<
        set -euo pipefail

        wget https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/scripts/process_metal_outputs.R

        echo -e "MARKERLABEL ~{MARKERLABEL} \n \
ALLELELABELS Tested_Allele Other_Allele \n \
EFFECTLABEL BETA \n \
STDERRLABEL SE \n \
CUSTOMVARIABLE TotalSampleSize \n \
LABEL TotalSampleSize as N \n \
SCHEME STDERR \n \
GENOMICCONTROL ON \n \
PROCESSFILE ~{sep=',PROCESSFILE ' association_summary_files}  \n\
OUTFILE ~{prefix}.metal .tsv \n\
ANALYZE HETEROGENEITY \n\
QUIT" > metal_command_temp

        cat  metal_command_temp | tr ',' '\n' >metal_command

        metal metal_command

        echo -e "MarkerName\tCHR\tPOS\tAllele1\tAllele2\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal\tTotalSampleSize">~{prefix}.metal.tsv
        grep -v MarkerName ~{prefix}.metal1.tsv | tr ':' '\t' | awk '{print $1":"$2":"$3":"$4"\t"$0}' | cut -f1-3,6-16 >>~{prefix}.metal.tsv

        gzip ~{prefix}.metal.tsv
        mv ~{prefix}.metal1.tsv.info ~{prefix}.metal.tsv.info
        
    >>>

    runtime {
		docker: "quay.io/shukwong/metal:2011-03-25"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: "${preemptible_tries}"
	}

    output {
        File metal_output_file =  prefix + ".metal.tsv.gz"
        File metal_info_file = prefix + ".metal.tsv.info"
    }
}
