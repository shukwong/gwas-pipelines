
 import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/plink_workflow.wdl" as preprocess
 import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/tasks/saige_workflow.wdl" as saige
 import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/tasks/bolt_workflow.wdl" as bolt

#  import "plink_workflow.wdl" as preprocess
#  import "tasks/saige_workflow.wdl" as saige
#  import "tasks/bolt_workflow.wdl" as bolt

workflow run_association_test {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    File genotype_samples_to_keep_file
    File imputed_samples_to_keep_file

    File covariate_tsv_file
    File variable_info_tsv_file
    File sample_sets_json_file

   # String phenoCol
   # String covar_sampleID_colname

    Boolean? useBOLT
    Boolean? useSAIGE

    #for bolt
    File genetic_map_file
    File ld_scores_file
    
    File? imputed_list_of_vcf_file
    File? imputed_list_of_bgen_file
    File? chain_file

    String? id_delim #delim character for vcf file, if not defined, double ID is assumed
    
    call get_covar_subsets {
        input:
        covariate_tsv_file = covariate_tsv_file, 
        variable_info_tsv_file = variable_info_tsv_file, 
        sample_sets_json_file = sample_sets_json_file
    }

    Array[String] binary_covar_list_lines  = read_lines(get_covar_subsets.binary_covar_list_file)
    String binary_covar_list = binary_covar_list_lines[0]
    Array[String] continuous_covar_list_lines = read_lines(get_covar_subsets.continuous_covar_list_file)
    String continuous_covar_list = continuous_covar_list_lines[0]

    Array[String] phenoCol_lines = read_lines(get_covar_subsets.phenotype_line_file)
    String phenoCol = phenoCol_lines[0]

    Array[String] covar_sampleID_colname_lines = read_lines(get_covar_subsets.sampleid_line_file)
    String covar_sampleID_colname = covar_sampleID_colname_lines[0]

    Array[String] phenotype_type_lines = read_lines(get_covar_subsets.phenotype_type_file)
    String phenotype_type = phenotype_type_lines[0]

    scatter (covar_subset_file in get_covar_subsets.covar_subsets_files) {
        call preprocess.run_preprocess {
            input:
            genotype_bed = genotype_bed,
            genotype_bim = genotype_bim,
            genotype_fam = genotype_fam,

            genotype_samples_to_keep_file = genotype_samples_to_keep_file,
            imputed_samples_to_keep_file = imputed_samples_to_keep_file,
            covariate_tsv_file = covar_subset_file,
            covar_sampleID_colname = covar_sampleID_colname,
            imputed_list_of_vcf_file = imputed_list_of_vcf_file,
            imputed_list_of_bgen_file = imputed_list_of_bgen_file,
            chain_file = chain_file,
            id_delim = id_delim
        }

        Array[String] pcs_as_string_lines = read_lines(run_preprocess.pcs_as_string_file)
        String pcs_as_string = pcs_as_string_lines[0]

        if (defined(useSAIGE) && useSAIGE) {
            call saige.run_saige {
                input:
                    genotype_bed = run_preprocess.genotype_ready_bed,
                    genotype_bim = run_preprocess.genotype_ready_bim,
                    genotype_fam = run_preprocess.genotype_ready_fam,
                    #bgen_paths_file = run_preprocess.bgen_paths_file,
                    bgen_files_and_indices = run_preprocess.bgen_files_and_indices,
                    imputed_samples_file = run_preprocess.bgen_samples,
                    phenoCol = phenoCol,
                    covar_file = run_preprocess.covar_file,
                    covarColList = binary_covar_list + "," + continuous_covar_list + "," + pcs_as_string
            }
        }

        if (defined(useBOLT) && useBOLT ) {
            call bolt.bolt_workflow {
                input: 
                    genotype_bed = run_preprocess.genotype_ready_bed,
	                genotype_bim = run_preprocess.genotype_ready_bim,
	                genotype_fam = run_preprocess.genotype_ready_fam,
                    pheno_file = run_preprocess.covar_file,
                    ld_scores_file = ld_scores_file,
                    genetic_map_file = genetic_map_file,
                    imputed_samples_file = run_preprocess.bgen_samples,
                    covar_file = run_preprocess.covar_file,
                    #bgen_list_file = run_preprocess.bgen_paths_file,
                    bgen_files_and_indices = run_preprocess.bgen_files_and_indices,
                    pheno_col = phenoCol,
                    qCovarCol = continuous_covar_list + "," + pcs_as_string
            }
        }    
    }

	# output {
    #     Array[File] genotype_ready_bed = run_preprocess.genotype_ready_bed
    #     Array[File] genotype_ready_bim = run_preprocess.genotype_ready_bim
    #     Array[File] genotype_ready_fam = run_preprocess.genotype_ready_fam
    #     #Array[File] merged_saige_file = run_saige.merged_saige_file
 	# }
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Biobank scale association study with mutiple subsets (sample target variable)."
	}

    parameter_meta {
        genotype_bed : "plink bed file for direct genotypes"
        genotype_bim : "plink bim file for direct genotypes"
        genotype_fam : "plink fam file for direct genotypes"
        genotype_samples_to_keep_file :  "Listing the samples to be kept in the genotype file. Expected first two columns to be FID and IID."
        imputed_samples_to_keep_file : "Listing the samples to be kept in the imputed file. Expected just one column"
        covariate_tsv_file : "File that contains phenotype and covariates information of the samples"
        variable_info_tsv_file : "variable information"
        sample_sets_json_file : "sample json file to specify sample sets"
        phenoCol : "phenotype (target) column name"
        covar_sampleID_colname : "The sampleID column name for the covar file, default is IID"
        chain_file : "Optional: A liftover chain file to convert the genome build of the genotype files (in Plink) to another build"
        imputed_list_of_vcf_file : "Optional: A file that contains the file locations of imputed VCFs. They will be converted to bgen. Either this file or the list of bgen and bgen index files are required"
        imputed_list_of_bgen_file : "Optional: A file that contains the file locations of imputed bgen files"
        imputed_list_of_bgen_index_file : "Optional: A file that contains the file locations of imputed bgen index files. This needs to be in same order as the imputed_list_of_bgen_file. Will do it this way until the as_map() function works as expected."
        ld_scores_file : "ld scores file for BOLT"
        genetic_map_file : "genetic map file for BOLT"
        useBOLT : "true or false for whether to run BOLT"
        useSAIGE : "true or false for whether to run SAIGE"
    }

}



task get_covar_subsets {
    File covariate_tsv_file 
    File variable_info_tsv_file 
    File sample_sets_json_file

    Int? memory = 4
    Int? disk = 200
    Int? threads = 1
    Int? preemptible_tries = 3

#TODO, change this to git clone a release version when the pipeline is finalized
    command <<< 
       
        wget https://github.com/shukwong/gwas-pipelines/raw/master/scripts/create_covar_files_by_set.R

        Rscript create_covar_files_by_set.R ${covariate_tsv_file} ${variable_info_tsv_file} ${sample_sets_json_file} 
    >>>

    runtime {
		docker: "rocker/tidyverse:4.0.0"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: "${preemptible_tries}"
	}

    output {
        Array[File] covar_subsets_files = glob("covars*.tsv")
        Array[File] covar_subsets_log_files = glob("*.log")
        File binary_covar_list_file = "binary_covar_list.txt"
        File continuous_covar_list_file = "continuous_covar_list.txt"
        File phenotype_line_file = "phenotype_line.txt"
        File phenotype_type_file = "phenotype_type.txt"
        File sampleid_line_file = "sampleid_line.txt"
    }

}
