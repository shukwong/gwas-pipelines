import "plink_workflow.wdl" as preprocess
import "tasks/preprocess_tasks.wdl" as preprocess_tasks
import "tasks/saige_workflow.wdl" as saige

workflow run_assoication_test {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    File genotype_samples_to_keep_file
    File imputed_samples_to_keep_file

    File covariate_tsv_file
    File variable_info_tsv_file
    File sample_sets_json_file

    String phenoCol
    
    File? imputed_list_of_vcf_file
    File? imputed_list_of_bgen_file
    File? chain_file
    
    call  preprocess_tasks.get_covar_subsets {
        input:
        covariate_tsv_file = covariate_tsv_file, 
        variable_info_tsv_file = variable_info_tsv_file, 
        sample_sets_json_file = sample_sets_json_file,
        phenoCol = phenoCol
    }

    Array[String] binary_covar_list_lines  = read_lines(get_covar_subsets.binary_covar_list_file)
    String binary_covar_list = binary_covar_list_lines[0]
    Array[String] continuous_covar_list_lines = read_lines(get_covar_subsets.continuous_covar_list_file)
    String continuous_covar_list = continuous_covar_list_lines[0]

    scatter (covar_subset_file in get_covar_subsets.covar_subsets_files) {
        call preprocess.run_preprocess {
            input:
            genotype_bed = genotype_bed,
            genotype_bim = genotype_bim,
            genotype_fam = genotype_fam,

            genotype_samples_to_keep_file = genotype_samples_to_keep_file,
            imputed_samples_to_keep_file = imputed_samples_to_keep_file,
            #imputed_list_of_vcf_file = imputed_list_of_vcf_file,
            imputed_list_of_bgen_file = imputed_list_of_bgen_file,
            covariate_tsv_file = covar_subset_file,
            chain_file = chain_file
        }

        call saige.run_saige {
            input:
                genotype_bed = run_preprocess.genotype_ready_bed,
                genotype_bim = run_preprocess.genotype_ready_bim,
                genotype_fam = run_preprocess.genotype_ready_fam,
                bgen_paths_file = run_preprocess.bgen_paths_file,
                imputed_samples_file = imputed_samples_to_keep_file,
                phenoCol = phenoCol,
                covar_file = run_preprocess.covar_subset_file,
                covarColList = binary_covar_list + "," + continuous_covar_list

        }
    }

	output {
        Array[File] genotype_ready_bed = run_preprocess.genotype_ready_bed
        Array[File] genotype_ready_bim = run_preprocess.genotype_ready_bim
        Array[File] genotype_ready_fam = run_preprocess.genotype_ready_fam
        Array[File] merged_saige_file = run_saige.merged_saige_file
 	}
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Biobank scale association study with mutiple subsets (sample target variable)."
	}

}

