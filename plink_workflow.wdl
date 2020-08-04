import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/v0.01-alpha/tasks/preprocess_tasks.wdl" as preprocess_tasks

# import "tasks/preprocess_tasks.wdl" as preprocess_tasks

workflow run_preprocess {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    File genotype_samples_to_keep_file
    File imputed_samples_to_keep_file
    File covariate_tsv_file

    String covar_sampleID_colname
    
    File? chain_file
    File? imputed_list_of_vcf_file
    File? imputed_list_of_bgen_file
    #File? imputed_list_of_bgen_index_file

    call preprocess_tasks.get_cohort_samples {
        input: 
            covariate_tsv_file = covariate_tsv_file,
            genotype_samples_to_keep_file = genotype_samples_to_keep_file,
            imputed_samples_to_keep_file = imputed_samples_to_keep_file
    }

    call preprocess_tasks.plink_subset_sample {
        input:
             genotype_bed = genotype_bed, 
             genotype_bim = genotype_bim,
             genotype_fam = genotype_fam,
             samples_to_keep_file = get_cohort_samples.plink_subset_samples
    }

    call preprocess_tasks.run_genotype_qc_filter {
        input:
            genotype_bed = plink_subset_sample.genotype_subsetSample_bed,
            genotype_bim = plink_subset_sample.genotype_subsetSample_bim,
            genotype_fam = plink_subset_sample.genotype_subsetSample_fam

    }

    call preprocess_tasks.run_ld_prune {
        input:
            genotype_bed = run_genotype_qc_filter.gentoype_qc_filtered_bed,
            genotype_bim = run_genotype_qc_filter.gentoype_qc_filtered_bim,
            genotype_fam = run_genotype_qc_filter.gentoype_qc_filtered_fam
    }

    call preprocess_tasks.plink_pca {
        input:
    	    genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam
            
    }

    call preprocess_tasks.addPCs_to_covar_matrix {
         input:
             covar_file = covariate_tsv_file, 
             plink_pca_eigenvec_file = plink_pca.genotype_pruned_pca_eigenvec,
             covar_sampleID_colname = covar_sampleID_colname
    }
    


    if (defined(chain_file)) {

        call preprocess_tasks.liftover_plink {
            input:
    	    genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam,
            chain_file = select_first([chain_file, "null"])
        }
    }

    
    
    if (defined(imputed_list_of_vcf_file)) {

        Array[Array[File]] imputed_files = read_tsv(select_first([imputed_list_of_vcf_file,"null"]))

        scatter (imputed_file in imputed_files) {
	        call preprocess_tasks.vcf_to_bgen {
	            input:
                    vcf_file = imputed_file[0],
                    samples_to_keep_file = imputed_samples_to_keep_file
		        }
	    }

        call preprocess_tasks.cat_file{
            input:
                files = vcf_to_bgen.bgen_file_paths
        }
    } 


	output {
          File genotype_ready_bed = select_first([liftover_plink.output_bed, run_ld_prune.genotype_pruned_bed])
          File genotype_ready_bim = select_first([liftover_plink.output_bim, run_ld_prune.genotype_pruned_bim])
          File genotype_ready_fam = select_first([liftover_plink.output_fam, run_ld_prune.genotype_pruned_fam])
          #File genotype_pruned_pca_eigenvec = plink_pca.genotype_pruned_pca_eigenvec
          #Array[File] bgen_files = select_first([vcf_to_bgen.out_bgen, imputed_bgen_files])
          #Array[File] bgen_file_indices = select_first([index_bgen_file.bgen_file_index, imputed_bgen_index_files])
          #Array[Array[File]] bgen_files_and_indices = select_first([converted_bgen_file_list,read_tsv(select_first([imputed_list_of_bgen_file,"null"]))])          
          File bgen_paths_file = select_first([imputed_list_of_bgen_file,cat_file.merged_file])
          File bgen_subset_samples = get_cohort_samples.bgen_subset_samples
          File plink_subset_samples = get_cohort_samples.plink_subset_samples
          File covar_file = addPCs_to_covar_matrix.covar_file_with_pcs
          File pcs_as_string_file = addPCs_to_covar_matrix.pcs_as_string_file
 	}
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Preprocess Genotype files for biobank scale GWAS study."
	}

    

}

