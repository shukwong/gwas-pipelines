
import "tasks/preprocess_tasks.wdl" as preprocess_tasks

workflow run_preprocess {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    File genotype_samples_to_keep_file
    File imputed_samples_to_keep_file
    File covariate_tsv_file
    
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
            call preprocess_tasks.index_bgen_file {
                input:
                    bgen_file = vcf_to_bgen.out_bgen
            }
	    }

        Array[Array[File]] converted_bgen_file_list = transpose([vcf_to_bgen.out_bgen,  index_bgen_file.bgen_file_index])
        #write_tsv(converted_bgen_file_list)
    } 

    # Array[Array[File]] bgen_files_and_indices_out = if defined(imputed_list_of_bgen_file) then 
    #     read_tsv(select_first([imputed_list_of_bgen_file,"null"])) 
    #     else if (defined(imputed_list_of_vcf_file)) then 
    #     transpose([vcf_to_bgen.out_bgen,  index_bgen_file.bgen_file_index]) else "null"


    # if (defined(imputed_list_of_bgen_file)) {
    #     #Array[File] imputed_bgen_files = read_lines(select_first([imputed_list_of_bgen_file,"null"]))
    #     #Array[File] imputed_bgen_index_files = read_lines(select_first([imputed_list_of_bgen_index_file,"null"]))
    #     Array[Array[File]] bgen_files_and_indices = read_tsv(select_first([imputed_list_of_bgen_file, "null"])) 
    # } else {
    #     Array[File] bgen_files = select_first([vcf_to_bgen.out_bgen, imputed_bgen_files])
    #     Array[File] bgen_file_indices = select_first([index_bgen_file.bgen_file_index, imputed_bgen_index_files]) 
    # }     

    # Array[Array[File]] bgen_files_and_indices = if defined(imputed_list_of_bgen_file) then read_tsv(imputed_list_of_bgen_file)
    # else if (imputed_list_of_vcf_file) 

	output {
          File genotype_ready_bed = select_first([liftover_plink.output_bed, run_ld_prune.genotype_pruned_bed])
          File genotype_ready_bim = select_first([liftover_plink.output_bim, run_ld_prune.genotype_pruned_bim])
          File genotype_ready_fam = select_first([liftover_plink.output_fam, run_ld_prune.genotype_pruned_fam])
          File genotype_pruned_pca_eigenvec = plink_pca.genotype_pruned_pca_eigenvec
          #Array[File] bgen_files = select_first([vcf_to_bgen.out_bgen, imputed_bgen_files])
          #Array[File] bgen_file_indices = select_first([index_bgen_file.bgen_file_index, imputed_bgen_index_files])
          Array[Array[File]] bgen_files_and_indices = select_first([converted_bgen_file_list,read_tsv(select_first([imputed_list_of_bgen_file,"null"]))])          
          File bgen_subset_samples = get_cohort_samples.bgen_subset_samples
          File plink_subset_samples = get_cohort_samples.plink_subset_samples
 	}
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Preprocess Genotype files for biobank scale GWAS study."
	}

    parameter_meta {
        genotype_bed : "plink bed file for direct genotypes"
        genotype_bim : "plink bim file for direct genotypes"
        genotype_fam : "plink fam file for direct genotypes"

        genotype_samples_to_keep_file : "Listing the samples to be kept in the genotype file. Expected first two columns to be FID and IID."
        imputed_samples_to_keep_file : "Listing the samples to be kept in the imputed file. Expected just one column"
        covar_file : "File that contains phenotype and covariates information of the samples"
    
        chain_file : "Optional: A liftover chain file to convert the genome build of the genotype files (in Plink) to another build"
        imputed_list_of_vcf_file : "Optional: A file that contains the file locations of imputed VCFs. They will be converted to bgen. Either this file or the list of bgen and bgen index files are required"
        imputed_list_of_bgen_file : "Optional: A file that contains the file locations of imputed bgen files"
        imputed_list_of_bgen_index_file : "Optional: A file that contains the file locations of imputed bgen index files. This needs to be in same order as the imputed_list_of_bgen_file. Will do it this way until the as_map() function works as expected."
    }

}

