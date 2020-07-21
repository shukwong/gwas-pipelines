
import "tasks/preprocess_tasks.wdl" as preprocess_tasks

workflow run_preprocess {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    File genotype_samples_to_keep_file
    File imputed_samples_to_keep_file
    File imputed_list_of_vcfs_file
    File covar_file
    
    File chain_file

    call preprocess_tasks.get_cohort_samples {
        input: 
            covar_file = covar_file,
            genotype_samples_to_keep_file = genotype_samples_to_keep_file,
            imputed_samples_to_keep_file = imputed_samples_to_keep_file
    }

    call preprocess_tasks.plink_subset_sample {
        input:
             genotype_bed = genotype_bed, 
             genotype_bim = genotype_bim,
             genotype_fam = genotype_fam,
             samples_to_keep_file = genotype_samples_to_keep_file
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

    call preprocess_tasks.liftover_plink_bim {
        input:
    	    genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam,
            chain_file = chain_file
    }

    call preprocess_tasks.liftover_plink {
        input:
    	    genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam,
            liftover_mapped_ids_file = liftover_plink_bim.liftover_mapped_ids_file,
            liftover_mapped_new_bim_file = liftover_plink_bim.liftover_mapped_new_bim_file
           
    }

	Array[Array[File]] imputed_files = read_tsv(imputed_list_of_vcfs_file)
    #Array[File] imputed_files = glob(imputed_files_dir + "/*.dose.vcf.gz")

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

	output {
          File genotype_ready_bed = liftover_plink.output_bed
          File genotype_ready_bim = liftover_plink.output_bim
          File genotype_ready_fam = liftover_plink.output_fam
          File genotype_pruned_pca_eigenvec = plink_pca.genotype_pruned_pca_eigenvec
          Array[File] bgen_files = vcf_to_bgen.out_bgen
          Array[File] bgen_file_indices = index_bgen_file.bgen_file_index
 	}
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Preprocess Genotype files for biobank scale GWAS study."
	}

}

