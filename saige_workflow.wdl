
workflow run_saige {
    
    File imputed_bgen_list_file
    File imputed_sample_file
    File covar_file

    File gmmat_model_file
    File variance_ratio_file

    String chrom 

   Array[File] imputed_bgen_files = read_tsv(imputed_bgen_list_file)

    scatter (imputed_bgen_file in imputed_bgen_files) {
		call saige_step2_imputed_bgen {
			input:
                vcf_file = imputed_file
		}
	}
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Run SAIGE"
	}

}


task saige_step1 {

	File genotype_bed
	File genotype_bim
	File genotype_fam
	File samples_to_remove_file
    File pheno_file
    File ld_scores_file
    File genetic_map_file
    File imputed_bgen_file
    File imputed_sample_file
    File covar_file
    File output_dir

    String pheno_col
    String covarColList #covar list, separated by comma
    
    String genotype_plink_prefix 

	Int? memory = 360
	Int? disk = 500
    Int? threads = 64

	command {
        step1_fitNULLGLMM.R     \
            --plinkFile=${genotype_plink_prefix} \
            --phenoFile=${covar_file} \
            --phenoCol=j_pros_cancer \
            --covarColList=${covarColList} \
            --sampleIDColinphenoFile=IID \
            --traitType=binary        \
            --outputPrefix=${output_dir}/saige_step1_${phenoCol} \
            --nThreads=${threads}
	}

	runtime {
		docker: "wzhou88/saige:0.38"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

	output {
		File genotype_stats_file = ${genotype_stats_filename}.gz
        File imputed_stats_file = ${imputed_stats_filename}.gz
	}
}


task saige_step2_imputed_bgen {

    File imputed_bgen_file
    File imputed_bgen_file_index
    File imputed_sample_file
    File covar_file

    File gmmat_model_file
    File variance_ratio_file

    String chrom 

    String file_prefix = basename(imputed_bgen_file, ".bgen") 
    String saige_output_file = file_prefix + "." + chrom + ".txt"

	Int? memory = 64
	Int? disk = 500
    Int? threads = 64

	command {
      step2_SPAtests.R \
        --bgenFile=${imputed_bgen_file} \
        --bgenFileIndex=${imputed_bgen_file_index} \
        --IsDropMissingDosages=FALSE \
        --minMAF=0.01 \
        --minMAC=3 \
        --chrom=${chrom} \
        --sampleFile=${imputed_sample_file} \
        --GMMATmodelFile=${gmmat_model_file} \
        --varianceRatioFile=${variance_ratio_file} \
        --SAIGEOutputFile=${saige_output_file} \
        --numLinesOutput=2 \
        --IsOutputNinCaseCtrl=TRUE \
        --IsOutputHetHomCountsinCaseCtrl=TRUE \
         --IsOutputAFinCaseCtrl=TRUE
	}

	runtime {
		docker: "wzhou88/saige:0.38"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: 64
		gpu: false
	}

	output {
		File saige_output_file = "${saige_output_file}.gz"
	}
}