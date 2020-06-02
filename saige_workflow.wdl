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

    genotype_bed=/mnt/data/PLCO_GSA_ld_pruned.bed
    genotype_bim=/mnt/data/PLCO_GSA_ld_pruned.bim
    genotype_fam=/mnt/data/PLCO_GSA_ld_pruned.fam
    genotype_plink_prefix=/mnt/data/PLCO_GSA_ld_pruned
    covar_file=/mnt/data/Phenotype/gsa_qx_v5.with_na.augmented.18may2020.pheno

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
        cpu: 64
		gpu: false
	}

	output {
		File genotype_stats_file = ${genotype_stats_filename}.gz
        File imputed_stats_file = ${imputed_stats_filename}.gz
	}
}