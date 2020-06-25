workflow run_preprocess {
    
	File genotype_bed
	File genotype_bim
	File genotype_fam
    File pheno_file
    File ld_scores_file
    File genetic_map_file
    File imputed_sample_file
    File covar_file
    File imputed_bgen_filelist

    String pheno_col
    String qCovarCol 

    call run_bolt {
        input:
            genotype_bed = genotype_bed,
	        genotype_bim = genotype_bim,
	        genotype_fam = genotype_fam,
            pheno_file = pheno_file,
            ld_scores_file = ld_scores_file,
            genetic_map_file = genetic_map_file,
            imputed_sample_file = imputed_sample_file,
            covar_file = covar_file,
            imputed_bgen_filelist = imputed_bgen_filelist,
            pheno_col = pheno_col,
            qCovarCol = qCovarCol
    }

    output {
		File genotype_stats_file = run_bolt.genotype_stats_file
        File imputed_stats_file = run_bolt.imputed_stats_file
	}
}    

task run_bolt {

	File genotype_bed
	File genotype_bim
	File genotype_fam
    File pheno_file
    File ld_scores_file
    File genetic_map_file
    File imputed_sample_file
    File covar_file
    File imputed_bgen_filelist

    String pheno_col
    String qCovarCol 

    Array[Array[String]] imputed_bgen_files = read_tsv(imputed_bgen_filelist)
    Array[String] bgen_files_with_input_option = prefix("--bgenFile ", imputed_bgen_files[0]) 
    
	Int? memory = 32
	Int? disk = 200
    Int? threads = 32

	command {
        bolt \
            --bed=${genotype_bed} \
            --bim=${genotype_bim} \
            --fam=${genotype_fam} \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno_col} \
            --LDscoresFile=${ld_scores_file} \
            --geneticMapFile=${genetic_map_file} \
            --lmmForceNonInf \
            --numThreads=${threads} \
            --covarFile=${covar_file} \
            --qCovarCol=${qCovarCol} \
            --statsFile=bolt_genotype_stats_${pheno_col}.gz \
            ${sep=" " bgen_files_with_input_option} \
            --bgenMinMAF=0.01 \
            --bgenMinINFO=0.3 \
            --sampleFile=${imputed_sample_file} \
            --statsFileBgenSnps=bolt_imputed_stats_${pheno_col}.gz \
            --noBgenIDcheck \
            --verboseStats
	}

	runtime {
		docker: "quay.io/shukwong/bolt-lmm:2.3.4"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

	output {
		File genotype_stats_file = "bolt_genotype_stats_${pheno_col}.gz"
        File imputed_stats_file = "bolt_imputed_stats_${pheno_col}.gz"
	}
}