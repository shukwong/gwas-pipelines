task run_bolt {

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

    String pheno_col
    String genotype_stats_filename = bolt_genotype_stats_${pheno_col}.gz
    String imputed_stats_filename = bolt_imputed_stats_${pheno_col}.gz

	Int? memory = 360
	Int? disk = 500
    Int? threads = 64

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
            --qCovarCol=PC{1:10} \
            --statsFile=${genotype_stats_filename} \
            --bgenFile=${imputed_bgen_file} \
            --bgenMinMAF=0.01 \
            --bgenMinINFO=0.3 \
            --sampleFile=${imputed_sample_file} \
            --statsFileBgenSnps=${imputed_stats_filename} \
            --noBgenIDcheck \
            --verboseStats
	}

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: ${threads}
		gpu: false
	}

	output {
		File genotype_stats_file = "${genotype_stats_filename}"
        File imputed_stats_file = "${imputed_stats_filename}"
	}
}