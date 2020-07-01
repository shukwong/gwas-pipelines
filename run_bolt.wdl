#import "plink_workflow.wdl" as preprocess

workflow bolt_lmm_workflow {
    
	File genotype_bed
	File genotype_bim
	File genotype_fam
    File pheno_file
    File ld_scores_file
    File genetic_map_file
    File imputed_samples_file
    File covar_file
    File bgen_list_file

    String pheno_col
    String qCovarCol #need to figure out the best way to split string so for now the user would have to input things in the format "--qCovarCol=covar1 --qCovarCol=covar2"

    #preprocess 
    # call preprocess.match_genotype_and_imputed_samples {
    #      input:
    #         genotype_bed = genotype_bed,
    #         genotype_bim = genotype_bim,
    #         genotype_fam = genotype_fam,
    #         imputed_samples_file = imputed_samples_file
    # }

    Array[Array[String]] bgen_file_list = read_tsv(bgen_list_file)

    scatter (bgen_file_line in bgen_file_list) {
    
        call run_bolt {
            input:
            # genotype_bed = match_genotype_and_imputed_samples.matched_genotype_bed,
	        # genotype_bim = match_genotype_and_imputed_samples.matched_genotype_bim,
	        # genotype_fam = match_genotype_and_imputed_samples.matched_genotype_fam,
            genotype_bed = genotype_bed,
            genotype_bim = genotype_bim,
            genotype_fam = genotype_fam,
            pheno_file = pheno_file,
            ld_scores_file = ld_scores_file,
            genetic_map_file = genetic_map_file,
            imputed_samples_file = imputed_samples_file,
            covar_file = covar_file,
            bgen_file = bgen_file_line[0],
            pheno_col = pheno_col,
            qCovarCol = qCovarCol
        }
    }

    call combine_bolt_results { input: imputed_stats_files=run_bolt.imputed_stats_file }

    output {
		#File genotype_stats_file = run_bolt.genotype_stats_file
        File imputed_stats_file = combine_bolt_results.merged_stats_file
        Array[File] bolt_lmm_log_files = run_bolt.bolt_lmm_log_file
	}
}    

task combine_bolt_results {

    Array[File] imputed_stats_files
    String pheno_col

    Int? memory = 4
	Int? disk = 100
    Int? threads = 1
    Int? preemptible_tries = 3

    command {
        cat ${sep=' ' imputed_stats_files} >bolt_${pheno_col}_results_merged.gz
    }

    runtime {
		docker: "quay.io/shukwong/bolt-lmm:2.3.4"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: preemptible_tries
	}

    output {
        File merged_stats_file = "bolt_${pheno_col}_results_merged.gz"
    }
}

task run_bolt {

	File genotype_bed
	File genotype_bim
	File genotype_fam
    File pheno_file
    File ld_scores_file
    File genetic_map_file
    File imputed_samples_file
    File covar_file
    File bgen_file

    String pheno_col
    String qCovarCol #need to figure out the best way to split string so for now the user would have to input things in the format "--qCovarCol=covar1 --qCovarCol=covar2"

    #Array[File] imputed_bgen_files = read_lines(imputed_bgen_filelist)
    #Array[String] bgen_files_with_input_option = prefix("--bgenFile ", imputed_bgen_files) 
    
	Int? memory = 32
	Int? disk = 200
    Int? threads = 32
    Int? preemptible_tries = 3

	command {
        bolt \
            --bed=${genotype_bed} \
            --bim=${genotype_bim} \
            --fam=${genotype_fam} \
            --bgenFile ${bgen_file} \
            --phenoFile=${pheno_file} \
            --phenoCol=${pheno_col} \
            --LDscoresFile=${ld_scores_file} \
            --geneticMapFile=${genetic_map_file} \
            --lmmForceNonInf \
            --LDscoresMatchBp \
            --numThreads=${threads} \
            --covarFile=${covar_file} \
            ${qCovarCol} \
            --statsFile=bolt_genotype_stats_${pheno_col}.gz \
            --sampleFile=${imputed_samples_file} \
            --statsFileBgenSnps=bolt_imputed_stats_${pheno_col}.gz \
            --noBgenIDcheck \
            --verboseStats > ${pheno_col}.boltlmm.log 2>&1
	}

	runtime {
		docker: "quay.io/shukwong/bolt-lmm:2.3.4"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: preemptible_tries
	}

	output {
		File genotype_stats_file = "bolt_genotype_stats_${pheno_col}.gz"
        File imputed_stats_file = "bolt_imputed_stats_${pheno_col}.gz"
        File bolt_lmm_log_file = "${pheno_col}.boltlmm.log"
	}
}
