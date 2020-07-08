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
    call match_genotype_and_imputed_samples {
          input:
             genotype_bed = genotype_bed,
             genotype_bim = genotype_bim,
             genotype_fam = genotype_fam,
             imputed_samples_file = imputed_samples_file
    }

    Array[Array[String]] bgen_file_list = read_tsv(bgen_list_file)

    scatter (bgen_file_line in bgen_file_list) {
    
        call run_bolt {
            input:
            genotype_bed = match_genotype_and_imputed_samples.matched_genotype_bed,
	        genotype_bim = match_genotype_and_imputed_samples.matched_genotype_bim,
	        genotype_fam = match_genotype_and_imputed_samples.matched_genotype_fam,
            # genotype_bed = genotype_bed,
            # genotype_bim = genotype_bim,
            # genotype_fam = genotype_fam,
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

    call combine_bolt_results { 
        input: 
           imputed_stats_files = run_bolt.imputed_stats_file, 
           pheno_col = pheno_col
    }

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
        echo -e "SNP\tCHR\tBP\tGENPOS\tALLELE1\tALLELE0\tA1FREQ\tINFO\tCHISQ_LINREG\tP_LINREG\tBETA\tSE\tCHISQ_BOLT_LMM_INF\tP_BOLT_LMM_INF\tCHISQ_BOLT_LMM\tP_BOLT_LMM" > bolt_${pheno_col}_results_merged.tsv
        
        cat ${sep=' ' imputed_stats_files} | gzip -d | grep -v ^SNP >> bolt_${pheno_col}_results_merged.tsv
        
        gzip bolt_${pheno_col}_results_merged.tsv
    }

    runtime {
		docker: "ubuntu:20.10"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: preemptible_tries
	}

    output {
        File merged_stats_file = "bolt_${pheno_col}_results_merged.tsv.gz"
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

#TODO, generate PCs and add it to the covar file
task match_genotype_and_imputed_samples {
    File genotype_bed
    File genotype_bim
    File genotype_fam
    File imputed_samples_file

    Int? memory = 32
    Int? disk = 200
    Int? threads = 8

    command <<<
        awk '{print $1"\t"$1}' ${imputed_samples_file}  > samples_plink_format.txt

        /plink2 --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} --keep samples_plink_format.txt \
            --make-bed --out matched_genotype
    >>>

    runtime {
		docker: "quay.io/large-scale-gxe-methods/genotype-conversion"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

    output {
        File matched_genotype_bed = "matched_genotype.bed"
        File matched_genotype_bim = "matched_genotype.bim"
        File matched_genotype_fam = "matched_genotype.fam"
    }

}