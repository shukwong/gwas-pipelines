
workflow bolt_workflow {
    
	File genotype_bed
	File genotype_bim
	File genotype_fam
    File pheno_file
    File ld_scores_file
    File genetic_map_file
    #File imputed_samples_file
    File covar_file
    File bgen_list_file

    String pheno_col
    String qCovarCol #need to figure out the best way to split string so for now the user would have to input things in the format "--qCovarCol=covar1 --qCovarCol=covar2"

    #preprocess 
    # call match_genotype_and_imputed_samples {
    #       input:
    #          genotype_bed = genotype_bed,
    #          genotype_bim = genotype_bim,
    #          genotype_fam = genotype_fam,
    #          imputed_samples_file = imputed_samples_file
    # }

    #Array[Array[File]] bgen_files_and_indices
    Array[Array[File]] bgen_files_and_indices = read_tsv(bgen_list_file)

    scatter (bgen_file_line in bgen_files_and_indices) {

        # call subset_bgen_from_genotype {
        #     input:
        #         plink_fam_file = genotype_fam,
        #         bgen_file = bgen_file_line[0],
        #         imputed_samples_file = imputed_samples_file

        # }

        call bgen_get_samples_file {
            input: bgen_file = bgen_file_line[0]
        }

    
        call run_bolt_lmm {
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
            bgen_samples_file = bgen_get_samples_file.bgen_samples_file,
            covar_file = covar_file,
            bgen_file = bgen_file_line[0],
            pheno_col = pheno_col,
            #qCovars = qCovars
            qCovarCol = qCovarCol
        }
    }

    call combine_bolt_results { 
        input: 
           imputed_stats_files = run_bolt_lmm.imputed_stats_file, 
           pheno_col = pheno_col
    }

    output {
		#File genotype_stats_file = run_bolt.genotype_stats_file
        File imputed_stats_file = combine_bolt_results.merged_stats_file
        #Array[File] bolt_lmm_log_files = run_bolt_lmm.bolt_lmm_log_file
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
        #echo -e "SNP\tCHR\tBP\tGENPOS\tALLELE1\tALLELE0\tA1FREQ\tINFO\tCHISQ_LINREG\tP_LINREG\tBETA\tSE\tCHISQ_BOLT_LMM_INF\tP_BOLT_LMM_INF\tCHISQ_BOLT_LMM\tP_BOLT_LMM" > bolt_${pheno_col}_results_merged.tsv
        
        echo -e "CHR\tPOS\tSNP\tTested_Allele\tOther_Allele\tBETA\tSE\tP" > bolt_${pheno_col}_results_merged.tsv

        cat ${sep=' ' imputed_stats_files} | gzip -d | grep -v ^SNP | awk '{print $2"\t"$3"\t"$1"\t"$5"\t"$6"\t"$11"\t"$12"\t"$14}' >> bolt_${pheno_col}_results_merged.tsv
        
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

task run_bolt_lmm {

	File genotype_bed
	File genotype_bim
	File genotype_fam
    File pheno_file
    File ld_scores_file
    File genetic_map_file
    File bgen_samples_file
    File covar_file
    File bgen_file

    String pheno_col
    #Array[String] qCovars
    String qCovarCol 
    #Array[File] imputed_bgen_files = read_lines(imputed_bgen_filelist)
    #Array[String] bgen_files_with_input_option = prefix("--bgenFile ", imputed_bgen_files) 
    
	Int? memory = 32
	Int? disk = 200
    Int? threads = 32
    Int? preemptible_tries = 3

    #need to figure out the best way to split string so for now the user would have to input things in the format "--qCovarCol=covar1 --qCovarCol=covar2"
    #String qCovarOption = " --qCovarCol=" + qCovarCol + " "

	command <<<
        qCovars=$(echo ${qCovarCol} | awk '{print ","$0}' |  sed 's/,/ --qCovarCol=/g')

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
            $qCovars \
            --statsFile=bolt_genotype_stats_${pheno_col}.gz \
            --sampleFile=${bgen_samples_file} \
            --statsFileBgenSnps=bolt_imputed_stats_${pheno_col}.gz \
            --noBgenIDcheck \
            --verboseStats 
	>>>

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
	}
}

task subset_bgen_from_genotype {
    File plink_fam_file 
    File bgen_file
    File imputed_samples_file

    Int? memory = 16
    Int? disk = 200
    Int? threads = 4

    command <<<
        echo -e "ID_1 ID_2\n0 0" >imputed.bgen.samples
        awk '{print $1" "$1}' ${imputed_samples_file} >>imputed.bgen.samples       

        echo -e "ID_1 ID_2\n0 0" >plink.samples
        awk '{print $1" "$2}' ${plink_fam_file} >>plink.samples    

        qctool -g ${bgen_file} -s imputed.bgen.samples -incl-samples plink.samples  -bgen-bits 8 -og subsetted.bgen

        qctool -g subsetted.bgen -os subsetted_bgen.samples

    >>>

    runtime {
		docker: "quay.io/shukwong/qctool:v2.0.8"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

    output {
        File subsetted_bgen = "subsetted.bgen"
        File subsetted_bgen_samples = "subsetted_bgen.samples"
    }

}


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


task bgen_get_samples_file {
    File bgen_file

    Int? memory = 16
    Int? disk = 200
    Int? threads = 4

    command <<<
       
        qctool -g ${bgen_file} -os bgen.samples

    >>>

    runtime {
		docker: "quay.io/shukwong/qctool:v2.0.8"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

    output {
        File bgen_samples_file = "bgen.samples"
    }

}