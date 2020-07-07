
workflow run_saige {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    File bgen_list_file
    File imputed_samples_file
    File covar_file
    #File plink_pca_eigenvec_file

    String phenoCol
    String covarColList #covar list, separated by comma
    #String covar_sampleID_colname
    

    Array[Array[String]] bgen_file_list = read_tsv(bgen_list_file)

    call match_genotype_and_imputed_samples {
        input:
            genotype_bed = genotype_bed,
            genotype_bim = genotype_bim,
            genotype_fam = genotype_fam,
            imputed_samples_file = imputed_samples_file
    }

  
    # call addPCs_to_covar_matrix {
    #     input:
    #         covar_file = covar_file, 
    #         plink_pca_eigenvec_file = plink_pca_eigenvec_file,
    #         phenoCol = phenoCol,
    #         covar_sampleID_colname = covar_sampleID_colname, 
    #         covarColList = covarColList,
    #         sampleFile = imputed_samples_file
    # }
    

    call saige_step1_fitNULL  {
        input:
            genotype_bed = match_genotype_and_imputed_samples.matched_genotype_bed,
	        genotype_bim = match_genotype_and_imputed_samples.matched_genotype_bim,
	        genotype_fam = match_genotype_and_imputed_samples.matched_genotype_fam,
	        covar_file = covar_file, 
            phenoCol = phenoCol,
            covarColList = covarColList
    }

    scatter (bgen_file_line in bgen_file_list) {

	call saige_step2_SPAtests {
	   input:
                bgen_file = bgen_file_line[0],
                bgen_file_index = bgen_file_line[1],
                imputed_samples_file = imputed_samples_file,
                chrom = bgen_file_line[2],
                gmmat_model_file = saige_step1_fitNULL.gmmat_model_file,
                variance_ratio_file = saige_step1_fitNULL.variance_ratio_file,
                sparse_sigma_file = saige_step1_fitNULL.sparse_sigma_file
	    }
    }

    call combine_saige_results {
        input: 
            saige_result_files = saige_step2_SPAtests.saige_output_file,
            pheno_col = phenoCol
    }

    output {
        File merged_saige_file = combine_saige_results.merged_saige_file
	}
    

    meta {
	     author : "Wendy Wong"
	     email : "wendy.wong@gmail.com"
	     description : "Run SAIGE"
    }

}

task addPCs_to_covar_matrix {
    File covar_file 
    File plink_pca_eigenvec_file
    
    String phenoCol
    String covar_sampleID_colname
    String covarColList

    File? sampleFile

    Int? memory = 8
    Int? disk = 20
    Int? threads = 1

    command {
        Rscript combine_covars.R ${covar_file} ${plink_pca_eigenvec_file} ${phenoCol} \
            ${covar_sampleID_colname} ${covarColList} pcs_${covar_file} ${sampleFile}
    }

    runtime {
		docker: "rocker/tidyverse:3.6.3"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

    output {
        File covar_file_with_pcs =  "pcs_${covar_file}"
    }
}

task match_genotype_and_imputed_samples {
    File genotype_bed
    File genotype_bim
    File genotype_fam
    File imputed_samples_file

    Int? memory = 32
    Int? disk = 200
    Int? threads = 32

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


task saige_step1_fitNULL {

    File genotype_bed
    File genotype_bim
    File genotype_fam
    File covar_file #covar_file needs to include covars

    String phenoCol
    String covarColList #covar list, separated by comma
    Float relatedness_cutoff

    Int? memory = 64
    Int? disk = 200
    Int? threads = 32

    command {
        step1_fitNULLGLMM.R     \
            --plinkFile=${sub(genotype_bed,'\\.bed$','')} \
            --phenoFile=${covar_file} \
            --phenoCol=${phenoCol} \
            --covarColList=${covarColList} \
            --sampleIDColinphenoFile=IID \
            --traitType=binary        \
            --IsSparseKin=TRUE \
            --relatednessCutoff=${relatedness_cutoff} \
            --outputPrefix=saige_step1_${phenoCol} \
            --nThreads=${threads}
	
	spaseSigma_file=( "*.sparseSigma.mtx" )
	mv $spaseSigma_file saige_step1_${phenoCol}.sparseSigma.mtx
	}

    runtime {
		docker: "wzhou88/saige:0.38"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

    output {
	File gmmat_model_file = "saige_step1_${phenoCol}.rda"
        File variance_ratio_file = "saige_step1_${phenoCol}.varianceRatio.txt"
        File sparse_sigma_file = "saige_step1_${phenoCol}.sparseSigma.mtx"
    }
}


task saige_step2_SPAtests {

    File bgen_file
    File bgen_file_index
    File imputed_samples_file

    File gmmat_model_file
    File variance_ratio_file
    File sparse_sigma_file

    String chrom 

    String file_prefix = basename(bgen_file, ".bgen") 
    String saige_output_file_name = file_prefix + "." + chrom + ".txt"

	Int? memory = 64
	Int? disk = 500
    Int? threads = 64

	command {
      step2_SPAtests.R \
        --bgenFile=${bgen_file} \
        --bgenFileIndex=${bgen_file_index} \
        --IsDropMissingDosages=FALSE \
        --minMAF=0.01 \
        --minMAC=1 \
        --chrom=${chrom} \
        --sampleFile=${imputed_samples_file} \
        --GMMATmodelFile=${gmmat_model_file} \
        --varianceRatioFile=${variance_ratio_file} \
        --sparseSigmaFile=${sparse_sigma_file} \
        --SAIGEOutputFile=${saige_output_file_name} \
        --numLinesOutput=2 \
        --IsOutputNinCaseCtrl=TRUE \
        --IsOutputHetHomCountsinCaseCtrl=TRUE \
        --IsOutputAFinCaseCtrl=TRUE

        gzip ${saige_output_file_name}
	}

	runtime {
		docker: "wzhou88/saige:0.38"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: 64
		gpu: false
	}

	output {
		File saige_output_file = "${saige_output_file_name}.gz"
	}
}


task combine_saige_results {

    Array[File] saige_result_files
    String pheno_col

    Int? memory = 4
	Int? disk = 100
    Int? threads = 1
    Int? preemptible_tries = 3

    command {
        echo -e "CHR\tPOS\trsid\tSNPID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\timputationInfo\tN\tBETA\tSE\tTstat\tp.value\tp.value.NA\tIs.SPA.converge\tvarT\tvarTstar\tAF.Cases\tAF.Co
ntrols\tN.Cases\tN.Controls\thomN_Allele2_cases\thetN_Allele2_cases\thomN_Allele2_ctrls\thetN_Allele2_ctrls" > saige_${pheno_col}_results_merged.tsv
        
        cat ${sep=' ' saige_result_files} | gzip -d | grep -v ^CHR >> saige_${pheno_col}_results_merged.tsv
        
        gzip saige_${pheno_col}_results_merged.tsv
    }

    runtime {
		docker: "alpine"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: preemptible_tries
	}

    output {
        File merged_saige_file = "saige_${pheno_col}_results_merged.tsv.gz"
    }
}