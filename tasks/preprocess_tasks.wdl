
task plink_pca {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    String? approx = "approx"

    Int? memory = 60
    Int? disk = 200

    command {
		/plink2 --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} --pca ${approx} --out genotype_pruned_pca
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/genotype-conversion"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	File genotype_pruned_pca_eigenvec = "genotype_pruned_pca.eigenvec"
	File genotype_pruned_pca_eigenval = "genotype_pruned_pca.eigenval"
        File genotype_pruned_pca_log = "genotype_pruned_pca.log"
    }
}

task run_genotype_qc_filter {
	
    File genotype_bed
    File genotype_bim
    File genotype_fam

    Int? memory = 60
    Int? disk = 200

    command <<<

	/plink2 --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} \
              --maf 0.01 --geno 0.02 --hwe 0.001 --make-bed --out gentoype_qc_filtered

    >>>

    runtime {
	docker: "quay.io/large-scale-gxe-methods/genotype-conversion"
        memory: "${memory} GB"
        disks: "local-disk ${disk} HDD"
        gpu: false
    }
	
    output {
	File gentoype_qc_filtered_bed = "gentoype_qc_filtered.bed"	
        File gentoype_qc_filtered_bim = "gentoype_qc_filtered.bim"
        File gentoype_qc_filtered_fam = "gentoype_qc_filtered.fam"
    } 	
}


task run_ld_prune {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    Int? memory = 32
    Int? disk = 200
    

    command <<<
	
        plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam}  --indep 50 5 2 --out ld_indep_check

        plink --keep-allele-order --bed ${genotype_bed} --bim ${genotype_bim} \
              --fam ${genotype_fam} --extract ld_indep_check.prune.in \
              --make-bed --out ld_indep_check.prune

        plink --bfile ld_indep_check.prune  --indep-pairwise 50 5 0.2 \
              --out ld_indep_pairwise_check

        plink --keep-allele-order --bfile ld_indep_check.prune \
              --extract ld_indep_pairwise_check.prune.in \
              --make-bed --out genotype_pruned_plink

        plink --bfile genotype_pruned_plink --het --out genotype_pruned_plink_het
        awk 'NR > 1 && sqrt($6^2) > sqrt(0.2^2) {print $1"\t"$1}' genotype_pruned_plink_het.het > genotype_pruned_plink_het.het.remove      
        ## make a version of the data with QCed autosomal data for later
        plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} \
              --remove genotype_pruned_plink_het.het.remove --autosome \
              --make-bed --out genotype.nohet.autosomes
    >>>

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	File genotype_pruned_bed = "genotype_pruned_plink.bed"
        File genotype_pruned_bim = "genotype_pruned_plink.bim"
        File genotype_pruned_fam = "genotype_pruned_plink.fam"
        File genotype_nohet_autosomes_bed = "genotype.nohet.autosomes.bed"
        File genotype_nohet_autosomes_bim = "genotype.nohet.autosomes.bim"
        File genotype_nohet_autosomes_fam = "genotype.nohet.autosomes.fam"
    }
}

task plink_to_vcf {

	File genotype_bed
	File genotype_bim
	File genotype_fam

        String prefix = basename (genotype_bed, ".bed")
	
	Int? memory = 32
	Int? disk = 200

	command {
        plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} --recode vcf --out ${prefix}.vcf   
	}

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

	output {
		File out_vcf = "${prefix}.vcf"
	}
}

task vcf_to_plink_bed {

	File vcf_file
        String prefix = basename(vcf_file, ".vcf")
	Int? memory = 32
	Int? disk = 200

	command {
		plink --vcf ${vcf_file}  --make-bed --out ${prefix}
	}

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

	output {
		File out_bed = "${prefix}.bed"
		File out_bim = "${prefix}.bim"
		File out_fam = "${prefix}.fam"
	}
}

task vcf_to_bgen {
    File vcf_file
    File samples_to_keep_file

    Int? bits=8
    Int? memory = 24
    Int? memory_in_MB = 24000
    Int? disk = 100

    String prefix = basename(vcf_file, ".vcf.gz")
    

	command <<<
        /plink2 --memory ${memory_in_MB} --vcf ${vcf_file} dosage=HDS --id-delim _ \
            --make-pgen erase-phase --out plink_out

        /plink2 --memory ${memory_in_MB} --pfile plink_out --keep ${samples_to_keep_file} \
            --export bgen-1.2 bits=${bits} --out ${prefix}

        ## mysterious file format discrepancy causes failures with FASTGWA; change NA values to 0 in sample files
        sed 's/NA$/0/' ${prefix}.sample > ${prefix}-noNAs.sample    

        awk 'NR > 2 {print $1}' ${prefix}-noNAs.sample  > ${prefix}_bgen._saige.sample #samples file for saige
	>>>

	runtime {
		docker: "quay.io/large-scale-gxe-methods/genotype-conversion"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

	output {
		File out_bgen = "${prefix}.bgen"
		File out_bgen_sample = "${prefix}-noNAs.sample"
		File out_bgen_log = "${prefix}.log"
        File out_bgen_saige_sample = "${prefix}_bgen._saige.sample"
	}
}

task index_bgen_file {
    File bgen_file

    String bgen_filename = basename(bgen_file)	

    Int? memory = 60
    Int? disk = 100

    command {
        bgenix -g ${bgen_file} -index -clobber
	mv ${bgen_file}.bgi ./
    }

    output {
        File bgen_file_index = "${bgen_filename}.bgi"
    }

    runtime {
		docker: "quay.io/shukwong/bgen:fe6d17aa6933"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}
}


task liftover_plink {
   	File genotype_bed
	File genotype_bim
	File genotype_fam

    File chain_file

    Int? memory = 32
    Int? disk = 200

    command <<<
		
        awk 'start=$4-1 {print "chr"$1"\t"start"\t"$4"\t"$2"\t"$4"\t"$5"\t"$6}' ${genotype_bim} >bim_as_bed.bed

        CrossMap.py bed ${chain_file} bim_as_bed.bed bim_as_bed.crossmap.bed 

        grep -v ^chrUn bim_as_bed.crossmap.bed | grep -v random \
            | grep -v alt | grep -v ^chrX | grep -v ^chrY | sort -k1,1 -k2,2g \
            | cut -f4 >bim_as_bed.mapped.ids

        plink \
            --bed ${genotype_bed} --bim ${genotype_bim} \
            --fam ${genotype_fam} --extract bim_as_bed.mapped.ids \
            --make-bed --out genotypes_updated    
    >>>        

	runtime {
		docker: "quay.io/shukwong/plinkcrossmap:v1"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	    File output_bed = "genotypes_updated.bed"
        File output_bim = "genotypes_updated.bim"
        File output_fam = "genotypes_updated.fam"
    }
}



task plink_subset_sample {

    File genotype_bed
    File genotype_bim
    File genotype_fam
    File samples_to_keep_file

    Int? memory = 32
    Int? disk = 500
	
    command {
		plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} --keep ${samples_to_keep_file} --make-bed --out genotype_subsetSample
    }

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	File genotype_subsetSample_bed = "genotype_subsetSample.bed"
        File genotype_subsetSample_bim = "genotype_subsetSample.bim"
        File genotype_subsetSample_fam = "genotype_subsetSample.fam"
    }
}

task add_pcs_to_covar_file {
    File eigenvec_file
    File covar_file
    File samples_to_keep_file
   
    String  phenotype 

    Int? memory = 32
    Int? disk = 20

    command {
    
        wget https://github.com/shukwong/gwas-pipelines/raw/master/scripts/construct_model_matrix.R
    
        Rscript construct_model_matrix.R ${covar_file} ${samples_to_keep_file} ${phenotype} ${phenotype}_model_matrix.tsv
    }

    runtime {
		docker: "rocker/tidyverse:3.6.3-ubuntu18.04"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	    File out_covar_file = "${phenotype}_model_matrix.tsv"
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

task convert_gen_to_bgen {
    Array[File] gen_files

    String chrom
    
    Float threshold

    Boolean? outputBgenOnePointTwo = true

    Int? memory = 32
    Int? disk = 200
    Int? threads = 32
    Int? preemptible_tries = 3

    command <<<
        cat ${sep=' ' gen_files} > ${chrom}.gen

        qctool -g ${chrom}.gen  -threshold ${threshold}  -filetype gen -ofiletype \
             ${true='bgen_v1.2 -bgen-bits 8' false='bgen' outputBgenOnePointTwo}  -og ${chrom}.bgen 
    >>>

    runtime {
		docker: "quay.io/shukwong/qctool:v2.0.8"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: "${preemptible_tries}"
	}

    output {
        File converted_bgen_file = "${chrom}.bgen"
    }

}


task get_covar_subsets {
    File covariate_tsv_file 
    File variable_info_tsv_file 
    File sample_sets_json_file

    Int? memory = 4
    Int? disk = 200
    Int? threads = 1
    Int? preemptible_tries = 3

#TODO, change this to git clone a release version when the pipeline is finalized
    command <<< 
       
        wget https://github.com/shukwong/gwas-pipelines/raw/master/scripts/create_covar_files_by_set.R

        Rscript create_covar_files_by_set.R ${covariate_tsv_file} ${variable_info_tsv_file} ${sample_sets_json_file}
    >>>

    runtime {
		docker: "rocker/tidyverse:4.0.0"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: "${preemptible_tries}"
	}

    output {
        Array[File] covar_subsets_files = glob("covars*.tsv")
        Array[File] covar_subsets_log_files = glob("*.log")
    }

}

task get_cohort_samples {
    File covariate_tsv_file
    File genotype_samples_to_keep_file
    File imputed_samples_to_keep_file

    Int? memory = 16
    Int? disk = 200
    Int? threads = 1
    Int? preemptible_tries = 3

    command <<<

        wget https://github.com/shukwong/gwas-pipelines/raw/master/scripts/get_cohort_samples.R
        Rscript get_cohort_samples.R ${covariate_tsv_file} ${genotype_samples_to_keep_file} ${imputed_samples_to_keep_file} 
    >>>

    runtime {
		docker: "rocker/tidyverse:4.0.0"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: "${preemptible_tries}"
	}

    output {
        File covar_subset_file = "covars_subsetted.tsv"
        File plink_subset_samples= "plink_subsetted.samples"
    }
}