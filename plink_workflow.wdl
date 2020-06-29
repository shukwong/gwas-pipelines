
workflow run_preprocess {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    File gsa_samples_to_keep_file
    File imputed_samples_to_keep_file
    File imputed_list_of_vcfs_file
    
    File chain_file

    call plink_subset_sample {
        input:
             genotype_bed = genotype_bed, 
             genotype_bim = genotype_bim,
             genotype_fam = genotype_fam,
             samples_to_keep_file = gsa_samples_to_keep_file
    }

    call run_genotype_qc_filter {
        input:
            genotype_bed = plink_subset_sample.genotype_subsetSample_bed,
            genotype_bim = plink_subset_sample.genotype_subsetSample_bim,
            genotype_fam = plink_subset_sample.genotype_subsetSample_fam

    }

    call run_ld_prune {
        input:
            genotype_bed = run_genotype_qc_filter.gentoype_qc_filtered_bed,
            genotype_bim = run_genotype_qc_filter.gentoype_qc_filtered_bim,
            genotype_fam = run_genotype_qc_filter.gentoype_qc_filtered_fam
    }

    call plink_pca {
        input:
    	    genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam
            
    }

    call liftover_plink_bim {
        input:
    	    genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam,
            chain_file = chain_file
    }

    call liftover_plink {
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
	call vcf_to_bgen {
	    input:
                vcf_file = imputed_file[0],
                samples_to_keep_file = imputed_samples_to_keep_file
		}
        call index_bgen_file {
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
		description : "Preprocess Genotype files for GWAS"
	}

}

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
    File liftover_mapped_ids_file
    File liftover_mapped_new_bim_file

    Int? memory = 32
    Int? disk = 200

    command {
        plink \
            --bed ${genotype_bed} --bim ${genotype_bim} \
            --fam ${genotype_fam} --extract ${liftover_mapped_ids_file} \
            --make-bed --out genotypes_updated
            
        cp ${liftover_mapped_new_bim_file} genotypes_updated.bim       
    }

    runtime {
		docker: "quay.io/h3abionet_org/py3plink"
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

task liftover_plink_bim {
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

        grep -v ^chrUn bim_as_bed.crossmap.bed | grep -v random | grep -v alt \
            | grep -v ^chrX | grep -v ^chrY | sort -k1,1 -k2,2g | sed 's/^chr//' \
            | awk '{print $1"\tchr"$1":"$3":"$6":"$7"\t"$3"\t"$6"\t"$7}' > \
            bim_as_bed.mapped.bim
    >>>        

	runtime {
		docker: "crukcibioinformatics/crossmap"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	    File liftover_mapped_ids_file = "bim_as_bed.mapped.ids"
        File liftover_mapped_new_bim_file = "bim_as_bed.mapped.bim"
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
