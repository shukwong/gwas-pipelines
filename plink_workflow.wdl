
workflow run_preprocess {
    
	File genotype_bed
	File genotype_bim
	File genotype_fam
	
    File chain_file

    File imputed_files_dir

	Int? memory = 60
	Int? disk = 500
	Int? threads = 16


	call liftover_plink_bim {
        input:
    		genotype_bed = genotype_bed,
    	    genotype_bim = genotype_bim,
    	    genotype_fam = genotype_fam,
            chain_file = chain_file
    }

    call subset_plink_and_update_bim {
        input:
            genotype_bed = genotype_bed,
            genotype_bim = genotype_bim,
            genotype_fam = genotype_fam,
            mapped_ids = liftover_plink_bim.mapped_ids,
            mapped_bim = liftover_plink_bim.mapped_bim
    }

 	call run_ld_prune {
        input:
            genotype_bed = subset_plink_and_update_bim.output_bed,
            genotype_bim = subset_plink_and_update_bim.output_bim,
            genotype_fam = subset_plink_and_update_bim.output_fam
    }

	call plink_pca {
        input:
    		genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam
            
    }

	output {
        File genotype_pruned_bed = run_ld_prune.genotype_pruned_bed
        File genotype_pruned_bim = run_ld_prune.genotype_pruned_bim
        File genotype_pruned_fam = run_ld_prune.genotype_pruned_fam
		File genotype_pruned_pca_eigenvec = plink_pca.genotype_pruned_pca_eigenvec
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
    Int? disk = 500

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

task run_ld_prune {
    
    File genotype_bed
    File genotype_bim
    File genotype_fam

    Int? memory = 32
    Int? disk = 500
    

    command <<<
		
        plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam}  --indep 50 5 2 --out ld_indep_check

        plink --keep-allele-order --bed ${genotype_bed} --bim ${genotype_bim} \
              --fam ${genotype_fam} --extract ld_indep_check.prune.in \
              --maf 0.01 --make-bed --out ld_indep_check.prune

        plink --bfile /mnt/data/munge/ld_indep_check.prune  --indep-pairwise 50 5 0.5 \
              --make-bed --out ld_indep_pairwise_check

        plink --keep-allele-order --bfile /mnt/data/munge/ld_indep_check.prune \
              --extract ld_indep_pairwise_check.prune.in \
              --make-bed --out genotype_pruned_plink
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
    }
}

task plink_to_vcf {

	File genotype_bed
	File genotype_bim
	File genotype_fam
	String prefix = basename(genotype_bed, ".bed")
	Int? memory = 32
	Int? disk = 500

	command {
        plink --bfile ${prefix} --recode vcf --out ${prefix}.vcf   
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
	Int? disk = 500

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
    String prefix = basename(vcf_file, ".vcf.gz")
    Int? bits=8

    Int? memory = 60
    Int? disk = 500

	command <<<
        /plink2 --vcf ${vcf_file} \
            --make-pgen erase-phase --out plink_out

        /plink2 --pfile plink_out --export bgen-1.2 bits=${bits} --out ${prefix}
	>>>

	runtime {
		docker: "quay.io/large-scale-gxe-methods/genotype-conversion"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

	output {
		File out_bgen = "${prefix}.bgen"
		File out_bgen_sample = "${prefix}.sample"
		File out_bgen_log = "${prefix}.log"
	}
}

task liftover_plink_bim {
   	File genotype_bed
	File genotype_bim
	File genotype_fam

    File chain_file

    Int? memory = 32
    Int? disk = 500

    command <<<
		
        awk 'start=$4-1 {print "chr"$1"\t"start"\t"$4"\t"$2"\t"$4"\t"$5"\t"$6}' ${genotype_bim} >bim_as_bed.bed

        CrossMap.py bed ${chain_file} bim_as_bed.bed bim_as_bed.crossmap.bed 

        cat bim_as_bed.crossmap.bed  | grep -v ^chrUn | grep -v random \
            | grep -v alt | grep -v ^chrX | grep -v ^chrY | sort -k1,1 -k2,2g \
            | cut -f4 >bim_as_bed.mapped.ids

        cat bim_as_bed.crossmap.bed | grep -v ^chrUn | grep -v random | grep -v alt \
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
	    File mapped_ids = "bim_as_bed.mapped.ids"
        File mapped_bim = "bim_as_bed.mapped.bim"
    }
}

task subset_plink_and_update_bim {

    File genotype_bed
    File genotype_bim
    File genotype_fam
    File mapped_ids
    File mapped_bim
   
    Int? memory = 32
    Int? disk = 500
    
 
    command <<<

        plink \
            --bed ${genotype_bed} --bim ${genotype_bim} \
            --fam ${genotype_fam} --extract ${mapped_ids} \
            --make-bed --out genotypes_updated
            
        cp ${mapped_bim} genotypes_updated.bim   

    >>>    

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

task plink_bed_subset_sample {

    File genotype_bed
    File genotype_bim
    File genotype_fam
    File samples_to_keep_file

    String plink_bed_prefix

    Int? memory = 32
    Int? disk = 500
	
    command {
		
	    plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} --keep ${samples_to_keep_file} --make-bed --out ${plink_bed_prefix}

    }


	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	    File plink_bed = "${plink_bed_prefix}.bed"
        File plink_bim = "${plink_bed_prefix}.bim"
        File plink_fam = "${plink_bed_prefix}.fam"
    }
}
