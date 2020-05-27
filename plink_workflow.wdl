
workflow run_preprocess {
    
	File genotype_bed
	File genotype_bim
	File genotype_fam
	
    File chain_file

    File imputed_files_dir

	Int? memory = 60
	Int? disk = 500
	Int? threads = 16


    call plink_pca {
 
        input:
    		genotype_bed = genotype_bed,
    	    genotype_bim = genotype_bim,
    	    genotype_fam = genotype_fam
            
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
