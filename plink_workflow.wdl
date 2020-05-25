#some tasks are modified from https://github.com/large-scale-gxe-methods/genotype-conversion/blob/master/genotype_conversion.wdl

task run_ld_prune {
  
	File genotype_bed
	File genotype_bim
	File genotype_fam

    String genotype_pruned_plink

	File genotype_pruned_pca

    command {
		
        plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam}  --indep-pairwise 100 20 0.2 --out check
	    plink --keep-allele-order --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} --extract check.prune.in --make-bed --out ${genotype_pruned_plink}
	    plink --threads ${threads} --bfile ${genotype_pruned_plink} --pca --out ${genotype_pruned_pca}
   	}


	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	    File genotype_pruned_pca = "genotype_pruned_pca"
    }
}

task plink_bed_subset_sample {
  
	File genotype_bed
	File genotype_bim
	File genotype_fam
    File samples_to_keep_file

    String plink_bed_prefix

	
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
	    File plink_bed = ${plink_bed_prefix}.bed
        File plink_bim = ${plink_bed_prefix}.bim
        File plink_fam = ${plink_bed_prefix}.fam
    }
}

task run_crossmap {
    File input_file
	File chain_file
    String type
    String prefix = basename(vcf_file, ".${type}")

    command {
		
       CrossMap.py ${type} ${chain_file} ${input_file} ${prefix}.${type}
   	}    

	runtime {
		docker: "crukcibioinformatics/crossmap"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	    File output = ${prefix}.${type}
    }
}

task vcf_to_plink_bed {

	File vcf_file
    String prefix = basename(vcf_file, ".vcf")
	Int? memory = 10
	Int? disk = 20

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

task plink_to_vcf {

	File genotype_bed
	File genotype_bim
	File genotype_fam
	String prefix = basename(genotype_bed, ".bed")
	Int? memory = 16
	Int? disk = 100

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

workflow run_plink {

	File genotype_bed
	File genotype_bim
	File genotype_fam
	
    String genotype_pruned_plink

	Int? memory = 16
	Int? disk = 100
	Int? threads = 8
 
 	call run_ld_prune {
 		input:
 		    genotype_bed = genotype_bed
	        genotype_bim = genotype_bim
	        genotype_fam = genotype_fam

            genotype_pruned_plink = genotype_pruned_plink
 	}

        output {
		File results = cat_results.all_results
 		Array[File] system_resource_usage = run_interaction.system_resource_usage
 		Array[File] process_resource_usage = run_interaction.process_resource_usage
 	}

	parameter_meta {
		genofiles_pgen: "Array of PLINK2 genotype (.pgen) filepaths."
		genofiles_psam: "Array of PLINK2 sample (.psam) filepaths."
		genofiles_pvar: "Array of PLINK2 variant (.pvar) filepaths."
		genotype_pruned_pca: "PCA of the pruned Genotype file"
		cpu: "Minimum number of requested cores."
		disk: "Requested disk space (in GB)."
	}

	meta {
		author: "Wendy Wong"
		email: "wendy.wong@gmail.com"
		description: "Preprocess Genotype files for GWAS"
	}
}