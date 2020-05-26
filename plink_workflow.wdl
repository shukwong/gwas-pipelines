version 1.0
#some tasks are modified from https://github.com/large-scale-gxe-methods/genotype-conversion/blob/master/genotype_conversion.wdl

task run_ld_prune {
  
	File genotype_bed
	File genotype_bim
	File genotype_fam

    String genotype_pruned_plink

	File genotype_pruned_pca

    command<<<
		
        plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam}  --indep-pairwise 100 20 0.2 --out check
	    plink --keep-allele-order --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} --extract check.prune.in --make-bed --out ${genotype_pruned_plink}
	    plink --threads ${threads} --bfile ${genotype_pruned_plink} --pca --out ${genotype_pruned_pca}


        plink --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam}  --indep 50 5 2 --out ld_indep_check

        plink --keep-allele-order --bed ${genotype_bed} --bim ${genotype_bim} \
              --fam ${genotype_fam} --extract ld_indep_check.prune.in \
              --make-bed --out ld_indep_check.prune

        plink --bfile /mnt/data/munge/ld_indep_check.prune  --indep-pairwise 50 5 0.5 \
              --make-bed --out ld_indep_pairwise_check

        plink --keep-allele-order --bfile /mnt/data/munge/ld_indep_check.prune \
              --extract ld_indep_pairwise_check.prune.in \
              --make-bed --out ${genotype_pruned_plink}
   	>>>


	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	    File genotype_pruned_bed = "${genotype_pruned_plink}.bed"
        File genotype_pruned_bed = "${genotype_pruned_plink}.bed"
        File genotype_pruned_bed = "${genotype_pruned_plink}.bed"

    }
}

task plink_pca {
    File genotype_bed
	File genotype_bim
	File genotype_fam

    String? approx="approx"

    String genotype_pruned_pca


    command {
		/plink2 --bed ${genotype_bed} --bim ${genotype_bim} --fam ${genotype_fam} --pca ${approx} --out ${genotype_pruned_pca}
   	}


	runtime {
		docker: "quay.io/large-scale-gxe-methods/genotype-conversion"
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
	    File plink_bed = "${plink_bed_prefix}.bed"
        File plink_bim = "${plink_bed_prefix}.bim"
        File plink_fam = "${plink_bed_prefix}.fam"
    }
}

task subset_plink_and_update_bim {
    File genotype_bed
	File genotype_bim
	File genotype_fam
    File mapped_ids
	File mapped_bim
    
    command {
		
        plink \
            --bed ${genotype_bed} --bim ${genotype_bim} \
            --fam ${genotype_fam} --extract ${mapped_ids} \
            --make-bed --out genotypes_updated

        cp ${mapped_bim} genotypes_updated.bim   
   	}    

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
		gpu: false
	}

    output {
	    File output_bed = genotypes_updated.bed
        File output_bim = genotypes_updated.bim
        File output_fam = genotypes_updated.fam
    }
}

task liftover_plink_bim {
   	File genotype_bed
	File genotype_bim
	File genotype_fam

    File chain_file


    command {
		
        awk 'start=$4-1 {print "chr"$1"\t"start"\t"$4"\t"$2"\t"$4"\t"$5"\t"$6}' ${genotype_bim} >bim_as_bed.bed

        CrossMap.py bed ${chain_file} bim_as_bed.bed bim_as_bed.crossmap.bed 

        cat bim_as_bed.crossmap.bed  | grep -v ^chrUn | grep -v random \
            | grep -v alt | grep -v ^chrX | grep -v ^chrY | sort -k1,1 -k2,2g \
            | cut -f4 >bim_as_bed.mapped.ids

        cat bim_as_bed.crossmap.bed | grep -v ^chrUn | grep -v random | grep -v alt \
            | grep -v ^chrX | grep -v ^chrY | sort -k1,1 -k2,2g | sed 's/^chr//' \
            | awk '{print $1"\tchr"$1":"$3":"$6":"$7"\t"$3"\t"$6"\t"$7}' > \
            bim_as_bed.mapped.bim

   	}    

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

task vcf_to_bgen {

	File vcf_file
    String prefix = basename(vcf_file, ".vcf.gz")
	Int? memory = 10
	Int? disk = 20
    Int? bits=8

	command {
        /plink2 --vcf ${vcf_file} \
            --make-pgen erase-phase --out plink_out

        /plink2 --pfile plink_out --export bgen-1.2 bits=${bits} --out ${prefix}
	}

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

workflow run_preprocess {

	File genotype_bed
	File genotype_bim
	File genotype_fam
	
    String genotype_pruned_plink

	Int? memory = 60
	Int? disk = 100
	Int? threads = 16
 
    call liftover_plink_bim {
        input:
 		    genotype_bed = genotype_bed,
	        genotype_bim = genotype_bim,
	        genotype_fam = genotype_fam,

            chain_file = chain_file
    }

    #call subset_plink_and_update_bim{}

 	#call run_ld_prune {}

    output {
		File mapped_ids = "liftover_plink_bim.mapped_ids"
        File mapped_bim = "liftover_plink_bim.mapped_bim"
 	}

	parameter_meta {
		genofiles_bed: "PLINK genotype filepath"
		genofiles_bim: "PLINK genotype filepath"
		genofiles_fam: "PLINK genotype filepath"
		cpu: "Minimum number of requested cores."
		disk: "Requested disk space (in GB)."
	}

	meta {
		author: "Wendy Wong"
		email: "wendy.wong@gmail.com"
		description: "Preprocess Genotype files for GWAS"
	}
}