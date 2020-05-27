version 1.0
workflow run_preprocess {
    
    input {
	    File genotype_bed
	    File genotype_bim
	    File genotype_fam
	
        File chain_file

        File imputed_files_dir

    }

    call liftover_plink_bim {
        input:
    		genotype_bed = genotype_bed,
    	    genotype_bim = genotype_bim,
    	    genotype_fam = genotype_fam,
            chain_file = chain_file
    }



}


task liftover_plink_bim {
    input {
   	    File genotype_bed
	    File genotype_bim
	    File genotype_fam

        File chain_file

        Int? memory = 32
        Int? disk = 500
    }

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

    meta {
		author: "Wendy Wong"
		email: "wendy.wong@gmail.com"
		description: "Preprocess Genotype files for GWAS"
	}
}

