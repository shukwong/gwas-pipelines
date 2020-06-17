task run_gemma_relatedness_matrix {

        File genotype_bed
        File genotype_bim
        File genotype_fam


        String genotype_prefix = basename(genotype_bed, ".bed")
    
	    Int? memory = 32
	    Int? disk = 200
        Int? threads = 32


	command {
        gemma \
        -bfile ${genotype_prefix} -gk 2 -o  ${genotype_prefix}"_relatedness_matrix.gemma.txt.sXX.txt"
	}

	runtime {
		docker: "quay.io/shukwong/gemma:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

	output {
		File gemma_relatedness_matrix_file="${genotype_prefix}_relatedness_matrix.gemma.txt.sXX.txt"
	}
}

task run_gemma_lmm {

        File imputed_bed
        File imputed_bim
        File imputed_fam
        File gemma_relatedness_matrix_file

        String imputed_prefix = basename (imputed_bed, ".bed")

        Int? memory = 16
	    Int? disk = 200
        Int? threads = 8

    command {
        gemma \
           -bfile ${imputed_prefix} \
           -k ${gemma_relatedness_matrix_file} \
           -lmm 4 \
           -o gemma_lmm4
    }

    runtime {
		docker: "quay.io/shukwong/gemma:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		gpu: false
	}

	output {
		File gemma_association_file="gemma_lmm4.assoc.txt"
        File gemma_log_file="gemma_lmm4.log.txt"
	}
}