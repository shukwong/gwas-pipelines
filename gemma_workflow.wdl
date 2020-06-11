task run_gemma_relatedness_matrix {
    input:
        File genotype_bed
        File genotype_bim
        File genotype_fam


        String genotype_prefix 
    
	    Int? memory = 32
	    Int? disk = 200
        Int? threads = 32

	command {
        gemma \
        -bfile ${genotype_prefix} -gk 2 -o ${gemma_relatedness_matrix}
	}

	runtime {
		docker: "quay.io/shukwong/gemma:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: ${threads}
		gpu: false
	}

	output {
		File gemma_relatedness_matrix="relatedness_matrix.gemma.txt.sXX.txt"
	}
}

task run_gemma_lmm {
    input:
        File gemma_out=gemma.imputed.out.txt
        File imputed_prefix=/mnt/data/gemma/chr21-filtered.matched 

        Int? memory = 16
	    Int? disk = 200
        Int? threads = 8

    command {
        gemma \
           -bfile ${imputed_prefix} \
           -k relatedness_matrix.sXX.txt \
           -lmm 4 \
           -o gemma_lmm4
    }

    runtime {
		docker: "quay.io/shukwong/gemma:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: ${threads}
		gpu: false
	}

	output {
		File gemma_association_file="gemma_lmm4.assoc.txt"
        File gemma_log_file="gemma_lmm4.log.txt"
	}
}