workflow fastgwa_workflow {
    
	File genotype_bed
	File genotype_bim
	File genotype_fam

    File imputed_samples_file
    File pheno_file
    File covar_file
    File bgen_list_file

    String pheno_col

    Int num_parts_make_grm_chunk = 4

    scatter (part_number in range(num_parts_make_grm_chunk)) {
        call gcta_make_grm_chunk {
            input:
                genotype_bed = genotype_bed,
	            genotype_bim = genotype_bim,
	            genotype_fam = genotype_fam,
                num_parts = num_parts_make_grm_chunk,
                part_number = part_number
        }
    }

    call gcta_merge_and_create_sparse_grm {
        input:
	        fastgwa_grm_part_bin_files = gcta_make_grm_chunk.fastgwa_grm_part_bin_file,
	        fastgwa_grm_part_id_files = gcta_make_grm_chunk.fastgwa_grm_part_id_file,
            fastgwa_grm_part_Nbin_files = gcta_make_grm_chunk.fastgwa_grm_part_Nbin_file

    }

    Array[Array[String]] bgen_file_list = read_tsv(bgen_list_file)
    scatter (bgen_file_line in bgen_file_list) {
        call run_fastgwa {
            input:
                 fastgwa_sparse_grm_bin_file = gcta_merge_and_create_sparse_grm.fastgwa_sparse_grm_bin_file,
	             fastgwa_sparse_grm_id_file = gcta_merge_and_create_sparse_grm.fastgwa_sparse_grm_id_file,
                 fastgwa_sparse_grm_Nbin_file = gcta_merge_and_create_sparse_grm.fastgwa_sparse_grm_Nbin_file,
                 bgen_file = bgen_file_line[0],
                 bgen_samples_file = imputed_samples_file,
                 pheno_file = pheno_file,
                 covar_file = covar_file,
                 pheno_col = pheno_col
        }
    }

    call merge_fastgwa_results {
        input:
            fastgwa_result_file = run_fastgwa.fastgwa_result_file,
            pheno_col = pheno_col
    }

    output {
		File merged_fastgwa_result_file = merge_fastgwa_results.merged_fastgwa_result_file
	}

}    


task gcta_make_grm_chunk {

	File genotype_bed
	File genotype_bim
	File genotype_fam
    
    Int num_parts #this is the number of of parts we want to parallelization 
    Int part_number 
    
	Int? memory = 16
	Int? disk = 100
    Int? threads = 4
    Int? preemptible_tries = 3

	command {
        gcta64 \
            --bfile ${sub(genotype_bed,".bed",'')} \
            --make-grm-part ${num_parts} ${part_number} --thread-num ${threads} \
            --out fastgwa_grm_part${part_number}_of_${num_parts}
     
	}

	runtime {
		docker: "quay.io/shukwong/gcta:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: preemptible_tries
	}

	output {
		File fastgwa_grm_part_bin_file = "fastgwa_grm_part${part_number}_of_${num_parts}.bin"
        File fastgwa_grm_part_id_file = "fastgwa_grm_part${part_number}_of_${num_parts}.id"
        File fastgwa_grm_part_Nbin_file = "fastgwa_grm_part${part_number}_of_${num_parts}.N.bin"
	}
}

task gcta_merge_and_create_sparse_grm {

	Array[File] fastgwa_grm_part_bin_files
	Array[File] fastgwa_grm_part_id_files
    Array[File] fastgwa_grm_part_Nbin_files
    
	Int? memory = 8
	Int? disk = 100
    Int? threads = 4
    Int? preemptible_tries = 3

	command {

        cat ${sep=' ' fastgwa_grm_part_bin_files} > fastgwa_grm.bin
        cat ${sep=' ' fastgwa_grm_part_id_files} > fastgwa_grm.id
        cat ${sep=' ' fastgwa_grm_part_Nbin_files} > fastgwa_grm.N.bin

        ## turn the results into a sparse GRM
        gcta64 --grm fastgwa_grm \
               --make-bK-sparse 0.05 
               --out fastgwa_sparse_grm
     
	}

	runtime {
		docker: "quay.io/shukwong/gcta:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: preemptible_tries
	}

	output {
		File fastgwa_sparse_grm_bin_file = "fastgwa_sparse_grm.bin"
        File fastgwa_sparse_grm_id_file = "fastgwa_sparse_grm.id"
        File fastgwa_sparse_grm_Nbin_file = "fastgwa_sparse_grm.N.bin"
	}
}


task run_fastgwa {

	File fastgwa_sparse_grm_bin_file
	File fastgwa_sparse_grm_id_file
    File fastgwa_sparse_grm_Nbin_file
    File bgen_file
    File bgen_samples_file
    File pheno_file
    File covar_file

    String pheno_col
    
	Int? memory = 32
	Int? disk = 100
    Int? threads = 8
    Int? preemptible_tries = 3

	command {

        gcta64 \
             --bgen ${bgen_file} \
             --sample ${bgen_samples_file} \
             --fastGWA-mlm \
             --grm-sparse ${sub(fastgwa_sparse_grm_bin_file,".bin",'')} \
             --pheno ${pheno_file} \
             --qcovar ${covar_file} \
             --threads ${threads} \
             --out fastgwa_${pheno_col}
     
	}

	runtime {
		docker: "quay.io/shukwong/gcta:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: preemptible_tries
	}

	output {
		File fastgwa_result_file = "fastgwa_${pheno_col}.fastGWA"
	}
}


task merge_fastgwa_results {

	Array[File] fastgwa_result_file

    String pheno_col
    
	Int? memory = 8
	Int? disk = 100
    Int? threads = 1
    Int? preemptible_tries = 3

	command {

        echo -e "CHR\tSNP\tPOS\tA1\tA2\tN\tAF1\tBETA\tSE\tP\tINFO" > fastgwa_${pheno_col}.fastGWA.tsv
        
        cat ${sep=' ' fastgwa_result_file} | gzip -d | grep -v ^CHR >> fastgwa_${pheno_col}.fastGWA.tsv
        
        gzip fastgwa_${pheno_col}.fastGWA.tsv
     
	}

	runtime {
		docker: "quay.io/shukwong/gcta:latest"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
        cpu: "${threads}"
		preemptible: preemptible_tries
	}

	output {
		File merged_fastgwa_result_file = "fastgwa_${pheno_col}.fastGWA.tsv.gz"
	}
}