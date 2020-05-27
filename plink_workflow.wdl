
workflow run_preprocess {
    
	File genotype_bed
	File genotype_bim
	File genotype_fam
	
    File chain_file

    File imputed_files_dir

	Int? memory = 60
	Int? disk = 500
	Int? threads = 16


    call task plink_pca {
 
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