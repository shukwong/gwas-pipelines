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

task lift_over_plink {
    

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