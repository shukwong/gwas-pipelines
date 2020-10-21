version development

workflow run_preprocess {
    input {
        File genotype_bed
        File genotype_bim
        File genotype_fam

        File genotype_samples_to_keep_file
        File imputed_samples_to_keep_file
        File covariate_tsv_file
        File imputed_list_of_files  
        String impute_file_format 

        String covar_sampleID_colname
    
        File? chain_file
        #File? imputed_list_of_vcf_file
        #File? imputed_list_of_bgen_file
        #File? imputed_list_of_bgen_index_file

        String? dosageField
        String? id_delim
    }    

    call get_cohort_samples {
        input: 
            covariate_tsv_file = covariate_tsv_file,
            genotype_samples_to_keep_file = genotype_samples_to_keep_file,
            imputed_samples_to_keep_file = imputed_samples_to_keep_file,
            covar_sampleID_colname = covar_sampleID_colname
    }

    call plink_subset_sample {
        input:
             genotype_bed = genotype_bed, 
             genotype_bim = genotype_bim,
             genotype_fam = genotype_fam,
             samples_to_keep_file = get_cohort_samples.plink_subset_samples
    }

    call run_genotype_qc_filter {
        input:
            genotype_bed = plink_subset_sample.genotype_subsetSample_bed,
            genotype_bim = plink_subset_sample.genotype_subsetSample_bim,
            genotype_fam = plink_subset_sample.genotype_subsetSample_fam

    }

    call run_ld_prune {
        input:
            genotype_bed = run_genotype_qc_filter.gentoype_qc_filtered_bed,
            genotype_bim = run_genotype_qc_filter.gentoype_qc_filtered_bim,
            genotype_fam = run_genotype_qc_filter.gentoype_qc_filtered_fam
    }

    call plink_pca {
        input:
    	    genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam
            
    }

    call addPCs_to_covar_matrix {
         input:
             covar_file = covariate_tsv_file, 
             plink_pca_eigenvec_file = plink_pca.genotype_pruned_pca_eigenvec,
             covar_sampleID_colname = covar_sampleID_colname
    }

    if (defined(chain_file)) {

        call liftover_plink {
            input:
    	    genotype_bed = run_ld_prune.genotype_pruned_bed,
    	    genotype_bim = run_ld_prune.genotype_pruned_bim,
    	    genotype_fam = run_ld_prune.genotype_pruned_fam,
            chain_file = select_first([chain_file])
        }
    }

    
    
    if (impute_file_format == "vcf") {

        Array[Array[File]] imputed_files = read_tsv(imputed_list_of_files)

        scatter (imputed_file in imputed_files) {
	        call vcf_to_bgen {
	            input:
                    vcf_file = imputed_file[0],
                    samples_to_keep_file = get_cohort_samples.plink_subset_samples,
                    dosageField = dosageField,
                    id_delim = id_delim
		        }
	    }

        Array[Array[File]] converted_bgen_file_list = transpose([vcf_to_bgen.bgen_file, vcf_to_bgen.bgen_file_index]) 

        Array[File] bgen_file_saige_samples = vcf_to_bgen.bgen_file_saige_sample

        File bgen_file_saige_sample = bgen_file_saige_samples[0]

    }  
    
    if (impute_file_format == "bgen") {
        Array[Array[File]]  bgen_file_list = read_tsv(imputed_list_of_files)
    }


    # Map[Array[String], Array[String]] bgen_files_and_indices = as_map(zip(converted_bgen_file_list,bgen_file_list))
    

	output {
          File genotype_ready_bed = select_first([liftover_plink.output_bed, run_ld_prune.genotype_pruned_bed])
          File genotype_ready_bim = select_first([liftover_plink.output_bim, run_ld_prune.genotype_pruned_bim])
          File genotype_ready_fam = select_first([liftover_plink.output_fam, run_ld_prune.genotype_pruned_fam])
          #File genotype_pruned_pca_eigenvec = plink_pca.genotype_pruned_pca_eigenvec
          #Array[File] bgen_files = select_first([vcf_to_bgen.out_bgen, imputed_bgen_files])
          #Array[File] bgen_file_indices = select_first([index_bgen_file.bgen_file_index, imputed_bgen_index_files])
          Array[Array[File]] bgen_files_and_indices = select_first([converted_bgen_file_list,bgen_file_list])          
          #File bgen_paths_file = select_first([imputed_list_of_bgen_file,get_bgen_file_paths.merged_file])
          File bgen_samples = select_first([bgen_file_saige_sample, imputed_samples_to_keep_file])
          File plink_subset_samples = get_cohort_samples.plink_subset_samples
          File covar_file = addPCs_to_covar_matrix.covar_file_with_pcs
          File pcs_as_string_file = addPCs_to_covar_matrix.pcs_as_string_file
 	}
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Preprocess Genotype files for biobank scale GWAS study."
	}

    

}

task plink_pca {
    input {
        File genotype_bed
        File genotype_bim
        File genotype_fam

        String? approx = "approx"

        Int? memory = 60
        Int? disk = 200
    }

    command <<<
        set -euo pipefail
        plink2 --bed ~{genotype_bed} --bim ~{genotype_bim} --fam ~{genotype_fam} --pca ~{approx} --out genotype_pruned_pca
	>>>

	runtime {
		docker: "quay.io/shukwong/plink_crossmap_bgen:8984373caf8b"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
		gpu: false
	}

    output {
	    File genotype_pruned_pca_eigenvec = "genotype_pruned_pca.eigenvec"
	    File genotype_pruned_pca_eigenval = "genotype_pruned_pca.eigenval"
        File genotype_pruned_pca_log = "genotype_pruned_pca.log"
    }
}

task run_genotype_qc_filter {
	input {
        File genotype_bed
        File genotype_bim
        File genotype_fam

        Int? memory = 60
        Int? disk = 200
    }

    command <<<
        set -euo pipefail
        plink2 --bed ~{genotype_bed} --bim ~{genotype_bim} --fam ~{genotype_fam} \
              --maf 0.01 --geno 0.02 --hwe 0.001 --make-bed --out gentoype_qc_filtered

    >>>

    runtime {
	    docker: "quay.io/shukwong/plink_crossmap_bgen:8984373caf8b"
        memory: memory + " GiB"
        disks: "local-disk " + disk + " HDD"
        gpu: false
    }
	
    output {
	    File gentoype_qc_filtered_bed = "gentoype_qc_filtered.bed"	
        File gentoype_qc_filtered_bim = "gentoype_qc_filtered.bim"
        File gentoype_qc_filtered_fam = "gentoype_qc_filtered.fam"
    } 	
}


task run_ld_prune {
    input {
        File genotype_bed
        File genotype_bim
        File genotype_fam

        Int? memory = 32
        Int? disk = 200
    }

    command <<<
	
        plink --bed ~{genotype_bed} --bim ~{genotype_bim} --fam ~{genotype_fam}  --indep 50 5 2 --out ld_indep_check

        plink --keep-allele-order --bed ~{genotype_bed} --bim ~{genotype_bim} \
              --fam ~{genotype_fam} --extract ld_indep_check.prune.in \
              --make-bed --out ld_indep_check.prune

        plink --bfile ld_indep_check.prune  --indep-pairwise 50 5 0.2 \
              --out ld_indep_pairwise_check

        plink --keep-allele-order --bfile ld_indep_check.prune \
              --extract ld_indep_pairwise_check.prune.in \
              --make-bed --out genotype_pruned_plink

        plink --bfile genotype_pruned_plink --het --out genotype_pruned_plink_het
        awk 'NR > 1 && sqrt($6^2) > sqrt(0.2^2) {print $1"\t"$1}' genotype_pruned_plink_het.het > genotype_pruned_plink_het.het.remove      
        ## make a version of the data with QCed autosomal data for later
        plink --bed ~{genotype_bed} --bim ~{genotype_bim} --fam ~{genotype_fam} \
              --remove genotype_pruned_plink_het.het.remove --autosome \
              --make-bed --out genotype.nohet.autosomes
    >>>

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
		gpu: false
	}

    output {
	    File genotype_pruned_bed = "genotype_pruned_plink.bed"
        File genotype_pruned_bim = "genotype_pruned_plink.bim"
        File genotype_pruned_fam = "genotype_pruned_plink.fam"
        File genotype_nohet_autosomes_bed = "genotype.nohet.autosomes.bed"
        File genotype_nohet_autosomes_bim = "genotype.nohet.autosomes.bim"
        File genotype_nohet_autosomes_fam = "genotype.nohet.autosomes.fam"
    }
}

task plink_to_vcf {
    input {
	    File genotype_bed
	    File genotype_bim
	    File genotype_fam

        String prefix = basename (genotype_bed, ".bed")
	
	    Int? memory = 32
	    Int? disk = 200
    }

	command <<<
        plink --bed ~{genotype_bed} --bim ~{genotype_bim} --fam ~{genotype_fam} --recode vcf --out ~{prefix}.vcf   
	>>>

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
		gpu: false
	}

	output {
		File out_vcf = prefix + ".vcf"
	}
}

task vcf_to_plink_bed {
    input {
	    File vcf_file
        String prefix = basename(vcf_file, ".vcf")
	    Int? memory = 32
	    Int? disk = 200
    }

	command <<<
		plink --vcf ~{vcf_file}  --make-bed --out ~{prefix}
	>>>

	runtime {
		docker: "quay.io/h3abionet_org/py3plink"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
		gpu: false
	}

	output {
		File out_bed = prefix + ".bed"
		File out_bim = prefix + ".bim"
		File out_fam = prefix + ".fam"
	}
}

task vcf_to_bgen {
    input {
        File vcf_file
        File samples_to_keep_file

        String? dosageField = "HDS"
        String? id_delim

        Int? bits=8
        Int? memory = 32
        Int? disk = 500

        Int? threads = 4

        String prefix = basename(vcf_file, ".vcf.gz")

        String plink_id_delim_option = if defined (id_delim) then "--id-delim " + id_delim else "--double-id"
    }

    #to figure out the disk size based on the vcf file size
    Float vcf_file_size = size(vcf_file, "GiB")
    Int disk_size = select_first([disk, ceil(vcf_file_size * 4)])

	command <<<
        set -euo pipefail

        plink2 --threads ~{threads} --vcf ~{vcf_file} dosage=~{dosageField} ~{plink_id_delim_option} \
            --make-pgen erase-phase --out plink_out

        plink2 --threads ~{threads} --pfile plink_out --keep ~{samples_to_keep_file} \
            --export bgen-1.2 bits=~{bits} id-paste=iid --out ~{prefix}   

        bgenix -g ~{prefix}.bgen -index -clobber

        ## mysterious file format discrepancy causes failures with FASTGWA; change NA values to 0 in sample files
        sed 's/NA$/0/' ~{prefix}.sample > ~{prefix}-noNAs.sample    

        awk 'NR > 2 {print $1}' ~{prefix}-noNAs.sample  > ~{prefix}_bgen.sample #samples file for saige

	>>>

	runtime {
		docker: "quay.io/shukwong/plink_crossmap_bgen:8984373caf8b"
		memory: memory + " GB"
		disks: "local-disk " + disk_size + " HDD"
        cpu: threads
		gpu: false
	}

	output {
		File bgen_file = prefix + ".bgen"
        File bgen_file_index = prefix + ".bgen.bgi"
		File bgen_file_sample = prefix + "-noNAs.sample"
		File bgen_file_log = prefix + ".log"
        File bgen_file_saige_sample = prefix + "_bgen.sample"
        # File bgen_file_paths = "bgen_file_paths.txt"
	}
}


task liftover_plink {
    input {
   	    File genotype_bed
	    File genotype_bim
	    File genotype_fam

        File chain_file

        Int? memory = 32
        Int? disk = 200
    }

    command <<<
		
        awk 'start=$4-1 {print "chr"$1"\t"start"\t"$4"\t"$2"\t"$4"\t"$5"\t"$6}' ~{genotype_bim} >bim_as_bed.bed

        CrossMap.py bed ~{chain_file} bim_as_bed.bed bim_as_bed.crossmap.bed 

        grep -v ^chrUn bim_as_bed.crossmap.bed | grep -v random \
            | grep -v alt | grep -v ^chrX | grep -v ^chrY | sort -k1,1 -k2,2g \
            | cut -f4 >bim_as_bed.mapped.ids

        plink2 \
            --bed ~{genotype_bed} --bim ~{genotype_bim} \
            --fam ~{genotype_fam} --extract bim_as_bed.mapped.ids \
            --make-bed --out genotypes_updated    
    >>>        

	runtime {
		docker: "quay.io/shukwong/plink_crossmap_bgen:8984373caf8b"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
		gpu: false
	}

    output {
	    File output_bed = "genotypes_updated.bed"
        File output_bim = "genotypes_updated.bim"
        File output_fam = "genotypes_updated.fam"
    }
}



task plink_subset_sample {
    input {
        File genotype_bed
        File genotype_bim
        File genotype_fam
        File samples_to_keep_file

        Int? memory = 32
        Int? disk = 500
    }
	
    command <<<
		plink2 --bed ~{genotype_bed} --bim ~{genotype_bim} --fam ~{genotype_fam} --keep ~{samples_to_keep_file} --make-bed --out genotype_subsetSample
    >>>

	runtime {
		docker: "quay.io/shukwong/plink_crossmap_bgen:8984373caf8b"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
		gpu: false
	}

    output {
	    File genotype_subsetSample_bed = "genotype_subsetSample.bed"
        File genotype_subsetSample_bim = "genotype_subsetSample.bim"
        File genotype_subsetSample_fam = "genotype_subsetSample.fam"
    }
}

task match_genotype_and_imputed_samples {
    input {
        File genotype_bed
        File genotype_bim
        File genotype_fam
        File imputed_samples_file

        Int? memory = 32
        Int? disk = 200
        Int? threads = 32
    }

    command <<<
        awk '{print $1"\t"$1}' ~{imputed_samples_file}  > samples_plink_format.txt

        /plink2 --bed ~{genotype_bed} --bim ~{genotype_bim} --fam ~{genotype_fam} --keep samples_plink_format.txt \
            --make-bed --out matched_genotype
    >>>

    runtime {
		docker: "quay.io/large-scale-gxe-methods/genotype-conversion"
		memory: memory + " GB"
		disks: "local-disk " + disk + " HDD"
        cpu: threads
		gpu: false
	}

    output {
        File matched_genotype_bed = "matched_genotype.bed"
        File matched_genotype_bim = "matched_genotype.bim"
        File matched_genotype_fam = "matched_genotype.fam"
    }

}

task convert_gen_to_bgen {
    input {
        Array[File] gen_files
        String chrom
        Float threshold
        Boolean? outputBgenOnePointTwo = true

        Int? memory = 32
        Int? disk = 200
        Int? threads = 32
        Int? preemptible_tries = 3
    }

    command <<<
        cat ~{sep=' ' gen_files} > ~{chrom}.gen

        qctool -g ~{chrom}.gen  -threshold ~{threshold}  -filetype gen -ofiletype \
             ~{true='bgen_v1.2 -bgen-bits 8' false='bgen' outputBgenOnePointTwo}  -og ~{chrom}.bgen 
    >>>

    runtime {
		docker: "quay.io/shukwong/qctool:v2.0.8"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: threads
		preemptible: preemptible_tries
	}

    output {
        File converted_bgen_file = chrom + ".bgen"
    }

}



task get_cohort_samples {
    input {
        File covariate_tsv_file
        File genotype_samples_to_keep_file
        File imputed_samples_to_keep_file  
        String covar_sampleID_colname

        Int? memory = 16
        Int? disk = 200
        Int? threads = 1
        Int? preemptible_tries = 3
    }

    command <<<

        wget https://github.com/shukwong/gwas-pipelines/raw/master/scripts/get_cohort_samples.R

        Rscript get_cohort_samples.R ~{covariate_tsv_file} ~{genotype_samples_to_keep_file} ~{imputed_samples_to_keep_file} ~{covar_sampleID_colname}
    >>>

    runtime {
		docker: "rocker/tidyverse:4.0.0"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: threads
		preemptible: preemptible_tries
	}

    output {
        File covar_subset_file = "covars_subsetted.tsv"
        File plink_subset_samples = "plink_subsetted.samples"
        File subset_samples = "subsetted.samples"
    }
}

task get_bgen_file_paths {
    input {
        Array [File] bgen_files
        Array [File] bgen_file_indices
        Array[Array[File]] bgen_zipped = transpose([bgen_files, bgen_file_indices]) 

        Int? memory = 2
        Int? disk = 10
        Int? threads = 1
        Int? preemptible_tries = 3

        #Array[Pair[File, File]] bgen_zipped = zip(bgen_files, bgen_file_indices) 
    
    }

    command <<<

        cat ~{write_tsv(bgen_zipped)} > bgen_file.txt

    >>>

    runtime {
		docker: "ubuntu:18.04"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: threads
		preemptible: preemptible_tries
	}

    output {
        File merged_file = "bgen_file.txt"
    }
}

task addPCs_to_covar_matrix {
    input {
        File covar_file 
        File plink_pca_eigenvec_file
        String covar_sampleID_colname
    
        Int? preemptible_tries = 3
        Int? memory = 8
        Int? disk = 20
        Int? threads = 1
    }

    command <<<
        wget https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/scripts/combine_covars.R

        Rscript combine_covars.R ~{covar_file} ~{plink_pca_eigenvec_file} \
            ~{covar_sampleID_colname}
    >>>

    runtime {
		docker: "rocker/tidyverse:4.0.0"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: threads
		preemptible: preemptible_tries
	}

    output {
        File covar_file_with_pcs =  "covar_with_pcs.tsv"
        File pcs_as_string_file = "pcs_as_string.txt"
    }
}

task make_summary_plots {

    input {
    File association_summary_file

    String? BP_column = "POS"
    String? CHR_column = "CHR"
    String? pval_col = "P"
    String? minrep_col = "SNP"
    Int? loglog_pval=10

    String prefix = basename(association_summary_file, ".tsv.gz")

    Int? memory = 32
    Int? disk = 20
    Int? threads = 2
    Int? preemptible_tries = 3
    }

    command <<<
        wget https://raw.githubusercontent.com/FINNGEN/saige-pipelines/master/scripts/qqplot.R

        R --vanilla -e 'install.packages("qqman",repos = "https://cloud.r-project.org/")'
        R --vanilla -e 'install.packages("optparse",repos = "https://cloud.r-project.org/")'
        R --vanilla -e 'install.packages("R.utils",repos = "https://cloud.r-project.org/")'

        Rscript qqplot.R -f ~{association_summary_file} -o ~{prefix} --chrcol ~{CHR_column} -b POS 
    >>>

    runtime {
		docker: "rocker/tidyverse:3.6.3"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: threads
		preemptible: preemptible_tries 
	}

    output {
        File manhattan_file =  prefix + "_P_manhattan.png"
        File manhattan_loglog_file = prefix + "_P_manhattan_loglog.png"
        File qqplot_file = prefix + "_P_qqplot.png"
    }

}