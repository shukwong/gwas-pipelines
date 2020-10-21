version development

 # import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/tasks/preprocess_workflow.wdl" as preprocess
 # import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/tasks/saige_workflow.wdl" as saige
 # import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/tasks/bolt_workflow.wdl" as bolt

import "run_meta_analysis_workflow.wdl" as meta_analysis

workflow gwas_subsets {
    input {
    #File genotype_bed
    #File genotype_bim
    #File genotype_fam

    #File genotype_samples_to_keep_file
    #File imputed_samples_to_keep_file

    File batch_tsv_file 
    File covariate_tsv_file
    File variable_info_tsv_file
    File sample_sets_json_file

   # String phenoCol
   # String covar_sampleID_colname

    Boolean? useBOLT
    Boolean? useSAIGE

    #for bolt
    File? genetic_map_file
    File? ld_scores_file
    
    Float? minMAF=0.001
    Float? minMAC=1

    File? chain_file

    String? id_delim #delim character for vcf file, if not defined, double ID is assumed
    String? dosageField
    }
    #read subset information
    call get_covar_subsets {
        input:
            covariate_tsv_file = covariate_tsv_file, 
            variable_info_tsv_file = variable_info_tsv_file, 
            sample_sets_json_file = sample_sets_json_file
    }
    Array[String] binary_covar_list_lines  = read_lines(get_covar_subsets.binary_covar_list_file)
    String binary_covar_list = binary_covar_list_lines[0]
    Array[String] continuous_covar_list_lines = read_lines(get_covar_subsets.continuous_covar_list_file)
    String continuous_covar_list = continuous_covar_list_lines[0]

    Array[String] phenoCol_lines = read_lines(get_covar_subsets.phenotype_line_file)
    String phenoCol = phenoCol_lines[0]

    Array[String] covar_sampleID_colname_lines = read_lines(get_covar_subsets.sampleid_line_file)
    String covar_sampleID_colname = covar_sampleID_colname_lines[0]

    Array[String] phenotype_type_lines = read_lines(get_covar_subsets.phenotype_type_file)
    String phenotype_type = phenotype_type_lines[0]

    scatter (covar_subset_file in get_covar_subsets.covar_subsets_files) {
        String setname = basename(covar_subset_file, "_covars.tsv")
        call meta_analysis.run_meta_analysis {
            input:
                batch_tsv_file = batch_tsv_file,
                covariate_tsv_file = covar_subset_file,
                variable_info_tsv_file = variable_info_tsv_file,
                binary_covar_list = binary_covar_list,
                continuous_covar_list = continuous_covar_list,
                phenoCol = phenoCol,
                covar_sampleID_colname = covar_sampleID_colname,
                phenotype_type = phenotype_type,
                setname  = setname,

                useBOLT = useBOLT,
                useSAIGE = useSAIGE,

                genetic_map_file = genetic_map_file,
                ld_scores_file = ld_scores_file,
    
                minMAF = minMAF,
                minMAC = minMAC,

                chain_file = chain_file,

                id_delim = id_delim,
                dosageField = dosageField
        }

    }

 
	output {
       Array[File?] merged_saige_file_list = run_meta_analysis.bolt_metal_output_file 
       Array[File?] merged_bolt_file_list = run_meta_analysis.saige_metal_output_file
 	}
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Biobank scale association study with mutiple subsets (sample target variable)."
	}

}



task get_covar_subsets {
    input {
    File covariate_tsv_file 
    File variable_info_tsv_file 
    File sample_sets_json_file

    Int? memory = 4
    Int? disk = 200
    Int? threads = 1
    Int? preemptible_tries = 3
    }

#TODO, change this to git clone a release version when the pipeline is finalized
    command <<< 
       
        wget https://github.com/shukwong/gwas-pipelines/raw/master/scripts/create_covar_files_by_set.R

        R --vanilla -e 'install.packages("fastDummies",repos = "https://cloud.r-project.org/")'

        Rscript create_covar_files_by_set.R ~{covariate_tsv_file} ~{variable_info_tsv_file} ~{sample_sets_json_file} 

    >>>

    runtime {
		docker: "rocker/tidyverse:4.0.0"
		memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: threads
		preemptible: preemptible_tries
	}

    output {
        Array[File] covar_subsets_files = glob("*_covars.tsv")
        Array[File] covar_subsets_log_files = glob("*.log")
        File binary_covar_list_file = "binary_covar_list.txt"
        File continuous_covar_list_file = "continuous_covar_list.txt"
        File phenotype_line_file = "phenotype_line.txt"
        File phenotype_type_file = "phenotype_type.txt"
        File sampleid_line_file = "sampleid_line.txt"
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