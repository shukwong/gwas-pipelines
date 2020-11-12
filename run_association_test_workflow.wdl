version development

 # import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/tasks/preprocess_workflow.wdl" as preprocess
 # import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/tasks/saige_workflow.wdl" as saige
 # import "https://raw.githubusercontent.com/shukwong/gwas-pipelines/master/tasks/bolt_workflow.wdl" as bolt

import "./tasks/preprocess_workflow.wdl" as gwas_tasks
import "./tasks/saige_workflow.wdl" as saige
import "./tasks/bolt_workflow.wdl" as bolt

workflow run_association_test {
    input {
    File genotype_bed
    File genotype_bim
    File genotype_fam

    File genotype_samples_to_keep_file
    File imputed_samples_to_keep_file

    File covariate_tsv_file
    
    String binary_covar_list
    String continuous_covar_list 
    String phenoCol 
    String covar_sampleID_colname
    String phenotype_type
    String setname 
    String batch_name

    Boolean? useBOLT
    Boolean? useSAIGE

    #for bolt
    File? genetic_map_file
    File? ld_scores_file
    
    Float? minMAF=0.001
    Float? minMAC=1

    File imputed_list_of_files
    String impute_file_format
    File? chain_file

    String? id_delim #delim character for vcf file, if not defined, double ID is assumed
    String? dosageField
    }

    # Array[String] binary_covar_list_lines  = read_lines(get_covar_subsets.binary_covar_list_file)
    # String binary_covar_list = binary_covar_list_lines[0]
    # Array[String] continuous_covar_list_lines = read_lines(get_covar_subsets.continuous_covar_list_file)
    # String continuous_covar_list = continuous_covar_list_lines[0]
    # Array[String] phenoCol_lines = read_lines(get_covar_subsets.phenotype_line_file)
    # String phenoCol = phenoCol_lines[0]
    # Array[String] covar_sampleID_colname_lines = read_lines(get_covar_subsets.sampleid_line_file)
    # String covar_sampleID_colname = covar_sampleID_colname_lines[0]
    # Array[String] phenotype_type_lines = read_lines(get_covar_subsets.phenotype_type_file)
    # String phenotype_type = phenotype_type_lines[0]


    call gwas_tasks.run_preprocess {
        input:
            genotype_bed = genotype_bed,
            genotype_bim = genotype_bim,
            genotype_fam = genotype_fam,

            genotype_samples_to_keep_file = genotype_samples_to_keep_file,
            imputed_samples_to_keep_file = imputed_samples_to_keep_file,
            covariate_tsv_file = covariate_tsv_file,
            covar_sampleID_colname = covar_sampleID_colname,
            imputed_list_of_files = imputed_list_of_files, 
            impute_file_format = impute_file_format,
            chain_file = chain_file,
            id_delim = id_delim,
            dosageField = dosageField
    }

    Array[String] pcs_as_string_lines = read_lines(run_preprocess.pcs_as_string_file)
    String pcs_as_string = pcs_as_string_lines[0]

    #String setname = basename(covar_subset_file, "_covars.tsv")

    if (defined(useSAIGE) && useSAIGE) {
        call saige.run_saige {
            input:
                genotype_bed = run_preprocess.genotype_ready_bed,
                genotype_bim = run_preprocess.genotype_ready_bim,
                genotype_fam = run_preprocess.genotype_ready_fam,
                bgen_files_and_indices = run_preprocess.bgen_files_and_indices,
                imputed_samples_file = run_preprocess.bgen_samples,
                phenoCol = phenoCol,
                phenotype_type = phenotype_type,
                covar_file = run_preprocess.covar_file,
                covarColList = binary_covar_list + "," + continuous_covar_list + "," + pcs_as_string,
                setname = setname,
                batch_name = batch_name,
                minMAF=minMAF,
                minMAC=minMAC
        }

        call gwas_tasks.make_summary_plots as make_saige_plots {
            input:
                association_summary_file = run_saige.merged_saige_file
            }
        }

        if (defined(useBOLT) && useBOLT ) {
            call bolt.bolt_workflow {
                input: 
                    genotype_bed = run_preprocess.genotype_ready_bed,
	                genotype_bim = run_preprocess.genotype_ready_bim,
	                genotype_fam = run_preprocess.genotype_ready_fam,
                    pheno_file = run_preprocess.covar_file,
                    ld_scores_file = select_first([ld_scores_file]),
                    genetic_map_file = select_first([genetic_map_file]),
                    #imputed_samples_file = run_preprocess.bgen_samples,
                    covar_file = run_preprocess.covar_file,
                    #bgen_list_file = run_preprocess.bgen_paths_file,
                    bgen_files_and_indices = run_preprocess.bgen_files_and_indices,
                    pheno_col = phenoCol,
                    qCovarCol = continuous_covar_list + "," + pcs_as_string,
                    setname = setname,
                    batch_name = batch_name,
                    minMAF = minMAF
            }

        call gwas_tasks.make_summary_plots as make_bolt_plots {
            input: 
                association_summary_file = bolt_workflow.imputed_stats_file
        }
                
    }

 
	output {
       File merged_saige_file = select_first([run_saige.merged_saige_file, "null"])
       File merged_bolt_file = select_first([bolt_workflow.imputed_stats_file, "null"])
    
       File? saige_manhattan_plot =  make_saige_plots.manhattan_file
       File? saige_manhattan_loglog_plot = make_saige_plots.manhattan_loglog_file
       File? saige_qqplot = make_saige_plots.qqplot_file


       File? bolt_manhattan_plots =  make_bolt_plots.manhattan_file
       File? bolt_manhattan_loglog_plots = make_bolt_plots.manhattan_loglog_file
       File? bolt_qqplots = make_bolt_plots.qqplot_file
 	}
    

    meta {
		author : "Wendy Wong"
		email : "wendy.wong@gmail.com"
		description : "Biobank scale association study with mutiple subsets (sample target variable)."
	}

}

