# This is a GWAS pipeline that is intended to be run with biobank size data

The workflow is designed to be modular with a main workflow file that calls the preprocessing subworkflow and (optionally) the GWAS tools.
Currently it supports Bolt-LMM, SAIGE and fastGWA from GCTA. 

## Inputs for the workflow

### genotype and sample files (required) 
* *genotype_bed*, *genotype_bim* and *genotype_fam* are directly genotype files in Plink format

* genotype_samples_to_keep_file: subset of samples to be used  from the genotype file. Need to have FID and IID columns.
* imputed_samples_to_keep_file: subset of samples to be used from the imputed genotype file. Need to have FID and IID columns. 

* covariate_tsv_file: File that contains the sample ID, phenotype, and the covariates used in the association study.
* variable_info_tsv_file: File that defines the columns in the covariate_tsv_file.
* sample_sets_json_file: File that defines the subsets used for the association study. 

### Which GWAS tool to use 
* useBOLT: true or false
* useSAIGE: true or false
* if bolt is used, it requires the genetic_map_file and ld_scores_file

### Imputed files    
The workflow can handle either VCF or bgen files as inputs for the imputed genotype data. The list of the file and their location can be specified as *imputed_list_of_vcf_file* or *imputed_list_of_bgen_file*.

In case that the VCF file is used, the user can provide the *id_delim* string if FID and IID are concat into one ID, if *id_delim* is not provided, double ID is assumed - both FID and IID are the same as sample ID in VCF.    

Finally, if the imputed data and the genotype data are of different genome build, the workflow can liftover the genotype files to match with the imputed data if the user provides the *chain_file*.

## workflow steps

* Determine the sample subsets and create covar files for each subset
* The tasks are scattered for each sample subset:
  * The preprocessing steps are summarized in the DAG. 
![DAG for preprocessing pipeline](workflow_diagrams/preprocessing-pipeline.png)
  * Run on each of the association tool specified
  * Gather the results
  * Make Manhattan and QQ plots 

## Limitations

* Right now this pipeline only works on the 22 human autosomes
* Variable transformation is not yet handled
* The workflow currently points to the latest version of subworkflow on GitHub, this needs to be changed to a particular release instead.
