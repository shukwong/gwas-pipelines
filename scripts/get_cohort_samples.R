#This is to create the sample overlap between samples with available phenotypes, 
#samples with available genotypes, and samples with available imputed genotypes
#may extend to more than 3 input sample files

require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 | length(args) >4 ) {
  stop("Arguments to create_covar_files_by_set.R: \
          get_cohort_samples.R covariate_tsv_file genotype_samples_to_keep_file imputed_samples_to_keep_file <sampleIDcol (if it's not sampleID)>" )
}

covariate_tsv_file <- args[1]
genotype_samples_to_keep_file <- args[2]
imputed_samples_to_keep_file <- args[3]

sampleID = 'sampleID'
if (length(args)==4) {
  sampleID = args[4]
}

covariates <- read_delim(covariate_tsv_file, delim=" ")
genotype_samples <- read_delim(genotype_samples_to_keep_file, delim=" ", col_names = F)
imputed_samples <- read_delim(imputed_samples_to_keep_file, delim=" ", col_names = F)

covariants_matched <- covariates %>% filter (!!sym(sampleID) %in% genotype_samples$X2)  %>%               
                      filter (!!sym(sampleID) %in% imputed_samples$X1 )

write.table(covariants_matched, "covariants_matched.tsv", quote=F, sep=" ",
            col.names = T, row.names = F)