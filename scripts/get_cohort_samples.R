#This is to create the sample overlap between samples with available phenotypes, 
#samples with available genotypes, and samples with available imputed genotypes
#may extend to more than 3 input sample files

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4 ) {
  stop("Arguments to get_cohort_sample.R: \
        covariate_tsv_file genotype_samples_to_keep_file imputed_samples_to_keep_file sampleIDcol" )
}

covariate_tsv_file <- args[1]
genotype_samples_to_keep_file <- args[2]
imputed_samples_to_keep_file <- args[3]
sampleID = args[4]


covariates <- fread(covariate_tsv_file)
genotype_samples <- fread(genotype_samples_to_keep_file)
imputed_samples <- fread(imputed_samples_to_keep_file)

covariants_matched <- covariates %>% filter (!!sym(sampleID) %in% genotype_samples$IID)  %>%               
                      filter (!!sym(sampleID) %in% imputed_samples$IID )


plink_subset_samples <- covariants_matched %>% mutate (FID = !!sym(sampleID), IID = !!sym(sampleID)) %>% select (FID,IID)

#bgen_subset_samples <- covariants_matched %>% pull (!!sym(sampleID)) 

write.table(covariants_matched, "covars_subsetted.tsv", quote=F, sep=" ",
            col.names = T, row.names = F)

write.table(plink_subset_samples, "plink_subsetted.samples", quote=F, sep=" ",
            col.names = T, row.names = F)

write.table(plink_subset_samples$IID, "subsetted.samples", quote=F, sep=" ",
            col.names = F, row.names = F)

# write_lines("ID_1 ID_2", "bgen_subsetted.samples")
# write_lines("0 0", "bgen_subsetted.samples", append = TRUE)
# write.table(plink_subset_samples, "bgen_subsetted.samples", quote=F, sep=" ",
#             col.names = F, row.names = F, append = TRUE)


