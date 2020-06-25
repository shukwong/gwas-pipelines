require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6 | length(args) > 7) {
   stop("Arguments to combine_covars.R: \
     covar_filename pc_filename phenotype covar_file_id_colname covariate_list output_filename [sample_filename]" )
}

covar_filename <- args[1]
pc_filename <- args[2]
phenotype_name <- args[3]
covar_file_id_colname <- args[4]
covariate_list <- unlist(strsplit(args[5], ","))
covariate_list <- unique(covariate_list[covariate_list != "NA"])
output_filename <- args[6]

     

## try to read phenotype data
covars_data <- read_delim(covar_filename, col_names =TRUE, delim="\t")
if (length(unique(colnames(covars_data))) != length(colnames(covars_data))) {
   stop(paste("duplicate column names detected in phenotype file \"", covar_filename, "\"", sep=""))
}
## enforce IDs present once
if (length(which(colnames(covars_data) == covar_file_id_colname)) != 1) {
   stop(paste("header of phenotype data does not contain exactly one ", covar_file_id_colname, " column", sep=""))
}
## enforce phenotype present once
if (length(which(colnames(covars_data) == phenotype_name)) != 1) {
   stop(paste("header of phenotype data does not contain requested phenotype \"", phenotype_name, "\" exactly once", sep=""))
}

## read PC file
pcs <- read_delim(pc_filename, col_names=TRUE, delim="\t") 
pcs <- pcs %>% select (-`#FID`)
# figure out how many PCs are there and name them
#num_pcs <- pcs %>% select (-V1, -V2) %>% head() %>% ncol
## we only want the IID and the PC columns
#pcs <- pcs %>% set_names(c("FID", "IID", paste0("PC", 1:num_pcs))) %>%
#      select (IID, starts_with("PC"))


covar_with_PCs <- covars_data %>% 
            select (covar_file_id_colname, phenotype_name, covariate_list) %>% 
            rename (IID = !!covar_file_id_colname) %>%
            inner_join(pcs, by=("IID")) %>%
            mutate (FID = IID) %>%
            select (FID, IID, everything())


## subset to complete cases, which are all that are used by SAIGE
covar_with_PCs <- covar_with_PCs[complete.cases(covar_with_PCs),]

#subset samples if sample file is given
if (length(args)==7) {
  sample_filename <- args[7]
  samples <- read_delim(sample_filename, delim="\t", col_names = F)
  covar_with_PCs <- covar_with_PCs %>% filter (IID %in% samples$X1)
}  


#uoutput
write.table(covar_with_PCs, output_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
