require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3 ) {
   stop("Arguments to combine_covars.R: \
     covar_filename pc_filename  covar_file_id_colname" )
}

covar_filename <- args[1]
pc_filename <- args[2] 
covar_file_id_colname <- args[3]

     
## try to read phenotype data
covars_data <- read_delim(covar_filename, col_names =TRUE, delim="\t")
if (length(unique(colnames(covars_data))) != length(colnames(covars_data))) {
   stop(paste("duplicate column names detected in phenotype file \"", covar_filename, "\"", sep=""))
}


## read PC file, we assume this came from plink
pcs <- read_delim(pc_filename, col_names=TRUE, delim="\t") 
pcs <- pcs %>% rename (FID=`#FID`, !!covar_file_id_colname := IID)


covar_with_PCs <- covars_data %>% 
            inner_join(pcs, by=covar_file_id_colname) %>%
            select(FID, everything())


## subset to complete cases, which are all that are used by SAIGE
#covar_with_PCs <- covar_with_PCs[complete.cases(covar_with_PCs),]

#uoutput
write.table(covar_with_PCs, "covar_with_pcs.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
