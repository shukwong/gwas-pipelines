require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4 ) {
   stop("Arguments to combine_covars.R: \
     covar_filename pc_filename  covar_file_id_colname output_filename" )
}

covar_filename <- args[1]
pc_filename <- args[2] 
covar_file_id_colname <- args[3]
output_filename <- args[4]

     
## try to read phenotype data
covars_data <- read_delim(covar_filename, col_names =TRUE, delim="\t")
if (length(unique(colnames(covars_data))) != length(colnames(covars_data))) {
   stop(paste("duplicate column names detected in phenotype file \"", covar_filename, "\"", sep=""))
}


## read PC file, we assume this came from plink
pcs <- read_delim(pc_filename, col_names=TRUE, delim="\t") 
pcs <- pcs %>% select (-`#FID`) %>% rename (!!covar_file_id_colname := IID)


covar_with_PCs <- covars_data %>% 
            inner_join(pcs, by=(!!covar_file_id_colname)) 


## subset to complete cases, which are all that are used by SAIGE
#covar_with_PCs <- covar_with_PCs[complete.cases(covar_with_PCs),]

#uoutput
write.table(covar_with_PCs, output_filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
