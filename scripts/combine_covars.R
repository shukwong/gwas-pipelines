library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3 ) {
   stop("Arguments to combine_covars.R: \
     covar_filename pc_filename  covar_file_id_colname" )
}

covar_filename <- args[1]
pc_filename <- args[2] 
covar_file_id_colname <- args[3]

     
## try to read phenotype data
covars_data <- fread(covar_filename, header = TRUE)
if (length(unique(colnames(covars_data))) != length(colnames(covars_data))) {
   stop(paste("duplicate column names detected in phenotype file \"", covar_filename, "\"", sep=""))
}


## read PC file, we assume this came from plink
pcs <- fread(pc_filename, header=TRUE) 
pcs <- pcs %>% rename (FID=`#FID`, !!covar_file_id_colname := IID)


covar_with_PCs <- covars_data %>% 
            inner_join(pcs, by=covar_file_id_colname) 

if("IID" %in% colnames(covar_with_PCs)) {
  covar_with_PCs <- covar_with_PCs %>% select (FID, IID, everything())
} else {
  covar_with_PCs <- covar_with_PCs %>%
                    mutate (IID = FID) %>%
                    select (FID, IID, everything())
}

#get PCs as string 
pcs_as_string <-  glue::glue_collapse(colnames(pcs)[grep("^PC",(colnames(pcs)))], ",")

#output
write.table(covar_with_PCs, "covar_with_pcs.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
writeLines(pcs_as_string, "pcs_as_string.txt")
