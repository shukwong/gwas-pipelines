#sampleID column is required in covariate_tsv_file, variable info is defined in variable_info_tsv_file

require(tidyverse)
require(jsonlite)

if (0) {
args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2 ) {
    stop("Arguments to combine_covars.R: \
          create_covar_files_by_set.R covariate_tsv_file variable-info_tsv_file sample_sets_json_file" )
  }

covariate_tsv_file <- args[1]
variable_info_tsv_file <- args[2]
sample_sets_json_file <- args[3]
}

covariate_tsv_file="variables/sample_covars.tsv"
variable_info_tsv_file="variables/variable-info.tsv" 
sample_sets_json_file="variables/samples.json"

covariates <- read_delim(covariate_tsv_file, delim="\t")
variable_info <- read_delim(variable_info_tsv_file, delim="\t")
sample_sets <- fromJSON(sample_sets_json_file)

#process sample sets
for (i in 1:length(sample_sets)) {
  sample_set_name = names(sample_sets)[i]
  sample_set <- sample_sets[[i]]
  covariates_current_set = covariates
  
  #start the log file
  covar_set_log <- paste0(sample_set_name, ".log")
  writeLines(paste0("Full covar file has ", nrow(covariates_current_set), "samples."),      
             covar_set_log)
  
  #filter on other criteria, use a for loop for now
  for (i in 1:nrow(sample_set)) {
    #if (varname=="samples_included") {next}
    
    varname = sample_set$varname[i]
    criteria = sample_set$criteria[i]
    value = sample_set$value[i]
    
    #if samples_included is specified, we would read in the samples file
    #we assume that sample names are unique and hence we are expecting a single column 
    #file without header
    if (varname == "samples") {
      samples_filename <- sample_set %>% 
            filter (varname == "samples_included") %>% 
            pull (criteria)
      samples_to_include <- read_delim(samples_filename, delim=" ", col_names = F) %>%
            pull (X1)
      covariates_current_set <- covariates_current_set %>% filter (sampleID %in% samples_to_include)
      writeLines(paste0("samples from ", samples_filename, " are included. There are ", nrow(covariates_current_set), "samples left. "), covar_set_log)
    } else if (criteria == "include") {
      covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] %in% value),]
      writeLines(paste0("samples with ", varname, "==", value, " are included. There are ", nrow(covariates_current_set), "samples left."), covar_set_log)
      
    } else if (criteria == ">") { # a bit stupid to separate them but not sure how to do this elegantly
      
      covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] > value),]
      writeLines(paste0("samples with ", varname, ">", value, " are included. There are ", nrow(covariates_current_set), "samples left."), covar_set_log)
      
    } else if (criteria == "<") {
      
      covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] < value),]
      writeLines(paste0("samples with ", varname, "<", value, " are included. There are ", nrow(covariates_current_set), "samples left."), covar_set_log)
      
    } else if (criteria == "==") {
      
      covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] == value),]
      writeLines(paste0("samples with ", varname, "==", value, " are included. There are ", nrow(covariates_current_set), "samples left."), covar_set_log)
      
    }
    
    #output the file
    outfile_name <- paste0("covars_", sample_set_name, ".tsv")
    write.table(covariates_current_set, outfile_name, sep="\t", col.names = T,
                row.names = F, sep="\t")
  }
  
      
}
