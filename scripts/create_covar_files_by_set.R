#This script takes in the files specified by the user and creates different sets of covar files for different runs

#sampleID column is required in covariate_tsv_file, variable info is defined in variable_info_tsv_file

#TODO: variables/variable-info.tsv to transform the variables, currently phenoCol is not defined here but in the main json file
#TODO: add transformation

require(tidyverse)
require(jsonlite)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4 ) {
    stop("Arguments to create_covar_files_by_set.R: \
          create_covar_files_by_set.R covariate_tsv_file variable-info_tsv_file sample_sets_json_file phenoCol" )
}

  
covariate_tsv_file <- args[1]
variable_info_tsv_file <- args[2]
sample_sets_json_file <- args[3]
phenoCol <- args[4]

covariates <- read_table2(covariate_tsv_file)
variable_info <- read_table2(variable_info_tsv_file)
sample_sets <- fromJSON(sample_sets_json_file)

if (ncol(variable_info)!=5) {
  stop ("Error: expecting column names variableName,excluded,variableType,na_values,transformation for the variable-info file\n" )
}

#phenoCol = sample_sets$phenoCol

#process sample sets
for (i in 1:length(sample_sets)) {
#  if (names(sample_sets[1]) == "phenoCol") {next}
  
  sample_set_name = names(sample_sets)[i]
  sample_set <- sample_sets[[i]]
  covariates_current_set = covariates
  
  #start the log file
  covar_set_log <- paste0(sample_set_name, ".log")
  write_lines(paste0("Full covar file has ", nrow(covariates_current_set), " samples."),      
             covar_set_log)
  
  #filter on other criteria, use a for loop for now
  for (j in 1:nrow(sample_set)) {
    #if (varname=="samples_included") {next}
    
    varname = sample_set$varname[j]
    criteria = sample_set$criteria[j]
    value = sample_set$value[j]
    
    #if samples_included is specified, we would read in the samples file
    #we assume that sample names are unique and hence we are expecting a single column 
    #file without header
    #TODO: add exclude option
    if (varname == "samples") { 
      samples_filename <- value
      samples <- read_delim(samples_filename, delim=" ", col_names = FALSE) %>%
            pull (X1)
      covariates_current_set <- covariates_current_set %>% filter (sampleID %in% samples_to_include)
      write_lines(paste0("samples from ", samples_filename, " are included. There are ", nrow(covariates_current_set), " samples left. "), covar_set_log, append = TRUE)
    } else if (criteria == "include") {
      covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] %in% value),]
      write_lines(paste0("samples with ", varname, "==", value, " are included. There are ", nrow(covariates_current_set), " samples left."), covar_set_log, append = TRUE)
      
    } else if (criteria == ">") { # a bit stupid to separate them but not sure how to do this elegantly
      
      covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] > value),]
      write_lines(paste0("samples with ", varname, ">", value, " are included. There are ", nrow(covariates_current_set), " samples left."), covar_set_log, append = TRUE)
      
    } else if (criteria == "<") {
      
      covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] < value),]
      write_lines(paste0("samples with ", varname, "<", value, " are included. There are ", nrow(covariates_current_set), " samples left."), covar_set_log, append = TRUE)
      
    } else if (criteria == "==") {
      
      covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] == value),]
      write_lines(paste0("samples with ", varname, "==", value, " are included. There are ", nrow(covariates_current_set), " samples left."), covar_set_log, append = TRUE)
      
    } else {
      stop (paste0 ("Criteria ", criteria, " is not recognized."))
    }
    
    #output the file
    outfile_name <- paste0("covars_", sample_set_name, ".tsv")
    write.table(covariates_current_set, outfile_name, col.names = T,
                row.names = F, quote=FALSE, sep="\t")
  }
  
      
}


#variable info file parse
#phenoCol
binary_covar_list = ""
continuous_covar_list = ""

for (i in 1:nrow(variable_info)) {
  
  if (variable_info$excluded[i]!="no") {next} else if (variable_info$variableName[i]==phenoCol) {next}
  
  if (toupper(variable_info$variableType[i])=="SAMPLEID") {
    next
  } else if (toupper(variable_info$variableType[i]) == "BINARY") {
    if (binary_covar_list=="") {
      binary_covar_list = variable_info$variableName[i]
    }
    else {binary_covar_list = paste0(binary_covar_list, ",", variable_info$variableName[i])}
  } else if (variable_info$variableType[i] == "continuous") {
    if (continuous_covar_list ==  "") {continuous_covar_list = variable_info$variableName[i]}
    else {continuous_covar_list = paste0(continuous_covar_list, ",", variable_info$variableName[i])}  
  } else { #implementing this as hard exit instead of try catch for now
    stop(paste0("variable Type ", variable_info$variableType[i], " is not known, exiting..."))
  }  
}


#write out the covar lists
write_lines(binary_covar_list, "binary_covar_list.txt")
write_lines(continuous_covar_list, "continuous_covar_list.txt")



