#This script takes in the files specified by the user and creates different sets of covar files for different runs

#variable info is defined in variable_info_tsv_file

#TODO: add transformation for quantitative traits

libraries_required <- c("tidyverse","data.table","jsonlite","fastDummies")

for (libs in libraries_required) {
  if( !require(libs, character.only = T)) {
    install.packages( libs,  repos = c(CRAN = "https://cloud.r-project.org/") )
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 ) {
    stop("Arguments to create_covar_files_by_set.R: \
          create_covar_files_by_set.R covariate_tsv_file variable-info_tsv_file sample_sets_json_file" )
}

  
covariate_tsv_file <- args[1]
variable_info_tsv_file <- args[2]
sample_sets_json_file <- args[3]

covariates <- fread(covariate_tsv_file)
variable_info <- fread(variable_info_tsv_file)
sample_sets <- fromJSON(sample_sets_json_file)


#variable info file parse
#phenoCol
binary_covar_string = ""
continuous_covar_string = ""
binary_covar_list = vector()
continuous_covar_list = vector()

for (i in 1:nrow(variable_info)) {
  
  if (variable_info$variableType[i]=="phenotype_quantitative") { 
    write_lines(variable_info$variableName[i], "phenotype_line.txt")
    write_lines("quantitative", "phenotype_type.txt")
    phenoCol = variable_info$variableName[i]
    next
  } else if (variable_info$variableType[i]=="phenotype_binary") {
    write_lines(variable_info$variableName[i], "phenotype_line.txt")
    write_lines("binary", "phenotype_type.txt")
    phenoCol = variable_info$variableName[i]
    next
  }   else if (toupper(variable_info$variableType[i])=="SAMPLEID") {
    write_lines(variable_info$variableName[i], "sampleid_line.txt")
    sampleID = variable_info$variableName[i]
    next
  } else if (toupper(variable_info$variableType[i]) == "BINARY") {
    binary_covar_list <- c(binary_covar_list, variable_info$variableName[i])
    if (binary_covar_string=="") {
       binary_covar_string = variable_info$variableName[i]
    }
    else {
      binary_covar_string = paste0(binary_covar_string, ",", variable_info$variableName[i])
    }
  } else if (variable_info$variableType[i] == "quantitative") {
    continuous_covar_list = c(continuous_covar_list, variable_info$variableName[i])
    if (continuous_covar_string ==  "") {continuous_covar_string = variable_info$variableName[i]}
    else {
      continuous_covar_string = paste0(continuous_covar_string, ",", variable_info$variableName[i])
    }  
  } else if (variable_info$variableType[i] == "categorical") {
    variable_name = variable_info$variableName[i]
    new_variable_name = paste0(variable_name, "_dummy")
    covariates <- covariates %>% 
                  mutate (!!new_variable_name := !!as.name(variable_name)) %>% 
                  dummy_cols(select_columns = new_variable_name , remove_first_dummy = TRUE, remove_selected_columns=TRUE) 
   binary_covar_list <- c(binary_covar_list, colnames(covariates)[grep("center_dummy" ,colnames(covariates))])  
   if (binary_covar_string=="") {
      binary_covar_string = paste(colnames(covariates)[grep("center_dummy" ,colnames(covariates))], collapse = ',', sep='')
    }
    else {
      binary_covar_string = paste0(binary_covar_string, ",", paste(colnames(covariates)[grep("center_dummy" ,colnames(covariates))], collapse = ',', sep=''))
    }
  }
  else { #implementing this as hard exit instead of try catch for now
    stop(paste0("variable Type ", variable_info$variableType[i], " is not known, exiting..."))
  }  
}


#write out the covar lists
write_lines(binary_covar_string, "binary_covar_list.txt")
write_lines(continuous_covar_string, "continuous_covar_list.txt")

#process sample sets
for (i in 1:length(sample_sets)) {
  if (names(sample_sets[1]) == "phenoCol") {next}
  
  sample_set_name = names(sample_sets)[i]
  sample_set <- sample_sets[[i]]
  covariates_current_set = covariates
  
  #start the log file
  covar_set_log <- paste0(sample_set_name, ".log")
  write_lines(paste0("Full covar file has ", nrow(covariates_current_set), " samples."),      
              covar_set_log)
  
  if (sample_set_name == "complete_set") {
    
    write_lines("Full set of samples are includes in current set", covar_set_log, append = TRUE)
    
  } else {
  
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
      if (varname == "samples") {  #this does not work yet, as cromwell needs to localize the sample file
        samples_filename <- value
        samples <- read_delim(samples_filename, delim=" ", col_names = FALSE) %>%
            pull (X1)
        covariates_current_set <- covariates_current_set %>% filter (sampleID %in% samples_to_include)
        write_lines(paste0("samples from ", samples_filename, " are included. There are ", 
                           nrow(covariates_current_set), " samples left. "), covar_set_log, append = TRUE)
      } else if (criteria == "include") {
        covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] %in% value),]
        write_lines(paste0("samples with ", varname, "==", value, " are included. There are ", 
                           nrow(covariates_current_set), " samples left."), covar_set_log, append = TRUE)
      
      } else if (criteria == ">") { # a bit stupid to separate them but not sure how to do this elegantly
      
        covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] > value),]
        write_lines(paste0("samples with ", varname, ">", value, " are included. There are ", 
                           nrow(covariates_current_set), " samples left."), covar_set_log, append = TRUE)
      
      } else if (criteria == "<") {
      
        covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] < value),]
        write_lines(paste0("samples with ", varname, "<", value, " are included. There are ", 
                           nrow(covariates_current_set), " samples left."), covar_set_log, append = TRUE)
      
      } else if (criteria == "==") {
      
        covariates_current_set <- covariates_current_set[which(covariates_current_set[[varname]] == value),]
        write_lines(paste0("samples with ", varname, "==", value, " are included. There are ", 
                           nrow(covariates_current_set), " samples left."), covar_set_log, append = TRUE)
      
      } else {
        stop (paste0 ("Criteria ", criteria, " is not recognized."))
      }
    
    }
  }  
  
  #get complete cases
  covariates_current_set <- covariates_current_set %>% select(c(!!sym(sampleID), !!sym(phenoCol), all_of(binary_covar_list), all_of(continuous_covar_list)))
  covariates_current_set <- covariates_current_set[complete.cases(covariates_current_set),]
  
  #output the file
  outfile_name <- paste0(sample_set_name, "_covars", ".tsv")
  write.table(covariates_current_set, outfile_name, col.names = T,
              row.names = F, quote=FALSE, sep="\t")
      
}





