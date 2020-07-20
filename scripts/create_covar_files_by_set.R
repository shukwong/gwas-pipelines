require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 ) {
  stop("Arguments to combine_covars.R: \
     create_covar_files_by_set.R variable-info_tsv_file sample_sets_json_file" )
}

variable-info_tsv_file <- args[1]
sample_sets_json_file <- args[2]

