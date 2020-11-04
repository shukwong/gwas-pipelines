library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 ) {
  stop("Arguments to process_metal_output.R: \
     metal_output_file" )
}

metal_filename <- args[1]

metal_output <- fread(metal_filename)

markers <- str_split_fixed(metal_output$MarkerName, ":", 4) %>% 
           as.data.frame() 
colnames(markers) <- c("CHR", "POS", "REF", "ALT")     


metal_output <- cbind (markers[,1:2], metal_output)


#write output
metal_prefix <- basename(metal_filename)

write_tsv(metal_output, paste0("./", metal_prefix))
