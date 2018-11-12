# merge all files ina folcer into one big .rds file
library(getopt)
args <- commandArgs(trailingOnly = TRUE)

# `folder` is a character denoting the folder in which to merge all files
# these are "truth/" and "results/"
merge_data <- function(folder){
  do.call('rbind', lapply(list.files(folder, full.names = TRUE), readRDS))
}

merged <- merge_data(args)

saveRDS(paste0(folder,"merged.RDS"))
  