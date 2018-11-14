# rscript to merge all .rds files to a chosen output 
# in a given folder with a given pattern
main <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  folder <- args[1]
  pattern <- args[2]
  output <- args[3]
  
  merged <- do.call('rbind', 
                    lapply(list.files(path = folder,
                                      pattern = pattern, 
                                      full.names = TRUE), readRDS))
  
  saveRDS(merged,paste0(folder,output))
  
}

main()


  




  