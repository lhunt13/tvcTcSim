# Parallelize the procedure of getting the truth.
# Simulate datasets of size 5000, 1000 times to get
# an effective sample size of 5,000,000 from which to 
# estimate the truth.
library(getopt)
library(data.table)
source('parVals.R')
source('getTruth.R')

# set seed for run SGE_TASK_ID
# grab value of SGE_TASK_ID 
# (which is given by "-t" in the corresponding bash script)
args <- commandArgs(trailingOnly = TRUE)
boot.index <- as.numeric(args[1])

# set initial seed for reproducibility 
set.seed(123)
boot.seed <- sample(1e6, size = 3, replace = F)[boot.index]
set.seed(boot.seed)
get_truth(500,15,30)

# store results
# save file in "truth" directory with file name "run-<boot.index>.rds"
results.file <- file.path("truth", paste0("run-", boot.index, ".rds"))
saveRDS(truth, results.file)
# quit R and don't save workspace
quit('no')

