# analyze data on cluster
library(getopt)
library(data.table)
library(splines)
library(rms)
source('parVals.R')
source('simulateData.R')
source('analyzeData.R')

# define some global variables
samp_size <- 5000
boot_num <- 100
tval <- 1000
bandwidth <- 365

# set seed for run SGE_TASK_ID
# grab value of SGE_TASK_ID 
# (which is given by "-t" in the corresponding bash script)
args <- commandArgs(trailingOnly = TRUE)
boot.index <- as.numeric(args[1])

# set initial seed for reproducibility 
set.seed(123)
boot.seed <- sample(1e6, size = tval, replace = F)[boot.index]
set.seed(boot.seed)

# simulate data
data <- sim_obs(samp_size,ssU,pR1,u_star)

# perform analysis
rmdiff <- analyze(DATA=data,BAND=bandwidth,NUMSIM=samp_size*4)

# perform bootstrap
ci <- bootstrap(data,boot_num,BAND=bandwidth,NUMSIM=samp_size*4)

results <- c(rmdiff,ci)

# store results
# save file in "truth" directory with file name "run-<boot.index>.rds"
results.file <- file.path("results", paste0("run-", boot.index, ".rds"))
saveRDS(results, results.file)
# quit R and don't save workspace
quit('no')






