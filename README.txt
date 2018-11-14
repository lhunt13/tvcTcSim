Simulate longitudinal survival data that exibits time varying confounding and temporal confounding. Use to perform simulation study of an analysis that adjusts for both of these.

First run analyzeDataParalBash. To change the number of simulations, open the bash script and change the line `#$ -t 1-1000`, as well as the last argument passed to the Rscript. To change the sample size, number of bootstrap samples, or bandwidth (measured in days) then change arguments 2-4, respectively. 

Second, run mergeBash. This merges all the output into one .rds file so you can examine coverage, bias, etc... In the bash script, the arguments passed to the Rscript are (1) the folder in which files will be merged, (2) the pattern determining all the files to be merged, (3) the name of the merged .rds file. 

Third, run getTruthParalBash. This simulates counterfactual data in order to get the true value.

Fourth, change the arguments in mergeBash in order to merge them into one file to compute the average. The average value will be the truth.