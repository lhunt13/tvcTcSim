# merge .rds files to one big .rds file

mergedat <- do.call('rbind', lapply(list.files("truth/", full.names = TRUE), readRDS))
grand_ave <- mean(mergedat)

saveRDS(grand_ave,"truth/grand_average.rds")