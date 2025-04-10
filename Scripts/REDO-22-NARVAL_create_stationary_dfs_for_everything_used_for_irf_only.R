### accompanying script: creates stationary dataframes on NARVAL

### this scripts writes stationary dataframes
### only non stationary time series are differentiated
### NO SMOOTHING is done
rm(list = ls())

library(vars)
library(tseries)
library(lmtest)
library(ggplot2)
library(pracma)

sites_df = data.frame(city_id = vector(),
                      site_number = vector())

sites_df[nrow(sites_df) + 1, ] = c("MTL", 1)
sites_df[nrow(sites_df) + 1, ] = c("MTL", 2)
sites_df[nrow(sites_df) + 1, ] = c("LVL", 1)
sites_df[nrow(sites_df) + 1, ] = c("LVL", 2)
sites_df[nrow(sites_df) + 1, ] = c("LVL", 5)
sites_df[nrow(sites_df) + 1, ] = c("GTN", 1)
sites_df[nrow(sites_df) + 1, ] = c("QC", 1)
sites_df[nrow(sites_df) + 1, ] = c("QC", 2)
sites_df[nrow(sites_df) + 1, ] = c("Riki", 1)
sites_df[nrow(sites_df) + 1, ] = c("Sept", 1)

### setdir
dir = "/home/sanna/projects/ctb-shapiro/sanna/sanna/CANMOD/19-REDO/4-irf_varselect"


### load everything and add columns of interest
everything = read.delim("Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv")


predictors = c("Number.of.Reads", 
               "Breadth.of.coverage.Genome", 
               "Breadth.of.coverage.Spike", 
               "Mean.depth.Genome", 
               "Mean.depth.ORFs",
               "Mean.depth.Amplicon", 
               "Percentage.of.SARS.CoV.2.reads", 
               "Number.of.Amplicons.mean.depth.greater.100x.", 
               "Number.of.ORFs.mean.depth.greater.100x",
               "nseg_relaxed",
               "vcopy_per_100k")

#predictors = c("Number.of.signficant.SNVs")
outcomes = c("daily_count", "daily_rate", "tests_taken")



write_stationary_dataframe = function(s, outcome){
  vars = c(predictors, outcome)
  ## get city data
  city_id = sites_df$city_id[s]
  site_number = sites_df$site_number[s]
  
  city_data = everything[which(everything$city_id == city_id), ]
  city_data = city_data[which(city_data$site_number == site_number), ]
  city_data$date = as.Date(city_data$date)
  city_data = city_data[order(city_data$date), ]
  
  ###### check for stationarity in the time series
  stationarity_df = data.frame(matrix(NA, nrow = nrow(city_data), ncol = (length(vars) + 3)))
  colnames(stationarity_df) = c("city_id", "site_number", "date", vars)
  
  stationarity_df$date = city_data$date
  stationarity_df$city_id = city_id
  stationarity_df$site_number = site_number
  
  for(var in vars){
    t = adf.test(city_data[, var])
    if(t$p.value < 0.05){
      stationarity_df[, var] = city_data[, var]
      next
    }
    
    ## take first gradient
    grad = diff(city_data[, var])
    t = adf.test(grad)
    if(t$p.value < 0.05){
      stationarity_df[(2:nrow(stationarity_df)), var] = grad
      next
    }
    
    ## take second gradient
    grad2 = diff(grad)
    t = adf.test(grad2)
    if(t$p.value < 0.05){
      stationarity_df[(3:nrow(stationarity_df)), var] = grad2
      next
    }
    
  }
  
  
  #system(paste0("mkdir ", dir,"stationary_dfs_seed_", seed, "_iterations_", its))
  
  write.table(stationarity_df, 
              file = paste0(dir, "/",city_id, "_", site_number,"_outcome_",outcome, "_stationary_df.tsv"),
              sep = "\t",
              row.names = FALSE)
  print("you can ignore the warnings!")
}

outcome = "daily_count"
lapply(1:10, write_stationary_dataframe, outcome = outcome)

outcome = "daily_rate"
lapply(1:10, write_stationary_dataframe, outcome = outcome)

outcome = "tests_taken"
lapply(1:10, write_stationary_dataframe, outcome = outcome)

