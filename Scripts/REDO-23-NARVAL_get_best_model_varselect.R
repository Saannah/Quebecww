### define sites df
rm(list = ls())
library(vars)
library(tseries)
library(lmtest)
library(ggplot2)
library(pracma)

### setdir
dir = "/home/sanna/projects/ctb-shapiro/sanna/sanna/CANMOD/19-REDO/4-irf_varselect/"


### load everything
everything = read.delim("Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv")

### define predictors and outcome
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

outcomes = c("daily_count", "daily_rate", "tests_taken")


### determine maximum order of models to be checked
maxorder=15


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


##### this has been removed because staionary df has been already written by script 22

# write_stationary_dataframe = function(s, outcome){
#   ## get city data
#   
#   vars = c(predictors, outcome)
#   city_id = sites_df$city_id[s]
#   site_number = sites_df$site_number[s]
#   
#   city_data = everything[which(everything$city_id == city_id), ]
#   city_data = city_data[which(city_data$site_number == site_number), ]
#   city_data$date = as.Date(city_data$date)
#   city_data = city_data[order(city_data$date), ]
#   
#   ###### check for stationarity in the time series
#   stationarity_df = data.frame(matrix(NA, nrow = nrow(city_data), ncol = (length(vars) + 3)))
#   colnames(stationarity_df) = c("city_id", "site_number", "date", vars)
#   
#   stationarity_df$date = city_data$date
#   stationarity_df$city_id = city_id
#   stationarity_df$site_number = site_number
#   
#   for(var in vars){
#     t = adf.test(city_data[, var])
#     if(t$p.value < 0.05){
#       stationarity_df[, var] = city_data[, var]
#       next
#     }
#     
#     ## take first gradient
#     grad = diff(city_data[, var])
#     t = adf.test(grad)
#     if(t$p.value < 0.05){
#       stationarity_df[(2:nrow(stationarity_df)), var] = grad
#       next
#     }
#     
#     ## take second gradient
#     grad2 = diff(grad)
#     t = adf.test(grad2)
#     if(t$p.value < 0.05){
#       stationarity_df[(3:nrow(stationarity_df)), var] = grad2
#       next
#     }
#     
#   }
#   
#   
#   
#   
#   write.table(stationarity_df, 
#               file = paste0(dir,city_id, "_", site_number,"_outcome_",outcome, "_stationary_df.tsv"),
#               sep = "\t",
#               row.names = FALSE)
#   print("you can ignore the warnings!")
# }
# 


get_varselect_grid = function(s, outcome, predictor){
  
  city_id = sites_df$city_id[s]
  site_number = sites_df$site_number[s]
  site_identifier = paste0(city_id, "_", site_number)
  
  stationarity_df = read.delim(paste0(dir,city_id, "_", site_number,"_outcome_",outcome, "_stationary_df.tsv"))
  
  ## remove first row
  city_data = stationarity_df[-1,]
  
  
  ### create model for nseg and daily count only
  data = city_data[, c(predictor, outcome)]
  data = sapply(data, scale)
  model_varselect = VARselect(data, lag.max = maxorder)
  
  ## take best model order
  u <<- model_varselect$selection[1]
  model = VAR(data, p=u, type = "const")
  
  t1 =cbind(
    data.frame(site = site_identifier, 
               predictor = predictor,
               outcome = outcome,
               best_model_order = model_varselect$selection[1],
               best_model_roots = paste0(vars::roots(model), collapse = "/"),
               best_roots_outsite_unit_circle = length(which(abs(vars::roots(model)) >= 1))),
    matrix(data = model_varselect$criteria["AIC(n)", ],
           nrow = 1)
  )
  colnames(t1)[7:length(colnames(t1))] = paste0("p=", colnames(t1)[7:length(colnames(t1))])
  
  return(t1)
}


combos_grid = expand.grid(outcomes, predictors)
varselect_grid = data.frame()

for(i in 1:nrow(combos_grid)){
  outcome <- as.character(combos_grid[i, 1])
  predictor <- as.character(combos_grid[i, 2])
  varselect_subgrid = do.call(rbind,
                           lapply(1:nrow(sites_df), get_varselect_grid, outcome = outcome, predictor = predictor))
  varselect_grid = rbind(varselect_grid, varselect_subgrid)
  
}

write.table(varselect_grid,
            file = paste0(dir, "varselect_grid_all_preds_and_outcomes.tsv"),
            sep = "\t",
            row.names = FALSE)

