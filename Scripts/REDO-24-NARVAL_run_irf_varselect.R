######## SCRIPT FOR NARVAL

### this script runs the iRF with different seeds and different iteration counts to. make sure the outcomes are stable
### the grid determines the values that are going to be run through

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

grid_index = as.numeric(args[1])

##### script to run IRF on narval

library(vars)
library(tseries)
library(lmtest)
library(ggplot2)
library(pracma)
library(stringr)


### defines sites df
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

### define grid of values
grid = expand.grid(seed = c(1), 
                   its = c(10000), 
                   outcome = c("daily_count", "daily_rate", "tests_taken"),
                   predictor = c("Number.of.Reads", 
                                 "Breadth.of.coverage.Genome", 
                                 "Breadth.of.coverage.Spike", 
                                 "Mean.depth.Genome", 
                                 "Mean.depth.ORFs",
                                 "Mean.depth.Amplicon", 
                                 "Percentage.of.SARS.CoV.2.reads", 
                                 "Number.of.Amplicons.mean.depth.greater.100x.", 
                                 "Number.of.ORFs.mean.depth.greater.100x",
                                 "nseg_relaxed",
                                 "vcopy_per_100k"))




### setdir
dir = "/home/sanna/projects/ctb-shapiro/sanna/sanna/CANMOD/19-REDO/4-irf_varselect/"

### load everything and add columns of interest
everything = read.delim("/lustre06/project/6032406/sanna/sanna/CANMOD/19-REDO/4-irf_varselect/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv")
varselect_grid = read.delim("/lustre06/project/6032406/sanna/sanna/CANMOD/19-REDO/4-irf_varselect/varselect_grid_all_preds_and_outcomes.tsv")
n_ahead = 28


get_irf_for_one_s = function(s,outcome, predictor, cumulative_bool){
  # print(s)
  # print(u)
  city_id = sites_df$city_id[s]
  site_number = sites_df$site_number[s]
  
  site_identifier = paste0(city_id, "_", site_number)
  
  stationarity_df = read.delim(paste0(dir, city_id, "_", site_number,"_outcome_",outcome, "_stationary_df.tsv"))
  
  ## remove first row
  city_data = stationarity_df[-1,]
  
  
  ### create model for nseg and daily count only
  data = city_data[, c(predictor, outcome)]
  data = sapply(data, scale)
  
  ### get the best model order from the grid previously created
  u <<- varselect_grid[(varselect_grid$site == site_identifier) & (varselect_grid$predictor == predictor) & (varselect_grid$outcome == outcome), 
                       "best_model_order"]
  
  
  model = VAR(data, p=u, type = "const")
  feir = irf(model, n.ahead = n_ahead, ortho = TRUE, runs = its, response = outcome, impulse = predictor, cumulative = cumulative_bool, seed = seed)
  
  
  out = data.frame(time_point = 0:n_ahead,
                         city_id = city_id,
                         site_number = site_number,
                         #permutation_number = "unshuffled data",
                         irf_vals = as.vector(feir$irf[[predictor]]),
                         upper_irf = as.vector(feir$Upper[[predictor]]),
                         lower_irf = as.vector(feir$Lower[[predictor]]),
                         shock_applied_to = predictor,
                         response = outcome,
                         model_type = "2 var model",
                         model_order = as.integer(u),
                         diff_flag = NA)
  
  out$sig = NA
  out$sig[out$lower_irf > 0] = out$irf_vals[out$lower_irf > 0]
  out$sig[out$upper_irf < 0] = out$irf_vals[out$upper_irf < 0]
  
  p = ggplot(data = out[out$shock_applied_to == predictor, ]) + 
    geom_line(aes(x=time_point, y=irf_vals), color = "blue") +
    geom_ribbon(aes(ymin=lower_irf, ymax=upper_irf, x=time_point), alpha = 0.3,  fill = "blue")+
    #geom_line(aes(x=time_point, y=upper_irf), color = "blue", linetype = 2)+
    #geom_line(aes(x=time_point, y=lower_irf), color = "blue", linetype = 2)+
    geom_point(aes(x=time_point, y=sig), shape = 8, color = "blue", size = 1, stroke = 1)+
    theme_bw()+
    ylim(c(-1.5,1.5))+
    ylab(str_wrap(paste0("Cumulative response in ", outcome), width = 20)) +
    xlab("n time points after shock")+
    ggtitle(paste0(sites_df$city_id[s], "-", sites_df$site_number[s], " - ", predictor))

  
  response_type = "cumulative"
  if(!cumulative_bool){
    response_type = "non_cumulative"
  }
  
  out$response_type = response_type
  out$seed = seed
  out$iterations = its
  
  ggsave(p, file = paste0(dir, "plots/",response_type, sites_df$city_id[s], "_", sites_df$site_number[s], "_2varmodel_", outcome, "_model_varselect","_seed_",seed, "_iterations_", its, ".pdf"),
         device = "pdf")
  # write.table(out,
  #             file = paste0(dir, response_type, sites_df$city_id[s], "_", sites_df$site_number[s], "_2varmodel_", outcome, "_model_varselect", "_seed_",seed, "_iterations_", its, ".tsv"),
  #             sep = "\t",
  #             row.names = FALSE)
  return(out)
  
}

outcome = as.character(grid$outcome[grid_index])
predictor = as.character(grid$predictor[grid_index])
seed = grid$seed[grid_index]
its = grid$its[grid_index]
# s = grid$s[grid_index]

all_s_out = do.call(rbind, lapply(1:10, get_irf_for_one_s, outcome = outcome, predictor = predictor, cumulative_bool = TRUE))
response_type = unique(all_s_out$response_type)

write.table(all_s_out,
            file = paste0(dir, response_type, "_all_sites", "_2varmodel_outcome_", outcome,"_pred_",predictor, "_model_varselect_seed_", seed, "_iterations_", its, ".tsv"),
            sep = "\t",
            row.names = FALSE)



