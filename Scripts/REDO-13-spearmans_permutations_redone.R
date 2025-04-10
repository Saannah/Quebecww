#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
s = as.numeric(args[1])

#### spearman's correlation on slimhq
#dir = "/home/sanna/projects/ctb-shapiro/sanna/sanna/CANMOD/19-REDO/2-spearmans"


lag_duration = 28
spearman_results = data.frame()


### load data
everything = read.delim("/lustre06/project/6032406/sanna/sanna/CANMOD/19-REDO/1-grangertest/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv")

predictors = c("nseg_relaxed", "vcopy_per_person", "vcopy_per_100k")
outcomes = c("daily_count", "daily_rate", "tests_taken")

vars = c(outcomes, predictors)

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




#create a funcion thar calculates cross correlation for a k shift in both directions of two vectors:
cross_correlation = function(k, vec1, vec2){
  #shift vec1 forward:
  if(k >= 0){
    vec1_shifted = vec1[1:(length(vec1) - k)]
    vec2_shifted = vec2[(k + 1):length(vec2)]
    #v1_forward_cc = (sum(vec1_shifted * vec2_shifted))/sqrt(sum(vec1 * vec1) * sum(vec2 * vec2))
    test_forward = cor.test(vec1_shifted, vec2_shifted, method = "spearman", exact = FALSE)
    return(data.frame(direction = "v1_forward",
                      lag = k,
                      #cc = v1_forward_cc,
                      estimate = test_forward$estimate,
                      pvalue = test_forward$p.value))
  }
  
  
  if(k < 0){
    lag = k
    k = abs(k)
    #shift vec1 backward:
    vec1_shifted = vec1[(k + 1):length(vec1)]
    vec2_shifted = vec2[1:(length(vec2) - k)]
    #v1_backward_cc = sum(vec1_shifted * vec2_shifted)/sqrt(sum(vec1 * vec1) * sum(vec2 * vec2))
    test_backward = cor.test(vec1_shifted, vec2_shifted, method = "spearman", exact = FALSE)
    return(data.frame(direction = "v1_backward",
                      lag = lag,
                      #cc = v1_backward_cc,
                      estimate = test_backward$estimate,
                      pvalue = test_backward$p.value))
  }
  
}


### shuffling the time series n times for 1 site
city_id = sites_df$city_id[s]
site_number = sites_df$site_number[s]

city_data = everything[which(everything$city_id == city_id), ]
city_data = city_data[which(city_data$site_number == site_number), ]
city_data$date = as.Date(city_data$date)
city_data = city_data[order(city_data$date), ]

## remove duplicated dates 
if(length(which(duplicated(city_data$date))) > 0){
  city_data = city_data[-which(duplicated(city_data$date)), ]
}
city_data = city_data[, vars]




##### run permutations

for(outcome in outcomes){ #outcome loop starts here
  print(outcome)
  for(pred in predictors){ ##predictor loop starts here
    print(pred)
    counter = 0
    set_of_500_permutations = data.frame()
    
    while(counter <= 500){
      print("counter:")
      print(counter)
      
      ## shuffle predictor
      city_data[, pred] = city_data[sample(1:nrow(city_data), nrow(city_data)), pred]
      
      cross_corr = as.data.frame(do.call(rbind, (lapply(-lag_duration:lag_duration, cross_correlation, 
                                                        city_data[, pred], city_data[, outcome]))))
      cross_corr$city_id = city_id
      cross_corr$site_number = site_number
      cross_corr$comparison = paste0(pred, " vs ", outcome)
      
      ## append this to the set of permutations
      set_of_500_permutations = rbind(set_of_500_permutations,
                                      cross_corr)
      
      
      counter = counter + 1
      
    }
    
    write.table(set_of_500_permutations, 
                file = paste0("spearmans_permutations_",city_id, "_", site_number, "_outcome_", outcome,"_pred_",pred,".tsv" ),
                sep = "\t",
                row.names = FALSE)
    
    
  } ##predictor loop ends here
  
} ## outcome loop ends here








