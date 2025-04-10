#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tseries)
library(lmtest)


s = as.numeric(args[1])
max_order = 28



### load data and create useful columns
everything = read.delim("/lustre06/project/6032406/sanna/sanna/CANMOD/19-REDO/1-grangertest/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv")

#everything$week_total_positive_count


predictors = c("nseg_relaxed", "vcopy_per_person", "vcopy_per_100k")
outcomes = c("daily_count", "daily_rate", "tests_taken")
vars = c(predictors, outcomes)

## remove any NAs in these variables (adf test cares about this)
for(var in vars){
  everything[, var] = as.numeric(everything[, var])
  everything = everything[which(!is.na(everything[, var])), ]
}



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



## check stationarity for outcome 
for(outcome in outcomes){
  t = adf.test(city_data[, outcome])
  if(t$p.value >= 0.05){
    city_data[2:nrow(city_data), outcome] = diff(city_data[, outcome])
    city_data[1, outcome] = NA
    #city_data = city_data[-1, ]
    t = adf.test(city_data[2:nrow(city_data), outcome])
    if(t$p.value >= 0.05){
      stop("time series is not stationary after 1 gradient")
    }
  }
  
}

for(outcome in outcomes){ 
  
  if(is.na(city_data[1, outcome])){
    city_data = city_data[-1, ]
  }
}


for(outcome in outcomes){ #outcome loop starts here
  print(outcome)
  for(pred in predictors){ ##predictor loop starts here
    print(pred)
    counter = 0
    set_of_500_permutations = data.frame(city_id = vector(),
                                         site_number = vector(),
                                         permutation_number = vector(),
                                         diff_flag = vector(),
                                         order = vector(),
                                         causality = vector(),
                                         pvalue = vector())
    while(counter <= 500){
      print("counter:")
      print(counter)
      permutations = data.frame(city_id = vector(),
                                site_number = vector(),
                                permutation_number = vector(),
                                diff_flag = vector(),
                                order = vector(),
                                causality = vector(),
                                pvalue = vector())
      diff_flag = "no"
      
      
      ### shuffle predictor
      stationary_df = data.frame(matrix(0, nrow = nrow(city_data), ncol = length(vars)))
      colnames(stationary_df) = vars
      
      stationary_df[, pred] = city_data[sample(1:nrow(city_data), nrow(city_data)), pred]
      stationary_df[, outcome] = city_data[, outcome]
      
      ### test for stationarity
      t = adf.test(stationary_df[, pred])
      if(t$p.value >= 0.05){
        stationary_df[2:nrow(stationary_df), pred] = diff(stationary_df[, pred])
        #city_data$postest[1] = NA
        diff_flag = "yes"
        stationary_df = stationary_df[-1, ]
        t = adf.test(stationary_df[, pred])
        if(t$p.value >= 0.05){
          next
        }
      }
      
      
      for(k in 1:max_order){ #order loop starts here
        
        result <- tryCatch({
          fmla = as.formula(paste0(outcome, "~", pred))
          #forward causality: y = postest, x = nseg
          gtest_forward = grangertest(fmla, order = k, data = stationary_df)
          
        }, error = function(err) {
          # This block will be executed if an error occurs
          break
        })
        
        
        if(!is.null(result)) {
          permutations = rbind(permutations,
                               data.frame(city_id = city_id,
                                          site_number = site_number,
                                          permutation_number = counter,
                                          diff_flag = diff_flag,
                                          order = k,
                                          causality = paste0("y=",outcome, " x=", pred),
                                          pvalue = gtest_forward$`Pr(>F)`[2]))
        }
        
        
        
      } #order loop ends here
      
      
      if(nrow(permutations) == max_order){
        counter = counter + 1
        # write.table(permutations,
        #             file = paste0(city_id, "_", site_number, "_outcome_",outcome,"_pred_",pred,"_", counter, ".tsv"),
        #             sep = "\t",
        #             row.names = FALSE)
        set_of_500_permutations = rbind(set_of_500_permutations,
                                        permutations)
      }
      
    }
    
    write.table(set_of_500_permutations, 
                file = paste0(city_id, "_", site_number, "_outcome_", outcome,"_pred_",pred,".tsv" ),
                sep = "\t",
                row.names = FALSE)
    
    
  } ##predictor loop ends here
  
} ## outcome loop ends here



