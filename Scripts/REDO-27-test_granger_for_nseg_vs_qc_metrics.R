rm(list = ls())
#dir = "/users/sana/Documents/CANMOD/7-RD-grangertests/7-freyjaxgranger"
dir = "/Users/sana/Documents/CANMOD/17-updated_vcf_files_2025/9-nseg_vs_qc_metrics/"

## load all data
everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv")

outcomes = c("nseg_relaxed")
predictors = c("Number.of.Reads", 
              "Breadth.of.coverage.Genome", 
              "Breadth.of.coverage.Spike", 
              "Mean.depth.Genome", 
              "Mean.depth.ORFs",
              "Mean.depth.Amplicon", 
              "Percentage.of.SARS.CoV.2.reads", 
              "Number.of.Amplicons.mean.depth.greater.100x.", 
              "Number.of.ORFs.mean.depth.greater.100x")


### granger causality test
library(lmtest)
library(tseries)
library(ggplot2)
max_order = 28

for(outcome in outcomes){ #outcome loop starts here
  
  print(outcome)
  vars = c(predictors, outcome)  
  
  ## remove any NAs in these variables (adf test cares about this)
  for(var in vars){
    everything[, var] = as.numeric(everything[, var])
    everything = everything[which(!is.na(everything[, var])), ]
  }
  
  
  
  #creates sites df
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
  
  ######################################################## section 1: only diffrentiate non-stationary series #########################################
  #setwd("/users/sana/Documents/CANMOD/7-RD-grangertests/6-qPCR/test/")
  for(s in 1:nrow(sites_df)){
    ## get city data
    city_id = sites_df$city_id[s]
    site_number = sites_df$site_number[s]
    
    city_data = everything[which(everything$city_id == city_id), ]
    city_data = city_data[which(city_data$site_number == site_number), ]
    city_data$date = as.Date(city_data$date)
    city_data = city_data[order(city_data$date), ]
    
    # ### convert character columns to numeric
    # city_data$tajimasd = as.numeric(city_data$tajimasd)
    # #city_data$filtered_shannon = as.numeric(city_data$filtered_shannon)
    # 
    # ## remove duplicated dates 
    # if(length(which(duplicated(city_data$date))) > 0){
    #   city_data = city_data[-which(duplicated(city_data$date)), ]
    # }
    # 
    # ##remove rows that tajima'sd and. shannon are NA
    # city_data = city_data[which(!is.na(city_data$tajimasd)), ]
    # #city_data = city_data[which(!is.na(city_data$filtered_shannon)), ]
    
    
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
      
      if(t$p.value >= 0.05){
        print(s)
        print("non stat")
      }
      # ## take second gradient
      # grad2 = diff(grad)
      # t = adf.test(grad2)
      # if(t$p.value < 0.05){
      #   stationarity_df[(3:nrow(stationarity_df)), var] = grad2
      #   next
      # }
      
    }
    
    write.table(stationarity_df, 
                file = paste0(dir,"/",city_id, "_", site_number, "_stationary_df.tsv"),
                sep = "\t",
                row.names = FALSE)
    print("you can ignore the warnings!")
  }
  
  
  
  ## remove postest from vars to avoid errors from granger:
  vars = vars[-which(vars == outcome)]
  
  
  ## dataframe to store all the granger test results
  all_sites_granger_df = data.frame(parameter = vector(),
                                    order = vector(),
                                    forward_causality = vector(),
                                    forward_pvalue = vector(),
                                    reverse_causality = vector(),
                                    reverse_pvalue = vector(),
                                    city_id = vector(),
                                    site_number = vector(),
                                    site_identifier = vector())
  
  
  
  
  for(s in 1:nrow(sites_df)){
    print(s)
    city_id = sites_df$city_id[s]
    site_number = sites_df$site_number[s]
    
    city_data = read.delim(paste0(dir,"/", city_id, "_", site_number, "_stationary_df.tsv"))
    city_data = city_data[-(1:2), ]
    
    for(i in 1:max_order){
      
      for(var in vars){ # var loop ends here
        
        
        #forward causality: 
        fmla_forward = as.formula(paste0(outcome, "~", var))  
        gtest_forward = grangertest(fmla_forward, order = i, data = city_data)
        
        #backward causality:
        fmla_backward = as.formula(paste0(var, "~", outcome))
        gtest_backward = grangertest(fmla_backward, order = i, data = city_data)
        
        ## add the values to the global dataframe
        all_sites_granger_df[nrow(all_sites_granger_df) + 1,] = data.frame(paste0(outcome, " vs ", var),
                                                                           i,
                                                                           paste0("y=",outcome, " x=", var),
                                                                           gtest_forward$`Pr(>F)`[2],
                                                                           paste0("y=", var, " x=",outcome),
                                                                           gtest_backward$`Pr(>F)`[2],
                                                                           city_id,
                                                                           site_number,
                                                                           paste0(city_id, "_0", site_number))
      } # var loop ends here
      
    }
    
    
  } ## loop ends here
  
  write.table(all_sites_granger_df,
              file = paste0(dir,"/all_sites_granger_df_outcome_", outcome, ".tsv"),
              sep = "\t",
              row.names = FALSE)
  
  
  ###fill out the counts
  counts_long_df = data.frame(order = vector(),
                              count = vector(),
                              arrangement = vector())
  
  forwards = unique(all_sites_granger_df$forward_causality)
  reverses = unique(all_sites_granger_df$reverse_causality)
  
  
  for(o in 1:max_order){
    for(c in forwards){
      sub = all_sites_granger_df[which(all_sites_granger_df$forward_causality == c), ]
      sub = sub[which(sub$order == o), ]
      counts_long_df[nrow(counts_long_df) + 1, ] = data.frame(order = o,
                                                              count = length(which(sub$forward_pvalue < 0.05)),
                                                              arrangement = c)
    }
    for(c in reverses){
      sub = all_sites_granger_df[which(all_sites_granger_df$reverse_causality == c), ]
      sub = sub[which(sub$order == o), ]
      counts_long_df[nrow(counts_long_df) + 1, ] = data.frame(order = o,
                                                              count = length(which(sub$reverse_pvalue < 0.05)),
                                                              arrangement = c)
    }
    
    
  }
  
  write.table(counts_long_df, 
              file = paste0(dir,"/counts_long_df_outcome_", outcome, ".tsv"),
              sep = "\t",
              row.names = FALSE)
  
  #counts_long_df$order = -1 * counts_long_df$order
  
  p = ggplot(data = counts_long_df) + 
    geom_tile(aes(x = order, y = arrangement, fill = count, group = arrangement), color = "black", linewidth = 0.95) +
    scale_fill_gradient2(low = "white",
                         high = "#752F45",
                         limits = c(0,10)) +
    geom_text(aes(x=order, y=arrangement, label = count), size=7) +
    theme_bw()+
    theme(text = element_text(size= 50))
  
  p
  ggsave(p, file = paste0(dir, "/granger_causality_outcome_",outcome, ".pdf"),
         width = 50,
         height = 50,
         limitsize = FALSE)
} #outcome loop ends here

