### get summary statistics on the correlation heatmaps:


#################################### Spearman's correlation
rm(list = ls())
sites_df = data.frame(city_id = vector(),
                      site_number = vector())
lag_duration = 28

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
sites_df$site_identifier = paste0(sites_df$city_id, "_", sites_df$site_number)

predictors = c("nseg_relaxed", "vcopy_per_100k")

## write a function that creates the summary dataframe for a given outcome
get_summary = function(outcome){
  all = data.frame()
  for(s in 1:nrow(sites_df)){
    all = rbind(all,
                read.delim(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/4-spearmans_permutations/spearmans_permutations_", sites_df$site_identifier[s], "_outcome_", outcome, "_pred_nseg_relaxed.tsv")),
                read.delim(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/4-spearmans_permutations/spearmans_permutations_", sites_df$site_identifier[s], "_outcome_", outcome, "_pred_vcopy_per_100k.tsv"))
    )
  }
  
  summary = data.frame(
    lag = vector(),
    sig_count = vector(),
    total_count = vector(),
    sig_proportion = vector(),
    Count = vector(),
    test = vector(),
    perm = vector())
  
  for(l in (-lag_duration:lag_duration)){
    
    for(pred in predictors){
      perms = all[which(all$lag == l), ]
      perms = perms[which(perms$comparison == paste0(pred, " vs ", outcome)), ]
      sigs = perms[perms$estimate > 0 & perms$pvalue < 0.05, ]
      
      
      summary = rbind(summary,
                      data.frame(
                        lag = l,
                        sig_count = nrow(sigs),
                        total_count = nrow(perms),
                        sig_proportion = (nrow(sigs)/nrow(perms)),
                        Count = (nrow(sigs)/nrow(perms))*10,
                        test = paste0(pred, " vs ", outcome),
                        perm = "permutated data"))
    }
  }
  ## get real data and append
  real_data = read.delim(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/2-spearmans/spearmans_outcome_", outcome,".tsv"))
  
  summary = rbind(summary,
                  data.frame(lag = real_data$lag,
                             sig_count = NA,
                             total_count = NA,
                             sig_proportion = NA,
                             Count = real_data$count,
                             test = real_data$comparison,
                             perm = real_data$comparison)
                  
                  
  )
  
  return(summary)
}


summary_daily_count = get_summary("daily_count")
daily_count_nseg = summary_daily_count[str_detect(summary_daily_count$perm, "nseg"), ]
daily_count_vcopy = summary_daily_count[str_detect(summary_daily_count$perm, "vcopy_per_100k"), ]

median(daily_count_nseg$Count[daily_count_nseg$lag %in% 0:7])
median(daily_count_vcopy$Count[daily_count_vcopy$lag %in% 0:7])


median(daily_count_nseg$Count[daily_count_nseg$lag %in% -7:0])
median(daily_count_vcopy$Count[daily_count_vcopy$lag %in% -7:0])



###################################### Compare strength of correlation ######################################
### create grid of correlation values
rm(list = ls())

library(lme4)
library(car)
library(MASS)
library(lme4)
library(knitr)
library(lmerTest)


#set directory
dir = "/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/2-spearmans/"


lag_duration = 28
spearman_results = data.frame()
everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv")

predictors = c("vcopy_per_100k", "nseg_relaxed")
outcomes = c("daily_count", "daily_rate", "tests_taken")



############################## Define dictionaries
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


outcome_name_dic = data.frame(outcome_var = c("daily_count", "daily_rate"),
                              outcome_name = c("Daily positive cases", "Daily positivity rate"))

predictor_name_dic = data.frame(predictor_var = c("vcopy_per_100k", "nseg_relaxed"),
                                predictor_name = c("Viral copy number", "# of segregating sites"))



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

calc_lagged_correlations = function(outcome){
  ## create dataframe to store this outcomes correlation results
  all_sites_correlation_df = data.frame()
  
  for(s in 1:nrow(sites_df)){ ##sites loop starts here
    
    ## get city data from everything
    city_id = sites_df$city_id[s]
    site_number = sites_df$site_number[s]
    
    city_data = everything[which(everything$city_id == city_id), ]
    city_data = city_data[which(city_data$site_number == site_number), ]
    city_data$date = as.Date(city_data$date)
    city_data = city_data[order(city_data$date), ]
    
    
    for(predictor in predictors){ ### var loop starts here
      cross_corr = as.data.frame(do.call(rbind, (lapply(-lag_duration:lag_duration, cross_correlation, 
                                                        city_data[, predictor], city_data[, outcome]))))
      cross_corr$city_id = city_id
      cross_corr$site_number = site_number
      cross_corr$comparison = paste0(predictor, " vs ", outcome)
      cross_corr$site_id = paste0(cross_corr$city_id, "_", cross_corr$site_number)
      
      all_sites_correlation_df = rbind(all_sites_correlation_df, cross_corr)
      
    } ### var loop ends here
    
  } ## sites loop ends here
  
  return(all_sites_correlation_df)
}

outcome = "daily_count"

df = calc_lagged_correlations(outcome)
df[df$pvalue >= 0.05, "estimate"] = 0
#df = df[df$estimate >= 0, ]
df = df[df$lag %in% -7:7, ]

model = lmer(estimate ~ comparison + (1|site_id),
             data = df)
summary(model)
