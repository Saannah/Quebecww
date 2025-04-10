#### script to generage new heatmap plots for Granger test results baased on 2025 version
rm(list = ls())

library(ggplot2)
library(stringr)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)


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
sites_df$site_identifier = paste0(sites_df$city_id, "_", sites_df$site_number)

## write a function that creates the summary dataframe for a given outcome
get_summary = function(outcome){
  max_order = 28
  
  predictors = c("nseg_relaxed", "vcopy_per_100k")
  all = data.frame()
  for(s in 1:nrow(sites_df)){
    all = rbind(all,
                read.delim(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/3-granger_permutations/", sites_df$site_identifier[s], "_outcome_", outcome, "_pred_vcopy_per_100k.tsv")),
                read.delim(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/3-granger_permutations/", sites_df$site_identifier[s], "_outcome_", outcome, "_pred_nseg_relaxed.tsv"))
    )
  }
  
  summary = data.frame(
    order = vector(),
    sig_count = vector(),
    total_count = vector(),
    sig_proportion = vector(),
    Count = vector(),
    causality = vector(),
    perm = vector())
  
  for(l in 1:max_order){
    
    for(pred in predictors){
      perms = all[which(all$order == l), ]
      perms = perms[which(perms$causality == paste0("y=",outcome, " x=", pred)), ]
      sigs = perms[perms$pvalue < 0.05, ]
      
      
      summary = rbind(summary,
                      data.frame(
                        order = l,
                        sig_count = nrow(sigs),
                        total_count = nrow(perms),
                        sig_proportion = nrow(sigs)/nrow(perms),
                        Count = (nrow(sigs)/nrow(perms))*10,
                        causality = paste0("y=",outcome, " x=", pred),
                        perm = "permutated data")
      )
    }
  }
  ## get real data and append
  real_data = read.delim(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/1-grangertest/counts_long_df_outcome_", outcome, ".tsv"))
  ## only keep forward causality rows
  real_data = real_data[which(str_detect(real_data$arrangement, paste0("y=", outcome))), ]
  
  summary = rbind(summary,
                  data.frame(
                    order = real_data$order,
                    sig_count = NA,
                    total_count = NA,
                    sig_proportion = NA,
                    Count = real_data$count,
                    causality = real_data$arrangement,
                    perm = real_data$arrangement)
                  
  )
  
  return(summary)
}




draw_complexheatmap = function(outcome){
  summary = get_summary(outcome)
  
  ### fix the inconsistency between vcopy and actual vcopy column names:
  #summary$causality = str_replace_all(summary$causality, "actual_vcopy", "vcopy")
  #summary$perm = str_replace_all(summary$perm, "actual_vcopy", "vcopy")
  
  summary$Prediction = summary$perm
  summary$Prediction[which(str_detect(summary$perm, "vcopy_per_100k"))] = "predictor viral copy number"
  summary$Prediction[which(str_detect(summary$perm, "nseg_relaxed"))] = "predictor number of segregating sites"
  
  summary$Prediction[(summary$perm == "permutated data" & 
                        (str_detect(summary$causality, "vcopy_per_100k")))] = "permutated data for viral copy number"
  
  summary$Prediction[(summary$perm == "permutated data" & 
                        (str_detect(summary$causality, "nseg_relaxed")))] = "permutated data for the number of segregating sites"
  
  
  # Create the matrix using xtabs
  matrix <- xtabs(Count ~ Prediction + order, data = summary)
  
  # Convert the xtabs object to a matrix
  matrix <- as.matrix(matrix)
  
  # assign column. names
  colnames(matrix) <- as.numeric(colnames(matrix))
  matrix <- matrix[c("predictor number of segregating sites",
                     "permutated data for the number of segregating sites",
                     "predictor viral copy number",
                     "permutated data for viral copy number"
  ), ]
  
  
  
  col_fun = colorRamp2(breaks = c(0, 10), colors = c("white", "#154EB7"))
  #col_fun = colorRamp2(breaks = c(0, 10), colors = c("white", "#AE1E1E"))
  #### try complex heatmaps
  cheatmap = Heatmap(matrix, name = str_wrap("Number of WW sites with significant Granger predictivity", 30),
                     width = ncol(matrix)*unit(9, "mm"), 
                     height = nrow(matrix)*unit(9, "mm"),
                     show_row_dend = FALSE,
                     show_column_dend = FALSE,
                     col = col_fun,
                     row_order = c("predictor number of segregating sites",
                                   "permutated data for the number of segregating sites",
                                   "predictor viral copy number",
                                   "permutated data for viral copy number"
                     ),
                     column_order = 1:28,
                     rect_gp = gpar(col = "black"),
                     row_labels = c("real data",
                                    "permuted data",
                                    "real data",
                                    "permuted data"),
                     column_labels = (1:28),
                     heatmap_legend_param = list(
                       title = str_wrap("Number of WW sites with significant Granger predictivity", 30),
                       #col_fun = color_palette,
                       at = (c(0,1,2,3,4,5,6,7,8,9,10)), 
                       legend_height = unit(5, "cm")),
                     row_split = c(str_wrap("a. # of segregating sites", width = 10), str_wrap("a. # of segregating sites", width = 10), str_wrap("b. Viral copy number", width = 10), str_wrap("b. Viral copy number",width = 10)),
                     row_gap = unit(5, "mm"),
                     row_title_rot = 0,
                     column_title = "Order",
                     column_names_rot = 0,
                     column_title_side = 'bottom',
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.rect(x, y, width, height,gp = gpar(col = "black", lwd = 2.7))
                       # if(matrix[i,j] >= 1 | matrix[i,j] == 0){
                       #   grid.text(sprintf(as.character(matrix[i, j])), x, y, gp = gpar(fontsize = 10))
                       # }
                       # if(matrix[i,j] > 0 & matrix[i,j] < 1){
                       #   grid.text(sprintf("< 1"), x, y, gp = gpar(fontsize = 10))
                       # }
                     },
                     column_names_centered = TRUE)
  
  outcome_name_dic = data.frame(outcome_var = c("daily_count", "daily_rate"),
                                outcome_name = c("Daily positive cases", "Daily positivity rate"))
  
  pdf(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/x-plots/", "granger_outcome_",outcome, "_new_heatmap.pdf"), width = 15, height = 5)
  draw(cheatmap, heatmap_legend_side = "right", 
       column_title = paste0(outcome_name_dic$outcome_name[outcome_name_dic$outcome_var == outcome], " ~ predictor"))
  dev.off()
  
}



draw_complexheatmap("daily_count")
draw_complexheatmap("daily_rate")


##################################################### Spearman's correlation heatmaps
### define get summary function for spearman's correlation'
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
                read.delim(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/4-spearmans_permutations//spearmans_permutations_", sites_df$site_identifier[s], "_outcome_", outcome, "_pred_vcopy_per_100k.tsv"))
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


draw_spearmans_complextheatmap = function(outcome, lag_range){
  summary = get_summary(outcome)
  
  summary = summary[summary$lag %in% lag_range, ]
  # ### fix the inconsistency between vcopy and actual vcopy column names:
  # summary$causality = str_replace_all(summary$causality, "actual_vcopy", "vcopy")
  # summary$perm = str_replace_all(summary$perm, "actual_vcopy", "vcopy")
  
  
  summary$Prediction = summary$perm
  summary$Prediction[which(str_detect(summary$perm, "vcopy_per_100k"))] = "viral copy number"
  summary$Prediction[which(str_detect(summary$perm, "nseg_relaxed"))] = "number of segregating sites"
  
  summary$Prediction[(summary$perm == "permutated data" & 
                        (str_detect(summary$test, "vcopy_per_100k")))] = "permutated data for viral copy number"
  
  summary$Prediction[(summary$perm == "permutated data" & 
                        (str_detect(summary$test, "nseg_relaxed")))] = "permutated data for the number of segregating sites"
  
  
  # Create the matrix using xtabs
  matrix <- xtabs(Count ~ Prediction + lag, data = summary)
  
  # Convert the xtabs object to a matrix
  matrix <- as.matrix(matrix)
  
  # assign column. names
  colnames(matrix) <- as.numeric(colnames(matrix))
  matrix <- matrix[c("number of segregating sites",
                     "permutated data for the number of segregating sites",
                     "viral copy number",
                     "permutated data for viral copy number"
  ), ]
  
  
  
  col_fun = colorRamp2(breaks = c(0, 10), colors = c("white", "#154EB7"))
  #col_fun = colorRamp2(breaks = c(0, 10), colors = c("white", "#AE1E1E"))
  #### try complex heatmaps
  cheatmap = Heatmap(matrix, name = str_wrap("Number of WW sites with positive, significant correlation", 30),
                     width = ncol(matrix)*unit(8, "mm"), 
                     height = nrow(matrix)*unit(8, "mm"),
                     show_row_dend = FALSE,
                     show_column_dend = FALSE,
                     col = col_fun,
                     row_order = c("number of segregating sites",
                                   "permutated data for the number of segregating sites",
                                   "viral copy number",
                                   "permutated data for viral copy number"
                     ),
                     column_order = 1:length(lag_range),
                     rect_gp = gpar(col = "black"),
                     row_labels = c("real data",
                                    "permuted data",
                                    "real data",
                                    "permuted data"),
                     column_labels = as.character(lag_range),
                     heatmap_legend_param = list(
                       title = str_wrap("Number of WW sites with positive, significant correlation", 30),
                       #col_fun = color_palette,
                       at = (c(0,1,2,3,4,5,6,7,8,9,10)), 
                       legend_height = unit(5, "cm")),
                     row_split = c(str_wrap("a. # of segregating sites", width = 10), str_wrap("a. # of segregating sites", width = 10), str_wrap("b. Viral copy number", width = 10), str_wrap("b. Viral copy number",width = 10)),
                     row_gap = unit(5.5, "mm"),
                     row_title_rot = 0,
                     column_title = "Time points lagged",
                     column_names_rot = 0,
                     column_title_side = 'bottom',
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.rect(x, y, width, height,gp = gpar(col = "black", lwd = 2.7))
                       # if(matrix[i,j] >= 1 | matrix[i,j] == 0){
                       #   grid.text(sprintf(as.character(matrix[i, j])), x, y, gp = gpar(fontsize = 10))
                       # }
                       # if(matrix[i,j] > 0 & matrix[i,j] < 1){
                       #   grid.text(sprintf("< 1"), x, y, gp = gpar(fontsize = 10))
                       # }
                     },
                     column_names_centered = TRUE)
  
  outcome_name_dic = data.frame(outcome_var = c("daily_count", "daily_rate"),
                                outcome_name = c("Daily positive cases", "Daily positivity rate"))
  
  pdf(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/x-plots/", "spearmans_outcome_",outcome, "_new_heatmap.pdf"), width = 16, height = 5)
  draw(cheatmap, heatmap_legend_side = "right", 
       column_title = paste0(outcome_name_dic$outcome_name[outcome_name_dic$outcome_var == outcome], " vs. metric"))
  dev.off()
}

draw_spearmans_complextheatmap("daily_count", -7:28)
draw_spearmans_complextheatmap("daily_rate", -7:28)

