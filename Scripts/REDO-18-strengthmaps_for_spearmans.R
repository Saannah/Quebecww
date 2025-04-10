### create grid of correlation values
rm(list = ls())
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


make_strengthmaps = function(predictor, outcome){
  df = calc_lagged_correlations(outcome)
  pred_df = df[str_detect(df$comparison, predictor), ]
  pred_df$strength = 0
  pred_df$strength[pred_df$pvalue < 0.05] = pred_df$estimate[pred_df$pvalue < 0.05]
  
  # Create the matrix using xtabs
  matrix <- xtabs(strength ~ site_id + lag, data = pred_df)
  
  # Convert the xtabs object to a matrix
  matrix <- as.matrix(matrix)
  
  # assign column. names
  colnames(matrix) <- as.numeric(colnames(matrix))
  # matrix <- matrix[c("predictor number of segregating sites",
  #                    "permutated data for the number of segregating sites",
  #                    "predictor viral copy number",
  #                    "permutated data for viral copy number"
  # ), ]
  col_fun = colorRamp2(breaks = c(-1, 0,1), colors = c("red", "white","#154EB7"))
  #col_fun = colorRamp2(breaks = c(0, 10), colors = c("white", "#AE1E1E"))
  #### try complex heatmaps
  cheatmap = Heatmap(matrix, name = str_wrap("Strength of correlation", 30),
                     width = ncol(matrix)*unit(22, "mm"), 
                     height = nrow(matrix)*unit(22, "mm"),
                     show_row_dend = FALSE,
                     show_column_dend = FALSE,
                     col = col_fun,
                     row_order = c(
                       "MTL_1",
                       "MTL_2",
                       "LVL_1",
                       "LVL_2",
                       "LVL_5",
                       "GTN_1",
                       "Riki_1",
                       "Sept_1",
                       "QC_1",
                       "QC_2"),
                     column_order = 1:ncol(matrix),
                     rect_gp = gpar(col = "black"),
                     # row_labels = c("real data",
                     #                "permuted data",
                     #                "real data",
                     #                "permuted data"),
                     #column_labels = (1:28),
                     heatmap_legend_param = list(
                       title = str_wrap("Strength of correlation", 30),
                       title_gp = gpar(fontsize = 30),
                       labels_gp = gpar(fontsize = 20),
                       #col_fun = color_palette,
                       #at = (c(0,1,2,3,4,5,6,7,8,9,10)), 
                       legend_height = unit(10, "cm"),
                       grid_width = unit(2, "cm")),
                     #row_split = c(str_wrap("# of segregating sites", width = 10), str_wrap("# of segregating sites", width = 10), str_wrap("Viral copy number", width = 10), str_wrap("Viral copy number",width = 10)),
                     #row_gap = unit(5, "mm"),
                     row_title_rot = 0,
                     row_names_gp = gpar(fontsize = 25),
                     column_title = "Order",
                     column_title_gp = gpar(fontsize = 40),
                     column_names_rot = 0,
                     column_names_gp = gpar(fontsize = 25),
                     column_title_side = 'bottom',
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.rect(x, y, width, height,gp = gpar(col = "black", lwd = 4))
                       if(matrix[i,j] != 0){
                         grid.text(sprintf("%.2f", matrix[i,j]), x, y, gp = gpar(fontsize = 20))
                       }
                       # if(matrix[i,j] > 0 & matrix[i,j] < 1){
                       #   grid.text(sprintf("< 1"), x, y, gp = gpar(fontsize = 10))
                       # }
                     },
                     column_names_centered = TRUE)
  
  
  pdf(paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/x-plots/strengthmaps_outcome_",
             outcome,
             "_predictor_",
             predictor,
             ".pdf"), 
      width = 58, height = 15)
  draw(cheatmap, heatmap_legend_side = "right")
  grid.text(paste0(predictor_name_dic$predictor_name[predictor_name_dic$predictor_var == predictor] ,
                   " vs ",
                   outcome_name_dic$outcome_name[outcome_name_dic$outcome_var == outcome]),
            x = unit(0.5, "npc"), y = unit(0.85, "npc"), # Position the title
            gp = gpar(fontsize = 40, fontface = "bold")) # Adjust the font size and style
  dev.off()
}

make_strengthmaps("nseg_relaxed", "daily_count")
make_strengthmaps("nseg_relaxed", "daily_rate")
make_strengthmaps("vcopy_per_100k", "daily_count")
make_strengthmaps("vcopy_per_100k", "daily_rate")

