## make plot of model orders
rm(list =ls())


dir = "/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/8-irf_w_varselect/1-varselect_grids/"

ordertab = read.delim(paste0("/users/sana/Documents/CANMOD/16-updated_results_for_manuscript/MS-09-IRF_with_95_CI_no_smoothing_varselect/varselect_grid_outcome_", outcome, "_RD68.tsv"))

nseg_matrix_form = matrix("0",
                          nrow = 10,
                          ncol = 15)

vcopy_matrix_form = matrix("0",
                           nrow = 10,
                           ncol = 15)

colnames(nseg_matrix_form) = (1:15)
rownames(nseg_matrix_form) = unique(ordertab$site)
colnames(vcopy_matrix_form) = (1:15)
rownames(vcopy_matrix_form) = unique(ordertab$site)



for(i in 1:nrow(ordertab)){
  if(ordertab$predictor[i] == "# of SNPs"){
    nseg_matrix_form[ordertab$site[i], ordertab$best_model_order[i]] = "1"
  }
  
  if(ordertab$predictor[i] == "Viral copy number"){
    vcopy_matrix_form[ordertab$site[i], ordertab$best_model_order[i]] = "1"
  }
  
}



#col_fun = colorRamp2(breaks = c(0,1), colors = c("white", "darkred"))
#col_fun = colorRamp2(breaks = c(0, 10), colors = c("white", "#AE1E1E"))
#### try complex heatmaps

matrix = rbind(nseg_matrix_form, vcopy_matrix_form)
colors = structure(c("white", "#618B25"), names = c("0", "1"))
cheatmap = Heatmap(matrix, name = str_wrap("IRF", 30),
                   width = ncol(matrix)*unit(22, "mm"), 
                   height = nrow(matrix)*unit(22, "mm"),
                   show_row_dend = FALSE,
                   show_column_dend = FALSE,
                   col = colors,
                   # row_order = c(
                   #   "MTL_1",
                   #   "MTL_2",
                   #   "LVL_1",
                   #   "LVL_2",
                   #   "LVL_5",
                   #   "GTN_1",
                   #   "Riki_1",
                   #   "Sept_1",
                   #   "QC_1",
                   #   "QC_2",
                   #   "MTL_1",
                   #   "MTL_2",
                   #   "LVL_1",
                   #   "LVL_2",
                   #   "LVL_5",
                   #   "GTN_1",
                   #   "Riki_1",
                   #   "Sept_1",
                   #   "QC_1",
                   #   "QC_2"),
                   column_order = 1:ncol(matrix),
                   rect_gp = gpar(col = "black"),
                   # row_labels = c("real data",
                   #                "permuted data",
                   #                "real data",
                   #                "permuted data"),
                   #column_labels = (1:28),
                   heatmap_legend_param = list(
                     title = str_wrap("Best model order (yes/no):", 30),
                     title_gp = gpar(fontsize = 30),
                     labels_gp = gpar(fontsize = 20),
                     labels = c("No" , "Yes"),
                     at = c("0", "1"),
                     #legend_gp = gpar(fill=col_fun(at)),
                     #col_fun = color_palette,
                     #at = (c(0,1,2,3,4,5,6,7,8,9,10)), 
                     legend_height = unit(10, "cm"),
                     grid_width = unit(2, "cm")),
                   row_split = c(rep("# of segregating sites", 10), rep("Viral copy number", 10)),
                   row_gap = unit(10, "mm"),
                   row_title_rot = 90,
                   row_title_gp = gpar(fontsize = 35),
                   row_names_gp = gpar(fontsize = 25),
                   column_title = "Model order",
                   column_title_gp = gpar(fontsize = 40),
                   column_names_rot = 0,
                   column_names_gp = gpar(fontsize = 25),
                   column_title_side = 'bottom',
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.rect(x, y, width, height,gp = gpar(col = "black", lwd = 4))
                     # if(matrix[i,j] != 0){
                     #   grid.text(sprintf("%.2f", matrix[i,j]), x, y, gp = gpar(fontsize = 20))
                     # }
                     # if(matrix[i,j] > 0 & matrix[i,j] < 1){
                     #   grid.text(sprintf("< 1"), x, y, gp = gpar(fontsize = 10))
                     # }
                   },
                   column_names_centered = TRUE)

#system(paste0("mkdir ", dir, "heatmaps"))
pdf(paste0(dir, "best_order_grids/Grid_best_order_", outcome, ".pdf"), 
    width = 25, height = 25)
draw(cheatmap, heatmap_legend_side = "right")
grid.text(str_wrap(paste0("Vector Autoregressive (VAR) models fitted for positive cases and each metric"), width = 90),
          x = unit(0.5, "npc"), y = unit(0.9, "npc"), # Position the title
          gp = gpar(fontsize = 35))
dev.off()













############################################################################################################
## make map of model orders
rm(list =ls())
#### daily count

dir = "/users/sana/Documents/CANMOD/16-updated_results_for_manuscript/MS-09-IRF_with_95_CI_no_smoothing_varselect/"

outcome = "daily_rate"
ordertab = read.delim(paste0("/users/sana/Documents/CANMOD/16-updated_results_for_manuscript/MS-09-IRF_with_95_CI_no_smoothing_varselect/varselect_grid_outcome_", outcome, "_RD68.tsv"))

nseg_matrix_form = matrix("0",
                          nrow = 10,
                          ncol = 15)

vcopy_matrix_form = matrix("0",
                           nrow = 10,
                           ncol = 15)

colnames(nseg_matrix_form) = (1:15)
rownames(nseg_matrix_form) = unique(ordertab$site)
colnames(vcopy_matrix_form) = (1:15)
rownames(vcopy_matrix_form) = unique(ordertab$site)



for(i in 1:nrow(ordertab)){
  if(ordertab$predictor[i] == "# of SNPs"){
    nseg_matrix_form[ordertab$site[i], ordertab$best_model_order[i]] = "1"
  }
  
  if(ordertab$predictor[i] == "Viral copy number"){
    vcopy_matrix_form[ordertab$site[i], ordertab$best_model_order[i]] = "1"
  }
  
}



#col_fun = colorRamp2(breaks = c(0,1), colors = c("white", "darkred"))
#col_fun = colorRamp2(breaks = c(0, 10), colors = c("white", "#AE1E1E"))
#### try complex heatmaps

matrix = rbind(nseg_matrix_form, vcopy_matrix_form)
colors = structure(c("white", "#618B25"), names = c("0", "1"))
cheatmap = Heatmap(matrix, name = str_wrap("IRF", 30),
                   width = ncol(matrix)*unit(22, "mm"), 
                   height = nrow(matrix)*unit(22, "mm"),
                   show_row_dend = FALSE,
                   show_column_dend = FALSE,
                   col = colors,
                   # row_order = c(
                   #   "MTL_1",
                   #   "MTL_2",
                   #   "LVL_1",
                   #   "LVL_2",
                   #   "LVL_5",
                   #   "GTN_1",
                   #   "Riki_1",
                   #   "Sept_1",
                   #   "QC_1",
                   #   "QC_2",
                   #   "MTL_1",
                   #   "MTL_2",
                   #   "LVL_1",
                   #   "LVL_2",
                   #   "LVL_5",
                   #   "GTN_1",
                   #   "Riki_1",
                   #   "Sept_1",
                   #   "QC_1",
                   #   "QC_2"),
                   column_order = 1:ncol(matrix),
                   rect_gp = gpar(col = "black"),
                   # row_labels = c("real data",
                   #                "permuted data",
                   #                "real data",
                   #                "permuted data"),
                   #column_labels = (1:28),
                   heatmap_legend_param = list(
                     title = str_wrap("Best model order (yes/no):", 30),
                     title_gp = gpar(fontsize = 30),
                     labels_gp = gpar(fontsize = 20),
                     labels = c("No" , "Yes"),
                     at = c("0", "1"),
                     #legend_gp = gpar(fill=col_fun(at)),
                     #col_fun = color_palette,
                     #at = (c(0,1,2,3,4,5,6,7,8,9,10)), 
                     legend_height = unit(10, "cm"),
                     grid_width = unit(2, "cm")),
                   row_split = c(rep("# of segregating sites", 10), rep("Viral copy number", 10)),
                   row_gap = unit(10, "mm"),
                   row_title_rot = 90,
                   row_title_gp = gpar(fontsize = 35),
                   row_names_gp = gpar(fontsize = 25),
                   column_title = "Model order",
                   column_title_gp = gpar(fontsize = 40),
                   column_names_rot = 0,
                   column_names_gp = gpar(fontsize = 25),
                   column_title_side = 'bottom',
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.rect(x, y, width, height,gp = gpar(col = "black", lwd = 4))
                     # if(matrix[i,j] != 0){
                     #   grid.text(sprintf("%.2f", matrix[i,j]), x, y, gp = gpar(fontsize = 20))
                     # }
                     # if(matrix[i,j] > 0 & matrix[i,j] < 1){
                     #   grid.text(sprintf("< 1"), x, y, gp = gpar(fontsize = 10))
                     # }
                   },
                   column_names_centered = TRUE)

#system(paste0("mkdir ", dir, "heatmaps"))
pdf(paste0(dir, "best_order_grids/Grid_best_order_", outcome, ".pdf"), 
    width = 25, height = 25)
draw(cheatmap, heatmap_legend_side = "right")
grid.text(str_wrap(paste0("Vector Autoregressive (VAR) models fitted for daily positivity rate and each metric"), width = 90),
          x = unit(0.5, "npc"), y = unit(0.9, "npc"), # Position the title
          gp = gpar(fontsize = 35))
dev.off()


