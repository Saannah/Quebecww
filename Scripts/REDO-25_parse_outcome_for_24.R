rm(list = ls())
library(stringr)
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

setwd("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/8-irf_w_varselect/3-irf_vals/")
filenames = read.delim("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/8-irf_w_varselect/3-irf_vals/tsvnames", header = FALSE)
filenames = filenames[, 1]


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

grid$identifier = paste0("_2varmodel_outcome_",grid$outcome, "_pred_", grid$predictor,"_model_varselect_seed_", grid$seed, "_iterations_", grid$its)


 
 
create_heatmap_of_significant_irf = function(out, shocked){
  
  out$significant_irf = out$sig
  out$significant_irf[is.na(out$sig)] = 0
  out$site_id = paste0(out$city_id, "_",out$site_number)
  
  
  # Create the matrix using xtabs
  matrix <- xtabs(significant_irf ~ site_id + time_point, data = out[out$shock_applied_to == shocked, ])
  
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
  cheatmap = Heatmap(matrix, name = str_wrap("IRF", 30),
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
                       title = str_wrap("Response", 30),
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
  
  #system(paste0("mkdir ", dir, "heatmaps"))
  pdf(paste0("heatmaps/hmaps_",
             identifier,
             "_",
             shocked,
             ".pdf"), 
      width = 58, height = 15)
  draw(cheatmap, heatmap_legend_side = "right")
  grid.text(paste0("shock applied to ", shocked, " - outcome: ", unique(out$response)),
            x = unit(0.5, "npc"), y = unit(0.85, "npc"), # Position the title
            gp = gpar(fontsize = 40))
  dev.off()
  
}

# system("mkdir heatmaps")

for(identifier in unique(grid$identifier)){
  out = read.delim(paste0("cumulative_all_sites", identifier, ".tsv"))
  create_heatmap_of_significant_irf(out, unique(out$shock_applied_to))
}

