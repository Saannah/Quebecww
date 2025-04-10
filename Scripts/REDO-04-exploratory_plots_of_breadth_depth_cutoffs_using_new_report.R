rm(list = ls())
data = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs.tsv")

breadth_table_long = data.frame()

for(breadth in seq(0.1, 0.9, 0.1)){
  sub_tab = data[which(data$Breadth.of.coverage.Genome >= breadth), ]
  sub_tab$breadth_cutoff = as.character(breadth)
  breadth_table_long = rbind(breadth_table_long,
                             sub_tab)
}


p = ggplot(data = breadth_table_long) +
  geom_boxplot(aes(x = breadth_cutoff, y = nseg_relaxed), alpha = 0.5) +
  theme_bw() +
  theme(text = element_text(size = 30))+
  ylim(c(0, 650))

ggsave(p,
       file = "/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/breadth.pdf",
       device = "pdf",
       limitsize = FALSE)


depth_table_long = data.frame()

for(depth in seq(1, 100, 10)){
  sub_tab = data[which(data$Mean.depth.Genome >= depth), ]
  sub_tab$depth_cutoff = (depth)
  depth_table_long = rbind(depth_table_long,
                           sub_tab)
}



p = ggplot(data = depth_table_long) +
  geom_boxplot(aes(x = depth_cutoff, y = nseg_relaxed, group = depth_cutoff)) +
  theme_bw() +
  theme(text = element_text(size = 20))+
  ylim(c(0, 2000))

ggsave(p,
       file = "/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/depth.pdf",
       device = "pdf",
       limitsize = FALSE,
       height = 20,
       width = 20)



############################### depth and breadth vs nseg

p = ggplot(data = data)+
  geom_point(aes(x = Breadth.of.coverage.Genome, y = nseg_relaxed), alpha = 0.5, size = 10)+
  theme_bw()+
  theme(text = element_text(size = 100))+
  ylab("Number of segregating sites called")+
  xlab("Breadth of coverage in the sample")

ggsave(p,
       file = "/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/corr_breadth.pdf",
       device = "pdf",
       limitsize = FALSE,
       height = 100,
       width = 100)


p = ggplot(data = data)+
  geom_point(aes(x = Mean.depth.Genome, y = nseg_relaxed), alpha = 0.5, size = 10)+
  theme_bw()+
  theme(text = element_text(size = 100))+
  ylab("Number of segregating sites called")+
  xlab("Mean depth")

ggsave(p,
       file = "/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/corr_depth.pdf",
       device = "pdf",
       limitsize = FALSE,
       height = 100,
       width = 100)

################ NA's in report
everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files/cases_estimates_combined_6cities_archdownload_qpcr_wFreyja_slimHQ_posrate_weektotal_hosp_newnseg.tsv")

na_ind = report[which(is.na(report$Breadth.of.coverage.Genome)), ]
