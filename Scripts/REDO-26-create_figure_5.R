####### create figure 5 with CI 95% taken into account + varselect
### plot IRF with 95% CI and 10k iterations calculated by MS-22 and MS-23
rm(list = ls())
library(vars)
library(tseries)
library(lmtest)
library(ggplot2)
library(pracma)
library(stringi)
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




### daily count
outcome = "daily_count"
s=1

out_nseg = read.delim("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/8-irf_w_varselect/3-irf_vals/cumulative_all_sites_2varmodel_outcome_daily_count_pred_nseg_relaxed_model_varselect_seed_1_iterations_10000.tsv")
out_vcopy = read.delim("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/8-irf_w_varselect/3-irf_vals/cumulative_all_sites_2varmodel_outcome_daily_count_pred_vcopy_per_100k_model_varselect_seed_1_iterations_10000.tsv")
out = rbind(out_nseg, out_vcopy)


out$sig = NA
out$sig[out$lower_irf > 0] = out$irf_vals[out$lower_irf > 0]
out$sig[out$upper_irf < 0] = out$irf_vals[out$upper_irf < 0]

out = out[out$city_id == "MTL" & out$site_number == 1, ]

p_nseg = ggplot(data = out[out$shock_applied_to == "nseg_relaxed", ]) + 
  geom_line(aes(x=time_point, y=irf_vals), color = "#3943B7") +
  geom_ribbon(aes(ymin=lower_irf, ymax=upper_irf, x=time_point), alpha = 0.3,  fill = "#3943B7")+
  #geom_line(aes(x=time_point, y=upper_irf), color = "blue", linetype = 2)+
  #geom_line(aes(x=time_point, y=lower_irf), color = "blue", linetype = 2)+
  geom_point(aes(x=time_point, y=sig), shape = 8, color = "#3943B7", size = 3, stroke = 1)+
  theme_bw()+
  ylim(c(-0.35,0.35))+
  ylab(str_wrap(paste0("Cumulative response in positive case counts"), width = 20)) +
  xlab("Time points after shock")+
  ggtitle("a. Montreal 1 - VAR model of order 5", "Response to a shock in the number of segregating sites")+
  theme(text = element_text(size = 20))

p_vcopy = ggplot(data = out[out$shock_applied_to == "vcopy_per_100k", ]) + 
  geom_line(aes(x=time_point, y=irf_vals), color = "#A51C30") +
  geom_ribbon(aes(ymin=lower_irf, ymax=upper_irf, x=time_point), alpha = 0.3,  fill = "#A51C30")+
  #geom_line(aes(x=time_point, y=upper_irf), color = "#A51C30", linetype = 2)+
  #geom_line(aes(x=time_point, y=lower_irf), color = "#A51C30", linetype = 2)+
  geom_point(aes(x=time_point, y=sig), shape = 8, color = "#A51C30", size = 1, stroke = 1)+
  theme_bw()+
  ylim(c(-0.35,0.35))+
  ylab(str_wrap(paste0("Cumulative response in positive case counts"), width = 20))+
  xlab("Time points after shock")+
  ggtitle("b. Montreal 1 - VAR model of order 9", "Response to a shock in viral copy number")+
  theme(text = element_text(size = 20))


both = ggarrange(p_nseg, p_vcopy, nrow = 2, ncol = 1)
response_type = out$response_type[1]


ggsave(both, 
       file = paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/x-plots/fig5_", response_type, sites_df$city_id[s], "_", sites_df$site_number[s], "_2varmodel_", outcome, "_model_varselect", "_10k_iterations.pdf"),
       device = "pdf",
       height = 10,
       width = 20)
ggsave(both, 
       file = paste0("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/x-plots/fig5_", response_type, sites_df$city_id[s], "_", sites_df$site_number[s], "_2varmodel_", outcome, "_model_varselect", "_10k_iterations.png"),
       device = "png",
       height = 10,
       width = 20)

