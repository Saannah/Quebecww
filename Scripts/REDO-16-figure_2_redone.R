#### plot time series for manuscript.
rm(list = ls())

legend_x = 0.5
legend_y = 0.8
#cases_color = "#A0AF84"
cases_color = "#D4AFB9"
cases_color = "#D8C0D7"
dir = "/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/x-plots/"

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

everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv")

predictors = c("nseg_relaxed", "vcopy_per_100k")
outcomes = c("daily_count", "daily_rate","tests_taken")

s=1
city_id = sites_df$city_id[s]
site_number = sites_df$site_number[s]

city_data = everything[which(everything$city_id == city_id), ]
city_data = city_data[which(city_data$site_number == site_number), ]
city_data$date = as.Date(city_data$date)
city_data = city_data[order(city_data$date), ]


fake_df = data.frame(metric = c("Number of segregating sites", "Positive cases", "Number of segregating sites", "Positive cases"),
                     value = NA)
# p = ggplot(data = fake_df) + geom_line(aes(x=metric, y=value, color = metric)) +
#   scale_color_manual(values = c("#274690", "#E69F00"))
# fake_leg = (get_legend(p))
# #ggsave(fake_leg, file = "nseg_legend.pdf", device = "pdf", height = 2, width = 4)

city_data$nseg_col_code = "The number of segregating sites"
city_data$nseg_col_code[1] = "Positive cases"

coeff1 = (max(city_data$nseg_relaxed, na.rm = TRUE)/max(city_data$daily_count, na.rm = TRUE))
p_nseg = ggplot(data = city_data) + 
  #geom_area(aes(x = date, y=(nseg_relaxed) ), fill = "#E69F00", alpha = 0.6) + 
  #geom_point(aes(x = date, y = daily_count * coeff1), alpha = 0.7, color = "#E69F00", size = 3) +
  #geom_line(aes(x = date, y = daily_count * coeff1), color = cases_color, linewidth = 2, alpha = 0.7) + 
  geom_area(aes(x = date, y = daily_count * coeff1), fill = cases_color, alpha = 0.6)+
  geom_line(aes(x = date, y=(nseg_relaxed) ), color = "#3943B7", linewidth = 2.5) + 
  geom_point(aes(x = date, y=(nseg_relaxed) ), color = "#3943B7", size = 2) + 
  geom_line(aes(x = date, y=12, color = nseg_col_code ), linewidth = 0, alpha =0)+
  
  scale_color_manual(values = c("#3943B7", cases_color), labels = c(str_wrap("The number of segregating sites", width = 15), str_wrap("Positive cases", width = 15))) + 
  labs(color = "Metric:")+
  guides(linetype = guide_legend(override.aes = list(size = 20)))+
  # scale_fill_manual(values = c("#E69F00", "#3943B7"))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.key.size = unit(3.5, "lines"),
        legend.position = c(legend_x,legend_y),
        text = element_text(size = 50)) + 
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 3, alpha = 0.8)  # Set line thickness in the legend keys
    )
  )+
  scale_x_date(date_breaks = "1 week") + 
  scale_y_continuous(name = "# of segregating sites", sec.axis = sec_axis(~.*(1/coeff1), name="Positive cases")) + 
  ggtitle("a. The number of segregating sites") +
  xlab("Date")


# ggsave(p_nseg, 
#        file = paste0(dir, "nseg_plot.pdf"), 
#        device = "pdf", 
#        height = 16, width = 48, 
#        limitsize = FALSE)

city_data$vcopy_col_code = "Viral copy number"
city_data$vcopy_col_code[1] = "Positive cases"

coeff2 = (max(city_data$vcopy_per_100k, na.rm = TRUE)/max(city_data$daily_count, na.rm = TRUE))
p_vcopy = ggplot(city_data) + 
  ### fake layer
  geom_line(aes(x = date, y=12, color = vcopy_col_code ), linewidth = 0, alpha = 0)+
  
  #geom_area(aes(x = date, y=(nseg_relaxed) ), fill = "#E69F00", alpha = 0.6) + 
  #geom_point(aes(x = date, y = daily_count * coeff2), alpha = 0.7, color = "#E69F00", size = 3) +
  #geom_line(aes(x = date, y = daily_count * coeff2), color = cases_color, linewidth = 2, alpha = 0.7) + 
  geom_area(aes(x = date, y = daily_count * coeff2), fill = cases_color, alpha = 0.6)+
  geom_line(aes(x = date, y=(vcopy_per_100k) ), color = "#A51C30", linewidth = 2.5) + 
  geom_point(aes(x = date, y=(vcopy_per_100k) ), color = "#A51C30", size = 2) + 
  
  
  scale_color_manual(values = c("#A51C30", cases_color), labels = c(str_wrap("Viral copy number", width = 15), str_wrap("Positive cases", width = 14))) + 
  labs(color = "Metric:")+
  guides(linetype = guide_legend(override.aes = list(size = 20)))+
  # scale_fill_manual(values = c("#E69F00", "#3943B7"))+
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.key.size = unit(5.5, "lines"),
        legend.position = c(legend_x,legend_y),
        text = element_text(size = 50)) + 
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 3, alpha = 0.8)  # Set line thickness in the legend keys
    )
  )+
  scale_x_date(date_breaks = "1 week") + 
  scale_y_continuous(name = "Viral copy number", sec.axis = sec_axis(~.*(1/coeff2), name="Positive cases")) + 
  ggtitle("b. Viral copy number") +
  xlab("Date")

# ggsave(p_vcopy, 
#        file = paste0(dir, "vcopy_plot.pdf"), 
#        device = "pdf", 
#        height = 16, width = 48, 
#        limitsize = FALSE)


both = ggarrange(p_nseg, p_vcopy, nrow=2, ncol=1)
ggsave(both, 
       file = paste0(dir, "FIG2-time_series_mtl1.pdf"), 
       device = "pdf", 
       height = 40, width = 48, 
       limitsize = FALSE)

ggsave(both, 
       file = paste0(dir, "FIG2-time_series_mtl1.jpeg"), 
       device = "jpeg", 
       height = 40, width = 48, 
       limitsize = FALSE)
