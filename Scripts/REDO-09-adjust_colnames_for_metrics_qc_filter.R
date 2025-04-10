rm(list = ls())
#### add columns with consistent names in the final data set
data = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated.tsv")

data$daily_count = data$Tests.positifs
data$daily_rate = (data$Tests.positifs/(data$Tests.negatifs + data$Tests.positifs)) * 100
data$tests_taken = (data$Tests.negatifs + data$Tests.positifs)
data$weekly_count = data$week_total_positive_count


colnames(data)[which(colnames(data) == "City_ID")] <- "city_id"
colnames(data)[which(colnames(data) == "Site_Number")] <- "site_number"

########## filter based on criteria
report = data
report = report[which(report$Number.of.Reads >= 100000), ]
report = report[which(report$Percentage.of.Human.Reads <= 20),]
report = report[which(report$Mean.depth.Genome >= 20), ]
report = report[which(report$Breadth.of.coverage.Genome >= 0.1), ]
report = report[which(report$Percentage.of.SARS.CoV.2.reads >= 1), ]


write.table(report,
            file = "/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated_colnames_edited_qc_filtered.tsv",
            sep = "\t",
            row.names = FALSE)
