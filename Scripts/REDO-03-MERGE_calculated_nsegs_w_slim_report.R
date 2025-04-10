rm(list = ls())

library(stringr)
library(ggplot2)

### load new analysis report
report = read.delim(file = "/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06.csv", sep = ",")


redonseg = read.delim("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/REDO_01_output_nseg_estimates_corrected.tsv")
redonseg$Sequence = str_remove(redonseg$filename, "_variantout_tabular.vcf")
redonseg = redonseg[, c("Sequence", "nseg_relaxed")]

length(which(report$Sequence %in% redonseg$Sequence))

redonseg$Sequence[which(redonseg$Sequence %in% report$Sequence)]


## loadd sequences that had nothing called
noseg = read.delim("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/nosegsites", header = FALSE)
noseg = noseg[, 1]
noseg = str_remove(noseg, "_variantout.parsed.NoSegSitesCalled")
length(noseg %in% report$Sequence)




redid_nseg_df = rbind(redonseg,
                      data.frame(Sequence = noseg,
                                 nseg_relaxed = 0))
length(which(report$Sequence %in% redid_nseg_df$Sequence))
length(which(redid_nseg_df$Sequence %in% report$Sequence))

redid_nseg_df = redid_nseg_df[which(redid_nseg_df$Sequence %in% report$Sequence), ]



missing_report = report[which(!report$Sequence %in% redid_nseg_df$Sequence), ]

data = merge(redid_nseg_df, report, "Sequence")

write.table(data,
            file = "/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs.tsv",
            sep = "\t",
            row.names = FALSE)

