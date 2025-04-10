library(stringr)

tab = read.csv("/users/sana/Downloads/Analysis_Report-Slim-HQ-2025-03-06.csv")
tab = read.csv("/users/sana/Downloads/Analysis_Report-2025-03-06.csv")
# tab = tab[tab$Breadth.of.coverage.Genome > 0.1, ]
# tab = tab[tab$Mean.depth.Genome >= 20, ]

everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files/cases_estimates_combined_6cities_archdownload_qpcr_wFreyja_slimHQ_posrate_weektotal_hosp_newnseg.tsv")


dim(tab)
dim(everything)

length(which(everything$sample %in% tab$Sequence))

vcfnames = read.delim("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/vcf_names", header = FALSE)

# vcfnames$sample = str_remove_all(vcfnames$V1, "_variantout.vcf")
# vcfnames$sample = str_remove_all(vcfnames$sample, "_variantout.snpEff.annotated.vcf")

names = vcfnames[which(str_detect(vcfnames$V1, "_variantout.vcf")), 1]
names = str_remove_all(names, "_variantout.vcf")
length(which(everything$sample %in% names))

length(which(tab$Sequence %in% names))
dim(tab)


everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files/cases_estimates_combined_6cities_archdownload_qpcr_wFreyja_slimHQ_posrate_weektotal_hosp_newnseg.tsv")
tab = read.csv("/users/sana/Downloads/Analysis_Report-Slim-2025-03-06.csv")
noseg = read.delim("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/nosegsites", header = FALSE)
noseg = noseg[, 1]
noseg = str_remove(noseg, "_variantout.parsed.NoSegSitesCalled")

length(which(noseg %in% tab$Sequence))



########## load nseg redo
rm(list = ls())
redonseg = read.delim("/users/sana/Documents/CANMOD/17-updated_vcf_files_2025/REDO_01_output_nseg_estimates_corrected.tsv")
everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files/cases_estimates_combined_6cities_archdownload_qpcr_wFreyja_slimHQ_posrate_weektotal_hosp_newnseg.tsv")
everything = everything[, c("Sequence", "nseg_relaxed")]


redonseg$Sequence = str_remove(redonseg$filename, "_variantout_tabular.vcf")
redonseg = redonseg[which(redonseg$Sequence %in% everything$Sequence), ]
redonseg = redonseg[, c("Sequence", "nseg_relaxed")]

all = merge(redonseg, everything, "Sequence")
length(which(!all$nseg_relaxed.x == all$nseg_relaxed.y))

sequences_different = all[which(!all$nseg_relaxed.x == all$nseg_relaxed.y),]
all_diff = all[, ]

tab = read.csv("/users/sana/Downloads/Analysis_Report-2025-03-06.csv")
sub_tab = tab[which(tab$Sequence %in% sequences_different), ]

