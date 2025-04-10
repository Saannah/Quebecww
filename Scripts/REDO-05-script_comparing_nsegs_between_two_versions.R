rm(list = ls())

## compare previously calculated nsegs with the new ones

new = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs.tsv")
previous = read.delim("/users/sana/Documents/CANMOD/0-reference_files/cases_estimates_combined_6cities_archdownload_qpcr_wFreyja_slimHQ_posrate_weektotal_hosp_newnseg.tsv")

common_seqs = new$Sequence[which(new$Sequence %in% previous$Sequence)]



### check if all previous sequences are in new
length(which(previous$Sequence %in% new$Sequence))
length(which(new$Sequence %in% previous$Sequence))

merged = merge(new[which(new$Sequence %in% common_seqs), c("Sequence", "nseg_relaxed")],
               previous[which(previous$Sequence %in% common_seqs), c("Sequence", "nseg_relaxed")], "Sequence")

length(which(merged$nseg_relaxed.x == merged$nseg_relaxed.y))
