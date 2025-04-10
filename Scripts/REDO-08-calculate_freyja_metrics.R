rm(list = ls())
library(stringr)

data = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr.tsv")

calculate_shannon_for_1_sample = function(rownum, all_tab){
  sample = all_tab[rownum, ]
  abundance_vector = as.numeric(unlist(str_split(sample$abundances, pattern = " ")))
  
  unfiltered_shannon = -sum(log(abundance_vector) * abundance_vector)
  unfiltered_richness = length(abundance_vector)
  
  if(length(which(abundance_vector < 0.05)) > 0){
    abundance_vector = abundance_vector[-which(abundance_vector < 0.05)]
  }
  
  if(length(abundance_vector) == 0){
    filtered_shannon = "lineage split"
    filtered_richness = "lineage split"
    return(list(unfiltered_shannon = unfiltered_shannon,
                filtered_shannon = filtered_shannon,
                unfiltered_richness = unfiltered_richness,
                filtered_richness = filtered_richness))
    
  }
  
  filtered_shannon = -sum(log(abundance_vector) * abundance_vector)
  filtered_richness = length(abundance_vector)
  
  return(data.frame(unfiltered_shannon = unfiltered_shannon,
                    filtered_shannon = filtered_shannon,
                    unfiltered_richness = unfiltered_richness,
                    filtered_richness = filtered_richness))
  
}


freyja_metrics = do.call(rbind, (lapply(1:nrow(data), calculate_shannon_for_1_sample, data)))
data_w_freyja = cbind(data, freyja_metrics)

##### sanity check
everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files/cases_estimates_combined_6cities_archdownload_qpcr_wFreyja_slimHQ_posrate_weektotal_hosp_newnseg.tsv")


old_fil_richness = everything$filtered_richness
new_fil_richness = data_w_freyja$filtered_richness[match(everything$Sequence, data_w_freyja$Sequence)]
length(which(old_fil_richness != new_fil_richness))

#### there are two sequences that have new freyja calls
write.table(data_w_freyja,
            file = "/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr_freyjacalculated.tsv",
            sep = "\t",
            row.names = FALSE)

