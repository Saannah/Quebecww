# #!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# ## first argument is the replicate name
# ## this script was run on the pooled cut offs so they are hard coded later in the script
rm(list = ls())
library(stringr)

table_names = read.delim("/lustre07/scratch/sanna/ww_sequences/tabular_names", header = FALSE)
table_names = table_names[, 1]

outdir = "/home/sanna/scratch/ww_sequences/"
read_dir = "/home/sanna/scratch/ww_sequences/vcfs_to_use_tabular/"

is_acceptable_snp_strict = function(freq, fixation_freq, min_alt_dp){
  if(freq >= fixation_freq){
    if(freq <= (1 - fixation_freq)){
      return(TRUE)
    }
  }
  return(FALSE)
}

is_acceptable_snp_relaxed = function(freq, sample_alt_dp, sample_ref_dp, fixation_freq, min_alt_dp){
  if(freq >= fixation_freq){
    if(freq <= (1 - fixation_freq)){
      return(TRUE)
    }
  }
  
  ## lowfreq supported by 25 alt reads
  if(freq < fixation_freq){
    if(sample_alt_dp >= min_alt_dp){
      return(TRUE)
    }
  }
  
  ## high freq supported by 25 ref reads
  if(freq > (1 - fixation_freq)){
    if(sample_ref_dp >= min_alt_dp){
      return(TRUE)
    }
  }
  
  return(FALSE)
}

out = data.frame(filename = vector(),
                 nseg_strict = vector(),
                 nonunique_strict = vector(),
                 noindels_strict = vector(),
                 nseg_relaxed = vector(),
                 nonunique_relaxed = vector(),
                 noindels_relaxed = vector())

min_depth = 50
fixation_freq = 0.015
min_alt_dp = 25


for(filename in table_names){
  print(filename)
  if(str_detect(filename, "NoSegSitesCalled")){
    out = rbind(out,
                data.frame(filename = filename,
                           nseg_strict = 0,
                           nonunique_strict = 0,
                           noindels_strict = 0,
                           nseg_relaxed = 0,
                           nonunique_relaxed = 0,
                           noindels_relaxed = 0)
    )
    next
    
  }
  
  vcf = read.delim(paste0(read_dir, filename))
  
  
  cleaned_vcf_strict = vcf[which(unlist(mapply(is_acceptable_snp_strict, vcf$ALT_FREQ, MoreArgs = list(fixation_freq = fixation_freq, min_alt_dp=min_alt_dp)))), ]
  cleaned_vcf_strict = cleaned_vcf_strict[which(cleaned_vcf_strict$TOTAL_DP >= 50), ]
  nseg_strict = length(unique(cleaned_vcf_strict$POS))
  nonunique_strict = nrow(cleaned_vcf_strict)
  
  ## remove indels from strict
  cleaned_vcf_strict = cleaned_vcf_strict[which(nchar(cleaned_vcf_strict$REF) == 1), ]
  cleaned_vcf_strict = cleaned_vcf_strict[which(nchar(cleaned_vcf_strict$ALT) == 1), ]
  noindels_strict = length(unique(cleaned_vcf_strict$POS))
  
  
  cleaned_vcf_relaxed = vcf[which(unlist(mapply(is_acceptable_snp_relaxed, vcf$ALT_FREQ, vcf$ALT_DP, vcf$REF_DP, MoreArgs = list(fixation_freq = fixation_freq, min_alt_dp=min_alt_dp)))), ]
  cleaned_vcf_relaxed = cleaned_vcf_relaxed[which(cleaned_vcf_relaxed$TOTAL_DP >= 50), ]
  nseg_relaxed = length(unique(cleaned_vcf_relaxed$POS))
  nonunique_relaxed = nrow(cleaned_vcf_relaxed)
  
  ##remove indels from relaxed
  cleaned_vcf_relaxed = cleaned_vcf_relaxed[which(nchar(cleaned_vcf_relaxed$REF) == 1), ]
  cleaned_vcf_relaxed = cleaned_vcf_relaxed[which(nchar(cleaned_vcf_relaxed$ALT) == 1), ]
  noindels_relaxed = length(unique(cleaned_vcf_relaxed$POS))
  
  
  out = rbind(out,
              data.frame(filename = filename,
                         nseg_strict = nseg_strict,
                         nonunique_strict = nonunique_strict,
                         noindels_strict = noindels_strict,
                         nseg_relaxed = nseg_relaxed,
                         nonunique_relaxed = nonunique_relaxed,
                         noindels_relaxed = noindels_relaxed)
  )
  
}


write.table(out, 
            file = "REDO_01_output_nseg_estimates_corrected.tsv",
            row.names = FALSE, 
            sep = "\t")
