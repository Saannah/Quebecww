rm(list = ls())

## laod the new analysis report with the newly calculated nsegs 
data = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs.tsv")

################################# This portion of the code adds case counts from the cumulative cases file
## add case city by city
cities = c("MTL", 
           "GTN",
           "LVL",
           "QC",
           "Sept",
           "Riki")

cases_estimates_combined = data.frame()

for(c in cities){
  est_df = data[which(data$City_ID == c), ]
  case_df = read.delim(paste0("/users/sana/Documents/CANMOD/1-RD-case_data/new_case_data_from_website/",c, "_perday.tsv"))
  est_df = est_df[which(est_df$Reference_Date %in% case_df$date), ]
  est_df = est_df[order(est_df$Reference_Date), ]
  est_df$daily_count = case_df$cases[match(est_df$Reference_Date, case_df$date)]
  cases_estimates_combined = rbind(cases_estimates_combined,
                                   est_df)
}

# everything = read.delim("/users/sana/Documents/CANMOD/0-reference_files/cases_estimates_combined_6cities_archdownload_qpcr_wFreyja_slimHQ_posrate_weektotal_hosp_newnseg.tsv")
# length(which(everything$Sequence %in% data$Sequence))
# everything_mtl = everything[which(everything$city_id == c), ]
# everything_mtl = everything_mtl[order(everything_mtl$date), ]

############################################################
##### dictionary for the sites
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

everything = data

all = data.frame()
for(s in 1:nrow(sites_df)){
  city_id = sites_df$city_id[s]
  site_number = sites_df$site_number[s]
  
  city_data = everything[which(everything$City_ID == city_id), ]
  city_data = city_data[which(city_data$Site_Number == site_number), ]
  
  ## check for duplicated dates in city data
  if(length(which(duplicated(city_data$date))) > 0){
    print("site has duplicated dates:")
    print(paste0(city_id, "-", site_number))
    city_data = city_data[-which(duplicated(city_data$date)), ]
  }
  
  ### pull up positivity rate data for this site
  pos_rate_data = read.delim(paste0("/users/sana/Documents/CANMOD/1-RD-case_data/positive_and_negative_case_data_downloaded_from_website/",city_id, "-tests_weeks_numbered.tsv"))
  
  ## take the substring for date column
  pos_rate_data$date = substr(pos_rate_data$Date.de.résultat.du.test, start = 1, stop = 10)
  
  ## check if all dates in city data are in the pos rate dataset
  if(!length(which(city_data$Reference_Date %in% pos_rate_data$date)) == nrow(city_data)){
    #print("some dates in this city dont have case data in posrate dataset")
    stop("some dates in this city dont have case data in posrate dataset")
  }
  
  
  pos_rate_data = pos_rate_data[which(pos_rate_data$date %in% city_data$Reference_Date), ]
  
  colnames(city_data)[colnames(city_data) == 'Reference_Date'] <- 'date'
  all = rbind(all, merge(city_data, pos_rate_data, "date"))
}


#### sanity check: comoare case counts from the two datasets (cumulative and the one with total tests)
for(s in 1:nrow(sites_df)){
  city_id = sites_df$city_id[s]
  site_number = sites_df$site_number[s]
  
  city_data_1 = all[which(all$City_ID == city_id), ]
  city_data_1 = city_data_1[which(city_data_1$Site_Number == site_number), ]
  

  city_data_2 = cases_estimates_combined[which(cases_estimates_combined$City_ID == city_id), ]
  city_data_2 = city_data_2[which(city_data_2$Site_Number == site_number), ]
  
  ### check 1
  if(nrow(city_data_1) != nrow(city_data_2)){
    stop("the two datasets are inconsistent")
  }
  
  ### check actual case counts
  city_data_1 = city_data_1[order(city_data_1$date), ]
  city_data_2 = city_data_2[order(city_data_2$Reference_Date), ]
  
  non_eq = (which(city_data_1$Tests.positifs != city_data_2$daily_count))
  non_eq2 = (which(city_data_2$daily_count != city_data_1$Tests.positifs))
  
  city_data_1_non_eqs = city_data_1[non_eq, ]
  city_data_2_non_eqs = city_data_2[non_eq2, ]
  
  }

######## the cumulative case counts have been removed from INSPQ website
#### there is minor descrepancy between the older dataset that had cumulative cases and the new dataset that still exista

##### I will write only the version with positive tests and negative tests that Ines had sent in March 2024
write.table(all,
            file = "/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases.tsv",
            sep = "\t",
            row.names = FALSE)



