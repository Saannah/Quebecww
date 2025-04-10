rm(list = ls())

#### add qPCR to the dataset


### load the viral copy table
vcopy = read.csv("/users/sana/Documents/CANMOD/9-RD-viral_copy_data/donnees-vigie-eaux-usees_archives_19072023.csv")
unique(vcopy$site_code)

sites = c("GTN_01",
          "LVL_01",
          "LVL_02",
          "LVL_05",
          "MTL_01",
          "MTL_02",
          "QC_01",
          "QC_02",
          "SI_01",
          "RIM_01")


df = data.frame()

for(site in sites){
  df = rbind(df,
             vcopy[which(vcopy$site_code == site), ])
}

## remove NA rows
df = df[which(!is.na(df$SRAS_ajust_flowPop_c)), ]
df$identifier = paste0(df$site_code, 
                       "_",
                       df$Date)

### load case data for 6cities and the associated metrics dataframe:

metrics = read.delim("/users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases.tsv")
metrics$city_id_new = metrics$City_ID
metrics$city_id_new[which(metrics$City_ID == "Sept")] = "SI"
metrics$city_id_new[which(metrics$City_ID == "Riki")] = "RIM"
metrics$identifier = paste0(metrics$city_id_new,
                            "_0",
                            metrics$Site_Number,
                            "_",
                            metrics$date)



#### some of the days for ww data do not have the qpcr data associated with it
length(which(!metrics$identifier %in% df$identifier))

metrics = metrics[which(metrics$identifier %in% df$identifier), ]
metrics$vcopy_per_person = df$SRAS_ajust_flowPop_c[match(metrics$identifier, df$identifier)]
metrics$vcopy_per_100k = df$SRAS_ajust_flowPop2_c[match(metrics$identifier, df$identifier)]

write.table(metrics,
            file = "/Users/sana/Documents/CANMOD/0-reference_files_redo_march2025/Analysis_Report-Slim-2025-03-06_merged_w_nsegs_w_cases_w_qpcr.tsv",
            sep = "\t",
            row.names = FALSE)



################ I omitted the following bit of code which used to multiply the viral copy by population of the station
###### because station population is fixed for every site, this will not affect the time series forecast for every site

# metrics$site_code = paste0(metrics$city_id_new,
#                            "_0",
#                            metrics$site_number)
# 
# ### load the city population dataframe
# city_population_table = read.delim("/users/sana/Documents/CANMOD/9-RD-viral_copy_data/CITY_POPULATION_TABLE.tsv")
# 
# metrics$station_population = city_population_table$Population[match(metrics$site_code, 
#                                                                     city_population_table$City_identifier)]
# 
# ### transform the vcopy to actual gene copy numbers
# metrics$actual_vcopy = metrics$vcopy * metrics$station_population 
# 
# 
# 
# dim(metrics)




