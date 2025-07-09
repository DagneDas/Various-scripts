###########################################################################################
#                              Script name: 'merge_betas_third_trimester_PART_3_14.12.22.R'
#                                                     Part 3
# This script merges the beta values of Ubiquitous, placenta-specific and KCNQ-BWS DMRs (Dave's data - BWS placenta tumour tissues) with the beta values from the paper (third trimester placenta samples- cell types).
# Input: BWS_DMRs_with_beta_values_x__29_11_22.csv, ubiquitous_DMRs_with_beta_values_x_29_11_22.csv', placenta-specific_DMRs_with_beta_values_x_29_11_22.csv - these tables contain beta values from each cell type for each placenta. These tables have to be manually edited(many the same duplicated columns have to be deleted before the use).
# Output: third_trimester_placenta_betas_merged.csv(KCNQ-BWS), ubiquitous_DMRS_third_trimester_placenta_betas_merged.csv, placenta_DMRS_third_trimester_placenta_betas_merged.csv - these tables have  all paper data.(2) Then you have to edit the files manually- delete many replicated columns (due to merging multiple times). (3) Import edited tables and merge them with Dave's five samples. (4) saves three tables: 'BWS_KCNQ_DMR_beta_comparison_table_05.12.22.csv', 'ubiquitous_DMRS_BWS_KCNQ_DMR_beta_comparison_table_15.12.22.csv', 'Placenta_DMR_final_table_KCNQ_DMR_beta_comparison_table_15.12.22.csv'
# * there is some code to change the orientation of those tables (transpose) to have the same order of cells across all placenta samples. 
# last updated on 14.12.22
# 
##########################################################################################


# import libraries 
library(xlsx)



library(dplyr)
library(ggpubr)
library(tidyverse)
library(stringr)

#set a working directory 
path_name <- 'D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/'
setwd(path_name)

# the IDs of all placentas from the paper
Placenta_IDs <- read.table(file = 'Third_trimester_placenta_IDS.txt', header = FALSE)

#set a working directory 
#path_name2 <- 'D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/PM374'
#setwd(path_name2)

# imports the beta values of the last placenta sample from the text file
#Placenta_374 <- read.csv('BWS_DMRs_with_beta_values_PM374_29_11_22.csv', header = TRUE)

# function
# need to provide placenta names/IDs
merge_placenta_betas <- function(placenta_1, placenta_2){

  # Placenta ID stored in a variable
  folder_name <- placenta_1
  
  # path into one of the placentas' folders
  path_name_temp <- paste('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/', folder_name, sep = "")
    
  setwd(path_name_temp)
    
  # import placenta table with betas
  # manually edited the tables deleted the same columns, left only beta values
  table1 <- read.csv(paste('BWS_DMRs_with_beta_values_', folder_name,'_29_11_22.csv', sep = ""), header = TRUE)
  
  # Stored other placenta ID into a variable
  folder_name <- placenta_2
  
  # directs the path to the folder of other placenta
  path_name_temp <- paste('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/', folder_name, sep = "")
    
  setwd(path_name_temp)
  
  # imports beta values of the other placenta sample
  # manually edited the tables deleted the same columns, left only beta values
  table2 <- read.csv(paste('BWS_DMRs_with_beta_values_', folder_name,'_29_11_22.csv', sep=""), header = TRUE)
  
  # merges both placenta tables  
  merged_table <- merge(table1, table2, by.x = 'IlmnID', by.y = 'IlmnID')
  
  # returns a table with merged betas of two placentas
  return(merged_table)
}

# merge all table of each placenta to a one big table
merged_table1 <- merge_placenta_betas(placenta_1 = 'PM365', placenta_2 = 'PM376')
merged_table2 <- merge_placenta_betas(placenta_1 = 'PM370', placenta_2 = 'PM368')  
merged_table3 <- merge_placenta_betas(placenta_1 = 'PM367', placenta_2 = 'PM371') 
merged_table4 <- merge_placenta_betas(placenta_1 = 'PM378', placenta_2 = 'PM369') 
merged_table5 <- merge_placenta_betas(placenta_1 = 'PM359', placenta_2 = 'PM379') 
merged_table6 <- merge_placenta_betas(placenta_1 = 'PM377', placenta_2 = 'PM381') 
merged_table7 <- merge_placenta_betas(placenta_1 = 'PM362', placenta_2 = 'PM372') 
merged_table8 <- merge_placenta_betas(placenta_1 = 'PM366', placenta_2 = 'PM375') 
merged_table9 <- merge_placenta_betas(placenta_1 = 'PM373', placenta_2 = 'PM374') 

merged_table <- merge(merged_table1, merged_table2, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table3, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table4, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table5, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table6, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table7, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table8, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table9, by.x = 'IlmnID', by.y = 'IlmnID')


setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')
# export the table and delete many columns, change some column names
#write.csv(merged_table, 'third_trimester_placenta_betas_merged.csv')


# import a fixed table 
setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')

fixed_table <- read.csv('third_trimester_placenta_betas_merged.csv', header = TRUE)

# update columns of these samples, were not downloaded properly previously
setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/PM359')

PM359_GSM4831954 = read.table('GSM4831954-118958.txt', header = TRUE)

setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/PM362')

PM362_GSM4831965 = read.table('GSM4831965-118969.txt', header = TRUE)

# merge a final table with tables of files that had to be replaced 
updated_fixed_table <- merge(fixed_table, PM359_GSM4831954, by.x = 'IlmnID', by.y = 'ID_REF', all.x = TRUE)
updated_fixed_table$Beta_value_GSM4831954.118958 <- updated_fixed_table$VALUE

updated_fixed_table <- merge(updated_fixed_table, PM362_GSM4831965, by.x = 'IlmnID', by.y = 'ID_REF', all.x = TRUE)
updated_fixed_table$Beta_value_GSM4831965.118969 <- updated_fixed_table$VALUE.y

# write a good table paper data
setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')
#write.csv(updated_fixed_table, 'third_trimester_placenta_betas_merged.csv')

## import data from Dave's methylationEpic array

# set a working directory to the folder there Dave's data is stored 
setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/Daves BWS data epic 850k array/BWS_analysis")

# imports Dave's data with 5 placenta samples
Daves_data <- read.csv('BWS_placentas_beta_values_03.12.22.csv', header = TRUE)

# merge paper data with Daves data
Final_table <- merge(updated_fixed_table, Daves_data,  by.x = 'IlmnID', by.y = 'IlmnID', all.x = TRUE)

# sort the order of columns
# first columns with general info, then columns with Dave's data, finally columns with paper data
BWS_final_table <- data.frame(Final_table[,1:7], Final_table[,105:109], Final_table[,8:104] )

# save final table
setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ")

write.csv(BWS_final_table, 'BWS_KCNQ_DMR_beta_comparison_table_05.12.22.csv')


#############################################################################################################################

                                                                      # Ubiquitous DMRs


#set a working directory 
path_name <- 'D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/'
setwd(path_name)

# the IDs of all placentas from the paper
Placenta_IDs <- read.table(file = 'Third_trimester_placenta_IDS.txt', header = FALSE)


merge_placenta_betas <- function(placenta_1, placenta_2){
  
  # Placenta ID stored in a variable
  folder_name <- placenta_1
  
  # path into one of the placentas' folders
  path_name_temp <- paste('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/', folder_name, sep = "")
  
  setwd(path_name_temp)
  
  # import placenta table with betas
  table1 <- read.csv(paste('ubiquitous_DMRs_with_beta_values_', folder_name,'_29_11_22.csv', sep = ""), header = TRUE)
  
  # Stored other placenta ID into a variable
  folder_name <- placenta_2
  
  # directs the path to the folder of other placenta
  path_name_temp <- paste('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/', folder_name, sep = "")
  
  setwd(path_name_temp)
  
  # imports beta values of the other placenta sample
  table2 <- read.csv(paste('ubiquitous_DMRs_with_beta_values_', folder_name,'_29_11_22.csv', sep=""), header = TRUE)
  
  # merges both placenta tables  
  merged_table <- merge(table1, table2, by.x = 'IlmnID', by.y = 'IlmnID')
  
  # returns a table with merged betas of two placentas
  return(merged_table)
}

# merge all table of each placenta to a one big table
merged_table1 <- merge_placenta_betas(placenta_1 = 'PM365', placenta_2 = 'PM376')
merged_table2 <- merge_placenta_betas(placenta_1 = 'PM370', placenta_2 = 'PM368')  
merged_table3 <- merge_placenta_betas(placenta_1 = 'PM367', placenta_2 = 'PM371') 
merged_table4 <- merge_placenta_betas(placenta_1 = 'PM378', placenta_2 = 'PM369') 
merged_table5 <- merge_placenta_betas(placenta_1 = 'PM359', placenta_2 = 'PM379') 
merged_table6 <- merge_placenta_betas(placenta_1 = 'PM377', placenta_2 = 'PM381') 
merged_table7 <- merge_placenta_betas(placenta_1 = 'PM362', placenta_2 = 'PM372') 
merged_table8 <- merge_placenta_betas(placenta_1 = 'PM366', placenta_2 = 'PM375') 
merged_table9 <- merge_placenta_betas(placenta_1 = 'PM373', placenta_2 = 'PM374') 

# after some time provides an error, due to not unique column names
# can ignore it, tables are merged. 
merged_table <- merge(merged_table1, merged_table2, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table3, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table4, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table5, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table6, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table7, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table8, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table9, by.x = 'IlmnID', by.y = 'IlmnID')


setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')
# export the table and delete many columns, change some column names
#write.csv(merged_table, 'ubiquitous_DMRS_third_trimester_placenta_betas_merged.csv')


# import a fixed table 
setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')

fixed_table <- read.csv('ubiquitous_DMRS_third_trimester_placenta_betas_merged.csv', header = TRUE)

# Fixed these columns (skip this step)
# update columns of these samples, were not downloaded properly previously
#setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/PM359')
#PM359_GSM4831954 = read.table('GSM4831954-118958.txt', header = TRUE)
#setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/PM362')
#PM362_GSM4831965 = read.table('GSM4831965-118969.txt', header = TRUE)
# merge a final table with tables of files that had to be replaced 
#updated_fixed_table <- merge(fixed_table, PM359_GSM4831954, by.x = 'IlmnID', by.y = 'ID_REF', all.x = TRUE)
#updated_fixed_table$Beta_value_GSM4831954.118958 <- updated_fixed_table$VALUE
#updated_fixed_table <- merge(updated_fixed_table, PM362_GSM4831965, by.x = 'IlmnID', by.y = 'ID_REF', all.x = TRUE)
#updated_fixed_table$Beta_value_GSM4831965.118969 <- updated_fixed_table$VALUE.y

# start from here 
# write a good table paper data
#setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')
#write.csv(updated_fixed_table, 'ubiquitous_DMRS_third_trimester_placenta_betas_merged.csv')

## import data from Dave's methylationEpic array

# set a working directory to the folder there Dave's data is stored 
setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/Daves BWS data epic 850k array/BWS_analysis")

# imports Dave's data with 5 placenta samples
Daves_data <- read.csv('BWS_placentas_beta_values_03.12.22.csv', header = TRUE)

# merge paper data with Daves data
Final_table <- merge(fixed_table, Daves_data,  by.x = 'IlmnID', by.y = 'IlmnID', all.x = TRUE)

# sort the order of columns
# first columns with general info, then columns with Dave's data, finally columns with paper data
ubiquitous_DMRS_final_table <- data.frame(Final_table[,1:7], Final_table[,103:107], Final_table[,8:102] )

# save final table
setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ")

#write.csv(ubiquitous_DMRS_final_table, 'ubiquitous_DMRS_BWS_KCNQ_DMR_beta_comparison_table_15.12.22.csv')



#############################################################################################################################

#                                                     placenta DMRs

#set a working directory 
path_name <- 'D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/'
setwd(path_name)

# the IDs of all placentas from the paper
Placenta_IDs <- read.table(file = 'Third_trimester_placenta_IDS.txt', header = FALSE)

#set a working directory 
#path_name2 <- 'D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/PM374'
#setwd(path_name2)

# imports the beta values of the last placenta sample from the text file
#Placenta_374 <- read.csv('BWS_DMRs_with_beta_values_PM374_29_11_22.csv', header = TRUE)

# function
# need to provide placenta names/IDs
merge_placenta_betas <- function(placenta_1, placenta_2){
  
  # Placenta ID stored in a variable
  folder_name <- placenta_1
  
  # path into one of the placentas' folders
  path_name_temp <- paste('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/', folder_name, sep = "")
  
  setwd(path_name_temp)
  
  # import placenta table with betas
  # manually edited the tables, deleted the same columns, left only beta values
  table1 <- read.csv(paste('placenta-specific_DMRs_with_beta_values_', folder_name,'_29_11_22.csv', sep = ""), header = TRUE)
  
  # Stored other placenta ID into a variable
  folder_name <- placenta_2
  
  # directs the path to the folder of other placenta
  path_name_temp <- paste('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/', folder_name, sep = "")
  
  setwd(path_name_temp)
  
  # imports beta values of the other placenta sample
  # manually edited the tables deleted the same columns, left only beta values
  table2 <- read.csv(paste('placenta-specific_DMRs_with_beta_values_', folder_name,'_29_11_22.csv', sep=""), header = TRUE)
  
  # merges both placenta tables  
  merged_table <- merge(table1, table2, by.x = 'IlmnID', by.y = 'IlmnID')
  
  # returns a table with merged betas of two placentas
  return(merged_table)
}

# merge all table of each placenta to a one big table
merged_table1 <- merge_placenta_betas(placenta_1 = 'PM365', placenta_2 = 'PM376')
merged_table2 <- merge_placenta_betas(placenta_1 = 'PM370', placenta_2 = 'PM368')  
merged_table3 <- merge_placenta_betas(placenta_1 = 'PM367', placenta_2 = 'PM371') 
merged_table4 <- merge_placenta_betas(placenta_1 = 'PM378', placenta_2 = 'PM369') 
merged_table5 <- merge_placenta_betas(placenta_1 = 'PM359', placenta_2 = 'PM379') 
merged_table6 <- merge_placenta_betas(placenta_1 = 'PM377', placenta_2 = 'PM381') 
merged_table7 <- merge_placenta_betas(placenta_1 = 'PM362', placenta_2 = 'PM372') 
merged_table8 <- merge_placenta_betas(placenta_1 = 'PM366', placenta_2 = 'PM375') 
merged_table9 <- merge_placenta_betas(placenta_1 = 'PM373', placenta_2 = 'PM374') 

merged_table <- merge(merged_table1, merged_table2, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table3, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table4, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table5, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table6, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table7, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table8, by.x = 'IlmnID', by.y = 'IlmnID')
merged_table <- merge(merged_table, merged_table9, by.x = 'IlmnID', by.y = 'IlmnID')


setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')
# export the table and delete many columns, change some column names
#write.csv(merged_table, 'placenta_DMRS_third_trimester_placenta_betas_merged.csv')


# import a fixed table 
setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')

fixed_table <- read.csv('placenta_DMRS_third_trimester_placenta_betas_merged.csv', header = TRUE)

# update columns of these samples, were not downloaded properly previously
#setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/PM359')
#PM359_GSM4831954 = read.table('GSM4831954-118958.txt', header = TRUE)
#setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/PM362')
#PM362_GSM4831965 = read.table('GSM4831965-118969.txt', header = TRUE)
# merge a final table with tables of files that had to be replaced 
#updated_fixed_table <- merge(fixed_table, PM359_GSM4831954, by.x = 'IlmnID', by.y = 'ID_REF', all.x = TRUE)
#updated_fixed_table$Beta_value_GSM4831954.118958 <- updated_fixed_table$VALUE
#updated_fixed_table <- merge(updated_fixed_table, PM362_GSM4831965, by.x = 'IlmnID', by.y = 'ID_REF', all.x = TRUE)
#updated_fixed_table$Beta_value_GSM4831965.118969 <- updated_fixed_table$VALUE.y

# write a good table paper data
#setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')
#write.csv(updated_fixed_table, 'third_trimester_placenta_betas_merged.csv')

## import data from Dave's methylationEpic array

# set a working directory to the folder there Dave's data is stored 
setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/Daves BWS data epic 850k array/BWS_analysis")

# imports Dave's data with 5 placenta samples
Daves_data <- read.csv('BWS_placentas_beta_values_03.12.22.csv', header = TRUE)

# merge paper data with Daves data
Final_table <- merge(fixed_table, Daves_data,  by.x = 'IlmnID', by.y = 'IlmnID', all.x = TRUE)

# sort the order of columns
# first columns with general info, then columns with Dave's data, finally columns with paper data
Placenta_DMR_final_table <- data.frame(Final_table[,1:7], Final_table[,103:107], Final_table[,8:102] )

# save final table
setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ")

#write.csv(Placenta_DMR_final_table, 'Placenta_DMR_final_table_KCNQ_DMR_beta_comparison_table_15.12.22.csv')





#############################################################################################################################

#                                                     Sort cell types in the table

#############################################################################################################################

## import the final beta table with cell names to reorder columns
#setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ")

#Final_table <- read.table("BWS_KCNQ_DMR_beta_comparison_table_05.12.22.csv", header =TRUE,  sep = ",")

## reverse axes, export, use excel to sort based on placenta ID and cell types
#Final_table_transpose <- t(Final_table)
#write.csv(Final_table_transpose, 'transpose_BWS_KCNQ_DMR_beta_comparison_table_15.12.22.csv')

##Final_table_inverted <- setNames(data.frame(t(data[ , - 1])), data[ , 1])  # Transpose data
##data_t  

## import sorted table, reverse the axes. 
#Final_table2 <- read.xlsx("transpose_BWS_KCNQ_DMR_beta_comparison_table_15.12.22.xlsx", sheetIndex = 1 , header =TRUE)
#transpose <-  t(Final_table2)

#write.xlsx(transpose, 'BWS_KCNQ_DMR_beta_comparison_final_table_15.12.22.xlsx', col.names = FALSE)

#############################################################################################################################

#                                                          ubiquitous DMRs


# import the final beta table with cell names to reorder columns
#setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ")

#Final_table <- read.csv("ubiquitous_DMRS_BWS_KCNQ_DMR_beta_comparison_table_15.12.22.csv", header =TRUE,  sep = ",")

# reverse axes, export, use excel to sort based on placenta ID and cell types
#Final_table_transpose <- t(Final_table)
#write.csv(Final_table_transpose, 'transpose_ubiquitous_DMRS_KCNQ_DMR_beta_comparison_table_15.12.22.csv')

#Final_table_inverted <- setNames(data.frame(t(data[ , - 1])), data[ , 1])  # Transpose data
#data_t  

# import sorted table, reverse the axes. 
#Final_table2 <- read.xlsx("transpose_ubiquitous_DMRS_KCNQ_DMR_beta_comparison_table_15.12.22.xlsx", sheetIndex = 1 , header =TRUE)
#transpose <-  t(Final_table2)

#write.xlsx(transpose, 'Ubiquitous_DMRs_KCNQ_DMR_beta_comparison_final_table_15.12.22.xlsx', col.names = FALSE)


#############################################################################################################################

#                                                          placenta DMRs


# import the final beta table with cell names to reorder columns
#setwd("D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ")

#Final_table <- read.csv("Placenta_DMR_final_table_KCNQ_DMR_beta_comparison_table_15.12.22.csv", header =TRUE,  sep = ",")

# reverse axes, export, use excel to sort based on placenta ID and cell types
#Final_table_transpose <- t(Final_table)
#write.csv(Final_table_transpose, 'transpose_Placenta_DMRS_KCNQ_DMR_beta_comparison_table_15.12.22.csv')

#Final_table_inverted <- setNames(data.frame(t(data[ , - 1])), data[ , 1])  # Transpose data
#data_t  

# import sorted table, reverse the axes.
#Final_table2 <- read.xlsx("transpose_Placenta_DMRS_KCNQ_DMR_beta_comparison_table_15.12.22.xlsx", sheetIndex = 1 , header =TRUE)
#transpose <-  t(Final_table2)

#write.xlsx(transpose, 'Placenta_DMRs_KCNQ_DMR_beta_comparison_final_table_15.12.22.xlsx', col.names = FALSE)
