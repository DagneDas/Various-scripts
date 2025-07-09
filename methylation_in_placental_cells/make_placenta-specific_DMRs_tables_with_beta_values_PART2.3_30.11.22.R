##############################################################################################################################################

#                                                       Script name: 'make_placenta-specific_DMRs_tables_with_beta_values_PART2.3_30.11.22.R'
#                                                                                   Part 2.3
# last updated 30.11.22
# Input: Have two provide 'placenta_DMRs_probes_with_CpG_probes_24_11_22.csv' - the list of my probes;
# Input:'Placenta_folder_list_5_datasets.txt placenta list with all folder names- for 5 cell types
# Input: or 'Placenta_folder_list_6_datasets.txt' - for 6 cell types
# output: gives a large table with 5  or 6 cell-types, their beta values
# output: returns a large list for each placenta sample with all cell-types


##############################################################################################################################################

# import libraries 
library(xlsx)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(stringr)

# set a working directory 
setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ')

# imports datasets with my probes
#kcnq_bws_probes <- read.csv('BWS_DMRs_with_CpG_probes_24_11_22.csv', header = TRUE)
#ubiquitous_DMRs_probes <- read.csv('ubiquitous_DMRs_with_CpG_probes_24_11_22.csv', header = TRUE)
placenta_DMRs_probes <- read.csv('placenta_DMRs_probes_with_CpG_probes_24_11_22.csv', header = TRUE)


# set a working directory for the placenta paper with all folders
setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester')

# import a list with placenta folders that have 5 cell types 
placenta_sample_list_5_datasets <- read.table('Placenta_folder_list_5_datasets.txt')

# import a list with placenta folders that have 6 cell types 
placenta_sample_list_6_datasets <- read.table('Placenta_folder_list_6_datasets.txt')


# function takes a name from the placenta list and directs a working directory to one placenta folder. 
# paper datasets are imported and cell-type specific beta values are extracted. 
# all tables are fused into a single list.
get_beta_values <- function(i) {
  
  # print the placenta sample name
  print(i)
  
  # the name of a working directory 
  folder_name <-  i
  
  # a full path to the working directory
  path_name <- paste('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/GSE159526_cell_line_methylation_placentas/Third trimester/', folder_name, sep="")
  
  # set a working directory - one placenta folder
  setwd(path_name)
  
  # imports a list, which contains the names of datasets of one placenta sample
  file_list_in_folder <- read.table('file_list.txt')
  
  # creates a new table to store data
  placenta_table <- data.frame()
  
  # function to extract placenta beta values
  fuse_placenta_tables <- function(y){
    
    # name of one datasets - methylation values of one cell type
    name <- y
    
    # full file name
    file_name <- paste(name,".txt", sep='')
    
    # imports table with beta values of one cell type
    probe_table <- read.table(file_name)
    
    # obtains the index of the last row in the table
    last_row_number <- as.numeric(nrow(probe_table))
    
    column_name1 = paste("Probe_ID_", name , sep = "")
    column_name2 = paste("Beta_value_", name, sep = "")
    
    # sorts a table - excludes some unwanted rows from the top 
    probe_table2 <- data.frame(column_name1 = probe_table[2:last_row_number,1], column_name2 = probe_table[2:last_row_number,2])
    
    #renames table's column names
    colnames(probe_table2)[1] = column_name1
    colnames(probe_table2)[2] = column_name2
    
    #removes not required data
    rm(probe_table)
    rm(last_row_number)
    
    # merges my probe IDs with beta values from the paper
    placenta_table <- merge(placenta_DMRs_probes, probe_table2, by.x = 'IlmnID', by.y = column_name1, all.x = TRUE)
    
    # returns a table that will be saves as a list
    return(placenta_table)
  }
  
  # calls a function and loops through all cell types, cell-type specific beta values are stored in a big list. This list has 5 list for 5 cell types.
  # to reach one list type  BWS_big_table[[X]]
  placenta_big_table <- lapply(file_list_in_folder[,1],  fuse_placenta_tables) 
  
  # fuses all list into one big table with cell-specific beta values
  placenta <- placenta_big_table[[1]] %>% 
    left_join(placenta_big_table[[2]], by = "IlmnID") %>% 
    left_join(placenta_big_table[[3]], by = "IlmnID") %>% 
    left_join(placenta_big_table[[4]], by = "IlmnID") %>%
    left_join(placenta_big_table[[5]], by = "IlmnID") 
  
  
  # makes a name for the file
  data_file_name = paste('placenta-specific_DMRs_with_beta_values_', folder_name, '_29_11_22.csv', sep = '')
  
  # writes a table for one placenta
  write.csv(placenta, file = data_file_name, row.names = FALSE)
  
  # saving an object in RData format
  save(placenta_big_table, file = paste(folder_name,"data.RData", sep=''))
  
}

# loops through all placenta folders
lapply(placenta_sample_list_5_datasets[,1], get_beta_values)


