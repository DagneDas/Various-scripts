#######################################################################################################################

#                         Script to obtain CpG IDs from the Illumina manifest file for BWS/KCNQ project: part 1
#
# For this script, I used tables with Ubiquitous and placenta-specific DMRs, new KCNQ or BWS related DMRs
# Output: get new tables for these three groups with CpG probe IDS. These probes should be used to obtain B values from the 'Cell-specific characterization of the placental methylome' paper. https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07186-6
#
# Script created 24.11.22
# Script last updated 24.11.22
# Script name: Obtain_CpG_probe_IDs_BWS-KCNQ_project_PART_1.R
#######################################################################################################################

# install packages
#install.packages('xlsx')
#install.packages('dplyr')
#install.packages('ggpubr')
#install.packages('tidyverse')
#install.packages("stringr")

# import libraries 
library(xlsx)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(stringr)

# set a working directory with files with the manifest file

setwd('D:/NRP DTP (PhD)/Data/3rd year/')

manifest_file = read.csv(file = 'infinium_manifest-file.csv', header = TRUE)


# set a working directory with other files

setwd('D:/NRP DTP (PhD)/Data/3rd year/BWS_related_data_KCNQ/')

# import files
ubiquitous_DMRs = read.xlsx(file = 'Probes_for_BWS_PL&_placenta_DMRs_file_for_analysis.xlsx', sheetName = 'Ubiquitous_DMRs', header = TRUE)
KCNQ_DMRs =  read.xlsx(file = 'Probes_for_BWS_PL&_placenta_DMRs_file_for_analysis.xlsx', sheetName = 'Novel_chr11_DMRs_KCNQ', header = TRUE)
placenta_DMRs =  read.xlsx(file = 'Probes_for_BWS_PL&_placenta_DMRs_file_for_analysis.xlsx', sheetName = 'placenta_specific_DMRs', header = TRUE)

# rename columns for clarity
manifest_file$hg19_start = manifest_file$MAPINFO
manifest_file$hg19_end = manifest_file$MAPINFO


# creates a vector containing 1
addition_vector = rep(1, times = length(manifest_file$MAPINFO))

# creates the hg19 genome end coordinates 
manifest_file$hg19_end = manifest_file$hg19_end + addition_vector
head(manifest_file)

# create a smaller manifest table 

smaller_manifest_table = data.frame('IlmnID' = manifest_file$IlmnID, 'Name'= manifest_file$Name, 'CHR'= manifest_file$CHR, 'hg19_start' = manifest_file$hg19_start, 'hg19_end' = manifest_file$hg19_end)

head(smaller_manifest_table)

# filer out rows with sex chromosomes
smaller_manifest_table <- smaller_manifest_table %>% filter(CHR != 'X')
smaller_manifest_table <- smaller_manifest_table %>% filter(CHR != 'Y')

# check if sex chromosomes were filtered 
y <- smaller_manifest_table %>% filter(CHR == 'Y')
x <- smaller_manifest_table %>% filter(CHR == 'X')

rm(x)
rm(y)

# remove characters from numbers - 'chromosome' column

# KCNQ_file
for (i in 1:6){
  x <- gsub('[chr]', '', KCNQ_DMRs$GRCh37_hg19_chr[i])
  KCNQ_DMRs$GRCh37_hg19_chr[i] <- x
}

# placenta_DMR file
for (i in 1:length(placenta_DMRs$GRCh37_hg19_chr)){
  x <- gsub('[chr]', '', placenta_DMRs$GRCh37_hg19_chr[i])
  placenta_DMRs$GRCh37_hg19_chr[i] <- x
}

# ubiquitous_DMR file
for (i in 1:length(ubiquitous_DMRs$GRCh37_hg19_chr)){
  x <- gsub('[chr]', '', ubiquitous_DMRs$GRCh37_hg19_chr[i])
  ubiquitous_DMRs$GRCh37_hg19_chr[i] <- x
}


# filter CpG probes - function
find_probe_IDs <- function(DMR, name){
  
  # creates an empty table
  DMR_table = data.frame()
  
  # a loop takes a single DMR region and looks for CpG probes at that DMR
  # at the end of the loop all DMRs with belonging CpG probes are fused in a single table
  for (i in 1:length(DMR$GRCh37_hg19_chr)){
    
    #selects one DMR row by row
    single_DMR = DMR[i,]
    
    #filter rows that are have the same chromosome as the DMR table
    single_DMR_chromosome = smaller_manifest_table %>% filter(CHR == single_DMR$GRCh37_hg19_chr)
    
    #filters out those rows that have CpG probes in the DMR region
    single_DMR_probes = single_DMR_chromosome %>% filter(hg19_start >= single_DMR$GRCh37_hg19_start & hg19_end <= single_DMR$GRCh37_hg19_end)
    
    #creates a vector the same length as the table with CpG probes. Vector contains the DMR name
    DMR_name_vector <- rep(single_DMR[1,1], times = length(single_DMR_probes$hg19_start))
    
    #In the table with CpG probes a new column is created with the DMR name, to know which probes belong to which DMRs.
    single_DMR_probes$name <- DMR_name_vector
    
    #All DMRs with their probes are fused in one table 
    DMR_table <- rbind(DMR_table, single_DMR_probes)
  }
  
  return(DMR_table) # returns a single table with DMRs and CpG probes
}

#table with CpG probes
KCNQ_BWS_DMRs <- find_probe_IDs(DMR = KCNQ_DMRs, name = 'Imprinted_DMR')
ubiquitous_DMRs_probes<- find_probe_IDs(DMR = ubiquitous_DMRs, name = 'Imprinted_DMR')
placenta_DMRs_probes <- find_probe_IDs(DMR = placenta_DMRs, name = 'Gene_name')



# finish by adding a column with the DMR coordinates to the DMR probe table
KCNQ_DMRs_smaller = data.frame('Imprinted_DMR' = KCNQ_DMRs$Imprinted_DMR, 'DMR_position_GRCh37_hg19' = KCNQ_DMRs$Position_GRCh37_hg19)

KCNQ_BWS_DMRs2 <- merge(KCNQ_BWS_DMRs, KCNQ_DMRs_smaller, by.x = 'name' , by.y = 'Imprinted_DMR')

#remove tables that are not needed
rm(KCNQ_BWS_DMRs)
rm(KCNQ_DMRs_smaller)



# finish by adding a column with the DMR coordinates to the DMR probe table
ubiquitous_DMRs_smaller = data.frame('Imprinted_DMR' = ubiquitous_DMRs$Imprinted_DMR, 'DMR_position_GRCh37_hg19'   =ubiquitous_DMRs$Position_GRCh37_hg19)

ubiquitous_DMRs_probes2 <- merge(ubiquitous_DMRs_probes, ubiquitous_DMRs_smaller, by.x = 'name' , by.y = 'Imprinted_DMR')

#remove tables that are not needed
rm(ubiquitous_DMRs_probes)
rm(ubiquitous_DMRs_smaller)



# finish by adding a column with the DMR coordinates to the DMR probe table
placenta_DMRs_smaller = data.frame('Gene_name' = placenta_DMRs$Gene_name, 'DMR_position_GRCh37_hg19' = placenta_DMRs$Position_GRCh37_hg19)

placenta_DMRs_probes2 <- merge(placenta_DMRs_probes, placenta_DMRs_smaller, by.x = 'name' , by.y = 'Gene_name')

#remove tables that are not needed
rm(placenta_DMRs_probes)
rm(placenta_DMRs_smaller)

#write tables with CpG probes
#write.csv(KCNQ_BWS_DMRs2, file = 'BWS_DMRs_with_CpG_probes_24_11_22.csv', row.names = FALSE)
#write.csv(ubiquitous_DMRs_probes2, file = 'ubiquitous_DMRs_with_CpG_probes_24_11_22.csv', row.names = FALSE)
#write.csv(placenta_DMRs_probes2, file = 'placenta_DMRs_probes_with_CpG_probes_24_11_22.csv', row.names = FALSE)
