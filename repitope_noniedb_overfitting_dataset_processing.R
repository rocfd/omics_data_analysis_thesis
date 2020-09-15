## NON-IEDB dataset obtained from Repitope package
## Overfitting control
## Data processing

## TCR Contact Region Extraction
setwd("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control")

## libraries
library(plyr)
library(dplyr)
library(gridExtra)
library(formattable)

## Upload Repitope dataset
filename <- "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200829_IEDB_dataset_processing/repitope_full_70543.csv"
repitope_full <- read.table(filename, sep = ",", header = T, skip = 0, fill = T, stringsAsFactors = T)
colnames(repitope_full)
glimpse(repitope_full)
dim(repitope_full) #70543

#### Filter Repitope dataset for IEDB database, human as host and MHC-I epitopes.
## Filter source column for contains IEDB
repitope_full_nonIEDB <- filter(repitope_full, !grepl('IEDB', Source)) #8104
dim(repitope_full_nonIEDB)
table(repitope_full_nonIEDB$Immunogenicity) # Neg 5505 # Pos 2599
write.table(repitope_full_nonIEDB,"repitope_nonIEDB_overfitting_8104.csv",quote=F,sep=",")
View(repitope_full_nonIEDB)

## Filter host column for contains homo sapiens
table(repitope_full_nonIEDB$Host_Organism) ## not indicated in non-IEDB data ## concerning ## keep all
#repitope_full_nonIEDB_human <- filter(repitope_full_nonIEDB, grepl('Homo sapiens', Host_Organism)) #41318
#dim(repitope_full_nonIEDB_human)
#write.table(repitope_full_nonIEDB,"repitope_IEDB_human_41318.csv",quote=F,sep=",")

## Filter MHCType column for MHC-I
repitope_full_nonIEDB_MHCI <- filter(repitope_full_nonIEDB, MHCType == "MHC-I") #11095
dim(repitope_full_nonIEDB_MHCI)
glimpse(repitope_full_nonIEDB_MHCI)
View(repitope_full_nonIEDB_MHCI)
write.table(repitope_full_nonIEDB_MHCI,"repitope_nonIEDB_mhcI_8104.csv",quote=F,sep=";", dec = ".")


## Filter MHCtype column for MHC-II
#repitope_full_nonIEDB_human_MHCII <- filter(repitope_full_nonIEDB_human, MHCType == "MHC-II") #30223
#dim(repitope_full_nonIEDB_human_MHCII)
#head(repitope_full_nonIEDB_human_MHCII)
#write.table(repitope_full_nonIEDB_human_MHCII,"repitope_IEDB_human_mhcII_30223.csv",quote=F,sep=",")

Dataset_Filters <- c("Repitope","non-IEDB","Homo sapiens (host)", "MHC-I", "MHC-II", "Negative", "Positive")
Count <- c(70543,8104,NA,8104,0,5505, 2599)
non_iedb_table <- as.data.frame(cbind(Dataset_Filters,Count))
png("non_iedb_table.pdf")
non_iedb_table_format <- formattable(non_iedb_table)
non_iedb_table_format
dev.off()

# Deconvolution performed at excel and python. Check descrioption at:
# https://benchling.com/rfarriol-duran/f/lib_RMu18YCf-immunogenicity-predictor/etr_vsbGTYw6-overfitting_control/edit

## Import deconvoluted object to merge with non-iedb dataset

## repitope_full_nonIEDB_MHCI

non_iedb_deconvoluted_data <- read.table("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control/non-iedb_deconvolution_out_nodups.csv", sep = ",", dec = ".", header = T)

non_iedb_deconvoluted <- left_join(x = non_iedb_deconvoluted_data, y = repitope_full_nonIEDB_MHCI, by = "Peptide", keep = T)
