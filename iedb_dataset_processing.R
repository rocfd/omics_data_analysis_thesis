## IEDB dataset obtained from Repitope package
## Data processing

## TCR Contact Region Extraction
setwd("")

## libraries
library(plyr)
library(dplyr)
library(gridExtra)
library(formattable)

## Upload Repitope dataset
filename <- "path/to/repitope_full_70543.csv"
repitope_full <- read.table(filename, sep = ",", header = T, skip = 0, fill = T)
colnames(repitope_full)
head(repitope_full)
dim(repitope_full) #70543

#### Filter Repitope dataset for IEDB database, human as host and MHC-I epitopes.
## Filter source column for contains IEDB
repitope_full_IEDB <- filter(repitope_full, grepl('IEDB', Source)) #62439
dim(repitope_full_IEDB)
write.table(repitope_full_IEDB,"repitope_IEDB_62439.csv",quote=F,sep=",")

## Filter host column for contains homo sapiens
repitope_full_IEDB_human <- filter(repitope_full_IEDB, grepl('Homo sapiens', Host_Organism)) #41318
dim(repitope_full_IEDB_human)
write.table(repitope_full_IEDB,"repitope_IEDB_human_41318.csv",quote=F,sep=",")

## Filter MHCType column for MHC-I
repitope_full_IEDB_human_MHCI <- filter(repitope_full_IEDB_human, MHCType == "MHC-I") #11095
dim(repitope_full_IEDB_human_MHCI)
head(repitope_full_IEDB_human_MHCI)
write.table(repitope_full_IEDB_human_MHCI,"repitope_IEDB_human_mhcI_11095.csv",quote=F,sep=",")


## Filter MHCtype column for MHC-II
repitope_full_IEDB_human_MHCII <- filter(repitope_full_IEDB_human, MHCType == "MHC-II") #30223
dim(repitope_full_IEDB_human_MHCII)
head(repitope_full_IEDB_human_MHCII)
write.table(repitope_full_IEDB_human_MHCII,"repitope_IEDB_human_mhcII_30223.csv",quote=F,sep=",")

Dataset_Filters <- c("Repitope","IEDB","Homo sapiens (host)", "MHC-I", "MHC-II")
Count <- c(70543,62439,41318,11095,30223)
iedb_table <- as.data.frame(cbind(Dataset_Filters,Count))
pdf("iedb_table.pdf")
iedb_table_format <- formattable(iedb_table)
dev.off()

