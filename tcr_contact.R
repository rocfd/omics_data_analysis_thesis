## TCR Contact Region Extraction

## Description: extract the TCR contact region (P4-P6) from any given epitope vector with a dataframe.

## Use the "substring" function of Base to generate
## two new features:
	## tcr_contact -> TCR Contact region (P4, P5, P6)


## Peptides
library(Peptides)

## Upload data

dataset <- read.table("/path/to/dataset.csv", sep = ";", header = T)

## Extract epitopes form dataframe to vector
epitope_vector <- as.vector(dataset$peptide)

## Extract residues 4-6 from mut peptides
tcr_contact <- substring(epitope_vector,4,6)


#Include contact region in dataframe
dataset_tcr_contact <- cbind(peptides,tcr_contact)
