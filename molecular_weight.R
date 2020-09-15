## Bulkiness or Molecular Weight

## Description: Use the "mw" function of peptides to generate
## two new features:
	## mw_epitope -> Epitope Molecular Weight
	## mw_tcr_contact -> TCR Contact Region Molecular Weight	

## Peptides
library(Peptides)

## Upload data
dataset <- read.table("path/to/dataset.csv", sep = ";", header = T)
colnames(peptides)

# Calculate the molecular weight of the Epitope
mw_peptide <- mw(dataset$peptide)


# Calculate the molecular weight of TCR contact region (surrogate of bulkyness)
mw_tcr_contact <- mw(dataset$tcr_contact)

# Append the new features to the original dataset

dataset_mw <- cbind(dataset, mw_peptide, mw_tcr_contact)
