## Net Charge Status

## Description: Use the "charge" function of peptides to generate
## two new features:
	## charge_epitope -> Epitope Molecular Weight
	## charge_tcr_contact -> TCR Contact Region Molecular Weight


## Peptides
library(Peptides)

## Upload data
dataset <- read.table("path/to/dataset.csv", sep = ";", header = T)
#colnames(dataset)

## netcharge

# Calculate the net charge of the epitope
charge_epitope <- charge(dataset$peptide)


# Calculate the net charge of the TCR contact region (surrogate of bulkyness)
charge_tcr_contact <- charge(dataset$tcr_contact)

# Append the new features to the original dataset
dataset_charge <- cbind(dataset, charge_epitope, charge_tcr_contact)

