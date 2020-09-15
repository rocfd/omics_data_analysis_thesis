# EPITOPE STABILITY

## Description: Use the "charge" function of peptides to generate
## two new features:
	## stab_epitope -> Epitope Stability
	## stab_tcr_contact -> TCR Contact Region Stability


## Peptides
library(Peptides)

## Upload data
peptides <- read.table("path/to/dataset.csv", sep = ";", header = T)

## Stability

# Calculate the net charge of the epitope
stab_epitope <- instaIndex(dataset$peptide)

# Calculate the net charge of the  TCR contact region (surrogate of bulkyness)
stab_tcr_contact <- instaIndex(dataset$tcr_contact)

# Append the new features to the original dataset
dataset_stab <- cbind(dataset, stab_epitope, stab_tcr_contact)


