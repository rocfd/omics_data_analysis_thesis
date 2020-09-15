## Hydrophobicity


## Description: Use the "charge" function of peptides to generate
## two new features:
	## hydroph_epitope -> Epitope Hydrophobicity
	## hydroph_tcr_contact -> TCR Contact Region Hydrophobicity


## Peptides
library(Peptides)

## Upload data
peptides <- read.table("path/to/dataset.csv", sep = ";", header = T)


## Hydrophobicity

# Calculate the hydrophbicity of the epitope
hydroph_epitope <- hydrophobicity(dataset$peptide)


# Calculate the hydrophbicity of the  TCR contact region (surrogate of bulkyness)
hydroph_tcr_contact <- hydrophobicity(dataset$tcr_contact)

# Append the new features to the original dataset
dataset_hydroph <- cbind(dataset, hydroph_epitope, hydroph_tcr_contact)

