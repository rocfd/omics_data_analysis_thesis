## IEDB dataset obtained from Repitope package
## Data filtering, feature extraction and statistical association to Immunogenicity

## TCR Contact Region Extraction
setwd("/Users/rocfarriolduran/Desktop/BSC/AA_pipeline/20200723_IEDB_peptide_properties")

## libraries
library(Peptides)
library(plyr)
library(dplyr)
library(sjmisc)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(gridExtra)

## Upload Repitope IEDB Human MHC-I dataset
filename <- "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200829_IEDB_dataset_processing/repitope_IEDB_human_mhcI_11095.csv"
iedb_original <- read.table(filename, sep = ",", header = T, skip = 0, fill = T)
colnames(iedb)
head(iedb)
dim(iedb) #11095


### Extract peptides
## Extract peptide
peptide <- iedb$Peptide
peptide_vector <- as.vector(peptide)
length(peptide_vector)

## Factorize Immunogenicity and reorder to have positive first 
iedb$Immunogenicity <- factor(iedb$Immunogenicity, levels = c("Negative","Positive"))
print(iedb$Immunogenicity)


### Extract TCR contating residues from every peptide. P4-P6
## Extract residues 4-6 from peptides
tcr_contact <- substring(peptide_vector,4,6)
tcr_contact
iedb <- cbind(iedb,tcr_contact)

## Peptide and TCR contacting region FEATURES

## Bulkiness or Molecular Weight
# Calculate the molecular weight at peptide level
mw_peptide <- mw(peptide_vector)
head(mw_peptide)
# append
iedb <- cbind(iedb,mw_peptide)
# save figure
pdf(file = "figures/hist_iedb_mw_peptide.pdf")
hist(mw_peptide, main = "Molecular Weight - Peptide level", xlab = "Molecular Weight")
dev.off()


# Calculate the molecular weight at TCR contact level
mw_tcr_contact <- mw(tcr_contact)
head(mw_tcr_contact)

# append
iedb <- cbind(iedb,mw_tcr_contact)

pdf(file = "figures/hist_iedb_mw_tcr_contact.pdf")
hist(mw_tcr_contact, main = "Molecular Weight - TCR contact level", xlab = "Molecular Weight")
dev.off()

## Hydrophobicity
# Calculate the hydrophbicity at peptide level
hydroph_peptide <- hydrophobicity(peptide_vector)
head(hydroph_peptide)

# append
iedb <- cbind(iedb,hydroph_peptide)

pdf(file = "figures/hist_iedb_hydrophobicity_peptide.pdf")
hist(hydroph_peptide, main = "Hydrophobicity - Peptide level", xlab = "GRAVY index")
dev.off()

# Calculate the hydrophbicity at TCR contact level
hydroph_tcr_contact <- hydrophobicity(tcr_contact)
head(hydroph_peptide)

# append
iedb <- cbind(iedb,hydroph_tcr_contact)

pdf(file = "figures/hist_iedb_hydrophobicity_tcr_contact.pdf")
hist(hydroph_tcr_contact, main = "Hydrophobicity - TCR contact level", xlab = "GRAVY index")
dev.off()

## Netcharge
# Calculate the net charge at peptide level
charge_peptide <- charge(peptide_vector)
head(charge_peptide)

# append
iedb <- cbind(iedb,charge_peptide)

pdf(file = "figures/hist_iedb_charge_peptide.pdf")
hist(charge_peptide, main = "Net charge - Peptide level", xlab = "Charge")
dev.off()

# Calculate the net charge at TCR contacting region
charge_tcr_contact <- charge(tcr_contact)
head(charge_tcr_contact)

# append
iedb <- cbind(iedb,charge_tcr_contact)

pdf(file = "figures/hist_iedb_charge_tcr_contact.pdf")
hist(charge_tcr_contact, main = "Net charge - TCR contact level", xlab = "Charge")
dev.off()


## Peptide instability index
# Calculate stability at peptide level
stab_peptide <- instaIndex(peptide_vector)

# append
iedb <- cbind(iedb,stab_peptide)

pdf(file = "figures/hist_iedb_stability_peptide.pdf")
hist(stab_peptide, main = "Stability - Peptide level", xlab = "Peptide instability index")
dev.off()


# Calculate the net charge at TCR contact level
stab_tcr_contact <- instaIndex(tcr_contact)

# append
iedb <- cbind(iedb,stab_tcr_contact)

pdf(file = "figures/hist_iedb_stability_tcr_contact.pdf")
hist(stab_tcr_contact, main = "Stability - TCR contact level", xlab = "Peptide instability index")
dev.off()


### ADD peptide features to dataset
## Append all the properties to iedb_human_MHCI dataset
#iedb_properties <- cbind(iedb,mw_peptide,hydroph_peptide,charge_peptide,stab_peptide,mw_tcr_contact,hydroph_tcr_contact,charge_tcr_contact,stab_tcr_contact)
head(iedb)

## Filter for positive properties
IEDB_human_MHCI_pos <- filter(iedb, Immunogenicity == "Positive")
IEDB_human_MHCI_neg <- filter(iedb, Immunogenicity == "Negative")

## Peptide properties Vs. Immunogenicity
# Color palette for plots
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Peptide Hydrophobicity Vs. Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$hydroph_peptide)
hist(IEDB_human_MHCI_neg$hydroph_peptide)
ks.test(IEDB_human_MHCI_pos$hydroph_peptide,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$hydroph_peptide,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
hydroph_pept_plot <- ggplot(iedb, aes(Immunogenicity,hydroph_peptide)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 5.5, method.args = list(alternative = "less"))
hydroph_pept_plot

## Peptide Bulkiness Vs. Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$mw_peptide)
hist(IEDB_human_MHCI_neg$mw_peptide)
ks.test(IEDB_human_MHCI_pos$mw_peptide,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$mw_peptide,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
mw_pept_plot <- ggplot(iedb, aes(Immunogenicity,mw_peptide)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 1600, method.args = list(alternative = "greater"))
mw_pept_plot

## Peptide Charge Vs. Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$charge_peptide)
hist(IEDB_human_MHCI_neg$charge_peptide)
ks.test(IEDB_human_MHCI_pos$charge_peptide,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$charge_peptide,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
charge_pept_plot <- ggplot(iedb, aes(Immunogenicity,charge_peptide)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 8, method.args = list(alternative = "two.sided"))
charge_pept_plot
median(IEDB_human_MHCI_pos$charge_peptide)
median(IEDB_human_MHCI_neg$charge_peptide)

## Peptide Stability Vs. Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$stab_peptide)
hist(IEDB_human_MHCI_neg$stab_peptide)
ks.test(IEDB_human_MHCI_pos$stab_peptide,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$stab_peptide,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
stab_pept_plot <- ggplot(iedb, aes(Immunogenicity,stab_peptide)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative","Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 350, method.args = list(alternative = "less"))
stab_pept_plot


## TCR Bulkiness Vs. Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$mw_tcr_contact)
hist(IEDB_human_MHCI_neg$mw_tcr_contact)
ks.test(IEDB_human_MHCI_pos$mw_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$mw_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
mw_tcr_plot <- ggplot(iedb, aes(Immunogenicity,mw_tcr_contact)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 600 ,method.args = list(alternative = "greater"))
mw_tcr_plot

## TCR contact hydrophobicity & Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$hydroph_tcr_contact)
hist(IEDB_human_MHCI_neg$hydroph_tcr_contact)
ks.test(IEDB_human_MHCI_pos$hydroph_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$hydroph_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test (development)
hydroph_tcr_plot <- ggplot(iedb, aes(Immunogenicity,hydroph_tcr_contact)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 6, method.args = list(alternative = "less"))
hydroph_tcr_plot

## TCR Charge Vs. Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$charge_tcr_contact)
hist(IEDB_human_MHCI_neg$charge_tcr_contact)
ks.test(IEDB_human_MHCI_pos$charge_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$charge_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
charge_tcr_plot <- ggplot(iedb, aes(Immunogenicity,charge_tcr_contact)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 4, method.args = list(alternative = "greater"))
charge_tcr_plot

## Filter for positive properties
IEDB_human_MHCI_pos <- filter(iedb_TAP, Immunogenicity == "Positive")
IEDB_human_MHCI_neg <- filter(iedb_TAP, Immunogenicity == "Negative")

## TAP score Vs. Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$TAP_efficiency)
hist(IEDB_human_MHCI_neg$TAP_efficiency)
ks.test(IEDB_human_MHCI_pos$TAP_efficiency,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$TAP_efficiency,"pnorm",alternative="two.sided") #Not normally distributed

pdf(file = "figures/hist_iedb_TAP_peptide.pdf")
hist(iedb_TAP$TAP_efficiency, main = "TAP transport efficiency - Peptide level", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures/hist_iedb_TAP_peptide_pos.pdf")
hist(IEDB_human_MHCI_pos$TAP_efficiency, main = "TAP transport efficiency - Immuno positive", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures/hist_iedb_TAP_peptide_neg.pdf")
hist(IEDB_human_MHCI_neg$TAP_efficiency, main = "TAP transport efficiency - Immuno negative", xlab = "TAP efficiency")
dev.off()

## Boxplot and Wilcoxon.test
TAP_plot <- ggplot(iedb_TAP, aes(Immunogenicity,TAP_efficiency)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 5, method.args = list(alternative = "greater"))
TAP_plot

# Proteasomal Cleavage Vs. Immunogenicity
hist(IEDB_human_MHCI_pos$Proteasomal_Cleavage)
hist(IEDB_human_MHCI_neg$Proteasomal_Cleavage)
ks.test(IEDB_human_MHCI_pos$TAP,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$TAP,"pnorm",alternative="two.sided") #Not normally distributed

pdf(file = "figures/hist_iedb_Cle_peptide.pdf")
hist(iedb_TAP$Proteasomal_Cleavage, main = "TAP transport efficiency - Peptide level", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures/hist_iedb_Cle_peptide_pos.pdf")
hist(IEDB_human_MHCI_pos$Proteasomal_Cleavage, main = "TAP transport efficiency - Immuno positive", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures/hist_iedb_cle_peptide_neg.pdf")
hist(IEDB_human_MHCI_neg$Proteasomal_Cleavage, main = "TAP transport efficiency - Immuno negative", xlab = "TAP efficiency")
dev.off()

TAP_Cle_plot <- ggplot(iedb_TAP, aes(Immunogenicity,Proteasomal_Cleavage)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 1.25, method.args = list(alternative = "greater"))
TAP_Cle_plot

#MHC_score Vs. Immunogenicity
hist(IEDB_human_MHCI_pos$MHC_score)
hist(IEDB_human_MHCI_neg$MHC_score)
ks.test(IEDB_human_MHCI_pos$TAP,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$TAP,"pnorm",alternative="two.sided") #Not normally distributed

pdf(file = "figures/hist_iedb_MHC_peptide.pdf")
hist(iedb_TAP$MHC_score, main = "TAP transport efficiency - Peptide level", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures/hist_iedb_MHC_peptide_pos.pdf")
hist(IEDB_human_MHCI_pos$MHC_score, main = "TAP transport efficiency - Immuno positive", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures/hist_iedb_MHC_peptide_neg.pdf")
hist(IEDB_human_MHCI_neg$MHC_score, main = "TAP transport efficiency - Immuno negative", xlab = "TAP efficiency")
dev.off()

TAP_MHC_plot <- ggplot(iedb_TAP, aes(Immunogenicity,MHC_score)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 1.15, method.args = list(alternative = "greater"))
TAP_MHC_plot

## Combined Score Vs. Immunogenicity
hist(IEDB_human_MHCI_pos$NetCTL_Combined_Score)
hist(IEDB_human_MHCI_neg$NetCTL_Combined_Score)
ks.test(IEDB_human_MHCI_pos$TAP,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$TAP,"pnorm",alternative="two.sided") #Not normally distributed

pdf(file = "figures/hist_iedb_comb_peptide.pdf")
hist(iedb_TAP$NetCTL_Combined_Score, main = "TAP transport efficiency - Peptide level", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures/hist_iedb_comb_peptide_pos.pdf")
hist(IEDB_human_MHCI_pos$NetCTL_Combined_Score, main = "TAP transport efficiency - Immuno positive", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures/hist_iedb_comb_peptide_neg.pdf")
hist(IEDB_human_MHCI_neg$NetCTL_Combined_Score, main = "TAP transport efficiency - Immuno negative", xlab = "TAP efficiency")
dev.off()

TAP_comb_plot <- ggplot(iedb_TAP, aes(Immunogenicity,NetCTL_Combined_Score)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 1.5, method.args = list(alternative = "greater"))
TAP_comb_plot

## Multiplot
hydroph_pept_plot
mw_pept_plot
charge_pept_plot
stab_pept_plot
hydroph_tcr_plot
mw_tcr_plot
charge_tcr_plot

pdf("figures/IEDB_peptide_properties_immunogenicity.pdf")
peptide_properties_plot<- grid.arrange(hydroph_pept_plot,mw_pept_plot,charge_pept_plot,stab_pept_plot, ncol=2, top = "Peptide properties Vs. Immunogenicity")
dev.off()

pdf("figures/IEDB_tcrcontact_properties_immunogenicity.pdf")
tcr_properties_plot <- grid.arrange(hydroph_tcr_plot,mw_tcr_plot,charge_tcr_plot, ncol=3, top = "TCR contacting residues (P4-P6) properties Vs. Immunogenicity")
dev.off()

pdf("figures/IEDB_peptide_processing_immunogenicity.pdf")
peptide_processing_plot <- grid.arrange(TAP_plot, TAP_Cle_plot, TAP_MHC_plot, TAP_comb_plot)
dev.off()
