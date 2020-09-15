##NON-IEDB dataset obtained from Repitope package
## Data filtering, feature extraction and statistical association to Immunogenicity

## TCR Contact Region Extraction
setwd("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control")

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
#filename <- ""
#iedb_original <- read.table(filename, sep = ",", header = T, skip = 0, fill = T)
# Unrequired, object already available.

## non_iedb_deconvoluted
colnames(non_iedb_deconvoluted)
glimpse(non_iedb_deconvoluted)
dim(non_iedb_deconvoluted) #8797



### Extract peptides
## Extract peptide
peptide <- non_iedb_deconvoluted$Peptide.x
peptide_vector <- as.vector(peptide)
length(peptide_vector)


## Factorize Immunogenicity and reorder to have positive first 
non_iedb_deconvoluted$Immunogenicity <- factor(non_iedb_deconvoluted$Immunogenicity, levels = c("Negative","Positive"))
print(non_iedb_deconvoluted$Immunogenicity)


### Extract TCR contating residues from every peptide. P4-P6
## Extract residues 4-6 from peptides
tcr_contact <- substring(peptide_vector,4,6)
tcr_contact
non_iedb_deconvoluted <- cbind(non_iedb_deconvoluted,tcr_contact)

## Peptide and TCR contacting region FEATURES

## Bulkiness or Molecular Weight
# Calculate the molecular weight at peptide level
mw_peptide <- mw(peptide_vector)
head(mw_peptide)
# append
non_iedb_deconvoluted <- cbind(non_iedb_deconvoluted,mw_peptide)
# save figure
pdf(file = "figures_overfitting/hist_non_iedbb_mw_peptide.pdf")
hist(mw_peptide, main = "Molecular Weight - Peptide level \n IEDB independent dataset", xlab = "Molecular Weight")
dev.off()


# Calculate the molecular weight at TCR contact level
mw_tcr_contact <- mw(tcr_contact)
head(mw_tcr_contact)

# append
non_iedb_deconvoluted  <- cbind(non_iedb_deconvoluted,mw_tcr_contact)

pdf(file = "figures_overfitting/hist_non_iedbb_mw_tcr_contact.pdf")
hist(mw_tcr_contact, main = "Molecular Weight - TCR contact level \n IEDB independent dataset", xlab = "Molecular Weight")
dev.off()

## Hydrophobicity
# Calculate the hydrophbicity at peptide level
hydroph_peptide <- hydrophobicity(peptide_vector)
head(hydroph_peptide)

# append
non_iedb_deconvoluted <- cbind(non_iedb_deconvoluted,hydroph_peptide)

pdf(file = "figures_overfitting/hist_non_iedbb_hydrophobicity_peptide.pdf")
hist(hydroph_peptide, main = "Hydrophobicity - Peptide level \n IEDB independent dataset", xlab = "GRAVY index")
dev.off()

# Calculate the hydrophbicity at TCR contact level
hydroph_tcr_contact <- hydrophobicity(tcr_contact)
head(hydroph_peptide)

# append
non_iedb_deconvoluted <- cbind(non_iedb_deconvoluted,hydroph_tcr_contact)

pdf(file = "figures_overfitting/hist_non_iedbb_hydrophobicity_tcr_contact.pdf")
hist(hydroph_tcr_contact, main = "Hydrophobicity - TCR contact level \n IEDB independent dataset", xlab = "GRAVY index")
dev.off()

## Netcharge
# Calculate the net charge at peptide level
charge_peptide <- charge(peptide_vector)
head(charge_peptide)

# append
non_iedb_deconvoluted <- cbind(non_iedb_deconvoluted,charge_peptide)

pdf(file = "figures_overfitting/hist_non_iedbb_charge_peptide.pdf")
hist(charge_peptide, main = "Net charge - Peptide level \n IEDB independent dataset", xlab = "Charge")
dev.off()

# Calculate the net charge at TCR contacting region
charge_tcr_contact <- charge(tcr_contact)
head(charge_tcr_contact)

# append
non_iedb_deconvoluted <- cbind(non_iedb_deconvoluted,charge_tcr_contact)

pdf(file = "figures_overfitting/hist_non_iedbb_charge_tcr_contact.pdf")
hist(charge_tcr_contact, main = "Net charge - TCR contact level \n IEDB independent dataset", xlab = "Charge")
dev.off()


## Peptide instability index
# Calculate stability at peptide level
stab_peptide <- instaIndex(peptide_vector)

# append
non_iedb_deconvoluted <- cbind(non_iedb_deconvoluted,stab_peptide)

pdf(file = "figures_overfitting/hist_non_iedbb_stability_peptide.pdf")
hist(stab_peptide, main = "Stability - Peptide level \n IEDB independent dataset", xlab = "Peptide instability index")
dev.off()


# Calculate the net charge at TCR contact level
stab_tcr_contact <- instaIndex(tcr_contact)

# append
non_iedb_deconvoluted <- cbind(non_iedb_deconvoluted,stab_tcr_contact)

pdf(file = "figures_overfitting/hist_non_iedbb_stability_tcr_contact.pdf")
hist(stab_tcr_contact, main = "Stability - TCR contact level \n IEDB independent dataset", xlab = "Peptide instability index")
dev.off()


### ADD peptide features to dataset
## Append all the properties to iedb_human_MHCI dataset
#iedb_properties <- cbind(iedb,mw_peptide,hydroph_peptide,charge_peptide,stab_peptide,mw_tcr_contact,hydroph_tcr_contact,charge_tcr_contact,stab_tcr_contact)
#head(iedb)

## Filter for positive properties
non_iedbb_human_MHCI_pos <- filter(non_iedb_deconvoluted, Immunogenicity == "Positive")
non_iedbb_human_MHCI_neg <- filter(non_iedb_deconvoluted, Immunogenicity == "Negative")

## Peptide properties Vs. Immunogenicity
# Color palette for plots
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Peptide Hydrophobicity Vs. Immunogenicity
#Normality testing
hist(non_iedbb_human_MHCI_pos$hydroph_peptide)
hist(non_iedbb_human_MHCI_neg$hydroph_peptide)
ks.test(non_iedbb_human_MHCI_pos$hydroph_peptide,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(non_iedbb_human_MHCI_neg$hydroph_peptide,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
hydroph_pept_plot <- ggplot(non_iedb_deconvoluted, aes(Immunogenicity,hydroph_peptide)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 5.5, method.args = list(alternative = "two.sided"))
hydroph_pept_plot

## Peptide Bulkiness Vs. Immunogenicity
#Normality testing
hist(non_iedbb_human_MHCI_pos$mw_peptide)
hist(non_iedbb_human_MHCI_neg$mw_peptide)
ks.test(non_iedbb_human_MHCI_pos$mw_peptide,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(non_iedbb_human_MHCI_neg$mw_peptide,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
mw_pept_plot <- ggplot(non_iedb_deconvoluted, aes(Immunogenicity,mw_peptide)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 1700, method.args = list(alternative = "two.sided"))
mw_pept_plot

## Peptide Charge Vs. Immunogenicity
#Normality testing
hist(non_iedbb_human_MHCI_pos$charge_peptide)
hist(non_iedbb_human_MHCI_neg$charge_peptide)
ks.test(non_iedbb_human_MHCI_pos$charge_peptide,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(non_iedbb_human_MHCI_neg$charge_peptide,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
charge_pept_plot <- ggplot(non_iedb_deconvoluted, aes(Immunogenicity,charge_peptide)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 8, method.args = list(alternative = "two.sided"))
charge_pept_plot
median(non_iedbb_human_MHCI_pos$charge_peptide)
median(non_iedbb_human_MHCI_neg$charge_peptide)

## Peptide Stability Vs. Immunogenicity
#Normality testing
hist(non_iedbb_human_MHCI_pos$stab_peptide)
hist(non_iedbb_human_MHCI_neg$stab_peptide)
ks.test(non_iedbb_human_MHCI_pos$stab_peptide,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(non_iedbb_human_MHCI_neg$stab_peptide,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
stab_pept_plot <- ggplot(non_iedb_deconvoluted, aes(Immunogenicity,stab_peptide)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative","Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 350, method.args = list(alternative = "two.sided"))
stab_pept_plot


## TCR Bulkiness Vs. Immunogenicity
#Normality testing
hist(non_iedbb_human_MHCI_pos$mw_tcr_contact)
hist(non_iedbb_human_MHCI_neg$mw_tcr_contact)
ks.test(non_iedbb_human_MHCI_pos$mw_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(non_iedbb_human_MHCI_neg$mw_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
mw_tcr_plot <- ggplot(non_iedb_deconvoluted, aes(Immunogenicity,mw_tcr_contact)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 600 ,method.args = list(alternative = "two.sided"))
mw_tcr_plot

## TCR contact hydrophobicity & Immunogenicity
#Normality testing
hist(non_iedbb_human_MHCI_pos$hydroph_tcr_contact)
hist(non_iedbb_human_MHCI_neg$hydroph_tcr_contact)
ks.test(non_iedbb_human_MHCI_pos$hydroph_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(non_iedbb_human_MHCI_neg$hydroph_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test (development)
hydroph_tcr_plot <- ggplot(non_iedb_deconvoluted, aes(Immunogenicity,hydroph_tcr_contact)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 6, method.args = list(alternative = "two.sided"))
hydroph_tcr_plot

## TCR Charge Vs. Immunogenicity
#Normality testing
hist(non_iedbb_human_MHCI_pos$charge_tcr_contact)
hist(non_iedbb_human_MHCI_neg$charge_tcr_contact)
ks.test(non_iedbb_human_MHCI_pos$charge_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(non_iedbb_human_MHCI_neg$charge_tcr_contact,"pnorm",alternative="two.sided") #Not normally distributed

## Boxplot and Wilcoxon.test
charge_tcr_plot <- ggplot(non_iedb_deconvoluted, aes(Immunogenicity,charge_tcr_contact)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 4, method.args = list(alternative = "greater"))
charge_tcr_plot
########################3 PQs fo far

## Export dataset to search using NETCTLpan 1.1
write.table(x = non_iedb_deconvoluted, file = "non_iedb_deconvoluted_toNetCTL.csv", quote = F, sep = ";", dec = ".")

## Export dataset to search using MHCflurry
write.table(x= non_iedb_deconvoluted, file = "non_iedb_deconvoluted_toMHCflurry.csv", quote = F, sep = ",", dec = ".")

## Import mhcflurry predicted non-iedb dataset
non_iedb_mhcflurry <- read.table(file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control/mhcflurry_non_iedb/non_iedb_deconvoluted_toMHCflurry_out_tomerge.csv", sep = ",", header = T, fill = T, dec = ".")
dim(non_iedb_mhcflurry) #8797
dim(non_iedb_deconvoluted) #8797

## Merge MHCflurry imported with non_iedb_deconvoluted
non_iedb_deconvoluted_mhcflurry <- cbind(non_iedb_deconvoluted, non_iedb_mhcflurry)
View(non_iedb_deconvoluted_mhcflurry)
dim(non_iedb_deconvoluted_mhcflurry)
glimpse(non_iedb_deconvoluted_mhcflurry)
class(non_iedb_deconvoluted$Immunogenicity)
print(non_iedb_deconvoluted_mhcflurry$Immunogenicity)

## Filter for positive properties
#non_iedbb_flurry_pos <- filter(non_iedb_deconvoluted_mhcflurry, Immunogenicity == "Positive")
#non_iedbb_flurry_neg <- filter(non_iedb_deconvoluted_mhcflurry, Immunogenicity == "Negative")

# Export and spilt in excel fuck
write.table(x= non_iedb_deconvoluted_mhcflurry, file = "non_iedb_deconvoluted_MHCflurry_tosplit.csv", quote = F, sep = ",", dec = ".")

## Import splitted
non_iedbb_flurry_pos <- read.table(file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control/non_iedb_deconvoluted_MHCflurry_positive.csv", sep = ",", header = T, fill = T, dec = ".")
non_iedbb_flurry_neg <- read.table(file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control/non_iedb_deconvoluted_MHCflurry_negative.csv", sep = ",", header = T, fill = T, dec = ".")

## Colnames
colnames(non_iedb_deconvoluted_mhcflurry)
##"mhcflurry_affinity"            "mhcflurry_best_allele"        
## [28] "mhcflurry_affinity_percentile" "mhcflurry_processing_score"    "mhcflurry_presentation_score" 

## Rename peptide x
colnames(non_iedb_deconvoluted_mhcflurry)[24] <- "peptide3"
colnames(non_iedb_deconvoluted_mhcflurry)[25] <- "deco_allele"

## Import formatted mhcflurry deconvoluted to merge with NETCTL pred
non_iedb_deconvoluted_mhcflurry_format <- read.table(file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control/non_iedb_deconvoluted_MHCflurry_toMergeNetCTL.csv", sep = ",", header = T, fill = T, dec = ".")
dim(non_iedb_deconvoluted_mhcflurry_format) #8797
View(non_iedb_deconvoluted_mhcflurry_format)

non_iedb_deconvoluted_netctl <- read.table(file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control/netctl_non_iedb/processing/noniedb_netctlpan_proc.csv", sep = ",", header = T, fill = T, dec = ".")
dim(non_iedb_deconvoluted_netctl) # 8797
View(non_iedb_deconvoluted_netctl)

## Join netctl and netmfclurry
non_iedb_deconvoluted_all <- left_join(x = non_iedb_deconvoluted_mhcflurry_format, y = non_iedb_deconvoluted_netctl, by = "Sequence.Name", keep = T)
dim(non_iedb_deconvoluted_all)
View(non_iedb_deconvoluted_all)

## Fixed merging problem with netCTL and MHCflurry data.

## lead to modelling

## MHCflurry binding
## Boxplot and Wilcoxon.test
binding_plot <- ggplot(non_iedb_deconvoluted_mhcflurry, aes(Immunogenicity,mhcflurry_affinity)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 5, method.args = list(alternative = "two.sided"))
binding_plot

## MHCflurry percentile
## Boxplot and Wilcoxon.test
percentile_plot <- ggplot(non_iedb_deconvoluted_mhcflurry, aes(Immunogenicity,mhcflurry_affinity_percentile)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 5, method.args = list(alternative = "two.sided"))
percentile_plot

## MHCflurry processing
## Boxplot and Wilcoxon.test
processing_plot <- ggplot(non_iedb_deconvoluted_mhcflurry, aes(Immunogenicity,mhcflurry_processing_score)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 1.2, method.args = list(alternative = "two.sided"))
processing_plot

## MHCflurry presentation
## Boxplot and Wilcoxon.test
presentation_plot <- ggplot(non_iedb_deconvoluted_mhcflurry, aes(Immunogenicity,mhcflurry_presentation_score)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 1.2, method.args = list(alternative = "two.sided"))
presentation_plot




## TAP score Vs. Immunogenicity
#Normality testing
hist(IEDB_human_MHCI_pos$TAP_efficiency)
hist(IEDB_human_MHCI_neg$TAP_efficiency)
ks.test(IEDB_human_MHCI_pos$TAP_efficiency,"pnorm",alternative="two.sided") #Not normally distributed
ks.test(IEDB_human_MHCI_neg$TAP_efficiency,"pnorm",alternative="two.sided") #Not normally distributed

pdf(file = "figures_overfitting/hist_iedb_TAP_peptide.pdf")
hist(iedb_TAP$TAP_efficiency, main = "TAP transport efficiency - Peptide level", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures_overfitting/hist_iedb_TAP_peptide_pos.pdf")
hist(IEDB_human_MHCI_pos$TAP_efficiency, main = "TAP transport efficiency - Immuno positive", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures_overfitting/hist_iedb_TAP_peptide_neg.pdf")
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

pdf(file = "figures_overfitting/hist_iedb_Cle_peptide.pdf")
hist(iedb_TAP$Proteasomal_Cleavage, main = "TAP transport efficiency - Peptide level", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures_overfitting/hist_iedb_Cle_peptide_pos.pdf")
hist(IEDB_human_MHCI_pos$Proteasomal_Cleavage, main = "TAP transport efficiency - Immuno positive", xlab = "TAP efficiency")
dev.off()

pdf(file = "figures_overfitting/hist_iedb_cle_peptide_neg.pdf")
hist(IEDB_human_MHCI_neg$Proteasomal_Cleavage, main = "TAP transport efficiency - Immuno negative", xlab = "TAP efficiency")
dev.off()

TAP_Cle_plot <- ggplot(iedb_TAP, aes(Immunogenicity,Proteasomal_Cleavage)) +
  geom_boxplot(color=1,fill=c("#999999", "#E69F00")) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE) + 
  stat_compare_means(method= "wilcox.test", label.y = 1.25, method.args = list(alternative = "greater"))
TAP_Cle_plot



## Multiplot
hydroph_pept_plot
mw_pept_plot
charge_pept_plot
stab_pept_plot
hydroph_tcr_plot
mw_tcr_plot
charge_tcr_plot

pdf("figures_overfitting/IEDB_peptide_properties_immunogenicity.pdf")
peptide_properties_plot<- grid.arrange(hydroph_pept_plot,mw_pept_plot,charge_pept_plot,stab_pept_plot, ncol=2, top = "Peptide properties Vs. Immunogenicity")
dev.off()

pdf("figures_overfitting/IEDB_tcrcontact_properties_immunogenicity.pdf")
tcr_properties_plot <- grid.arrange(hydroph_tcr_plot,mw_tcr_plot,charge_tcr_plot, ncol=3, top = "TCR contacting residues (P4-P6) properties Vs. Immunogenicity")
dev.off()

pdf("figures_overfitting/IEDB_peptide_processing_immunogenicity.pdf")
peptide_processing_plot <- grid.arrange(TAP_plot, TAP_Cle_plot, TAP_MHC_plot, TAP_comb_plot)
dev.off()
