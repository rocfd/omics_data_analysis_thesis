## FACET Allele Peptide Properties

## Goal: upload IEDB MHCrestriction faceted data "top_allele_facet" and join with iedb_TAP_mhcflurryrestricted data.
## Plot the peptide properties faceted per allele.

## DEVELOPMENT

setwd("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200831_IEDB_peptide_properties_per_Allele/facet_allele_properties/figures_perallele")

## libraries
library(dplyr)
library(ggplot2)
library(ggforce)
library(Peptides)
library(plyr)
library(sjmisc)
library(ggsignif)
library(ggpubr)
library(gridExtra)
library(gg.gap)

## Import IEDB faceted data:
iedb_faceted <- read.table("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200831_IEDB_peptide_properties_per_Allele/facet_allele_properties/iedb_alleles_table_dictionary.csv", header = T, sep = ",", dec = ".", stringsAsFactors = F)
dim(iedb_faceted) #25875
colnames(iedb_faceted)[which(names(iedb_faceted) == "Peptide")] <- "peptide"

## Filter HLA-C*0601
iedb_faceted <- filter(iedb_faceted,top_allele_facet !="HLA-C*0602")
dim(iedb_faceted) ## 25837

## Refer to iedb_tap_mhcflurryrestricted --> data to join
# iedb_TAP_mhcflurryrestricted
# rename peptide.x for peptide
colnames(iedb_TAP_mhcflurryrestricted)[which(names(iedb_TAP_mhcflurryrestricted) == "peptide.x")] <- "peptide"


## Join data by peptide column
#iedb_faceted_properties <- inner_join(iedb_faceted,iedb_TAP_mhcflurryrestricted, by = "peptide", keep = T)
#dim(iedb_faceted_properties) # 25965 
## excess due to duplicates in peptide

## remove peptide duplicates at iedb_TAP_mhcflurryrestricted
iedb_TAP_mhcflurryrestricted %>% distinct(peptide, .keep_all = TRUE) -> iedb_TAP_mhcflurryrestricted_nodup
colnames(iedb_TAP_mhcflurryrestricted_nodup)[which(names(iedb_TAP_mhcflurryrestricted_nodup) == "peptide.x")] <- "peptide"
dim(iedb_TAP_mhcflurryrestricted_nodup)

## remove dummy HLA-A0201
iedb_TAP_mhcflurryrestricted_nodup <- subset(x = iedb_TAP_mhcflurryrestricted_nodup, select = -c(allele))

## readd allele column by merging with contains mhc_restriction_genotype_NA
iedb_TAP_mhcflurryrestricted_nodup %>% mutate(allele = mhc_restriction_genotype_NA) -> iedb_TAP_mhcflurryrestricted_nodup
colnames(iedb_TAP_mhcflurryrestricted_nodup)

# remove biased A0201 prediction
iedb_TAP_mhcflurryrestricted_nodup <- subset(x = iedb_TAP_mhcflurryrestricted_nodup, select = -c(mhcflurry_affinity, mhcflurry_best_allele,mhcflurry_affinity_percentile, mhcflurry_processing_score, mhcflurry_presentation_score))

## export to recalculate using mhcflurry
write.table(x = iedb_TAP_mhcflurryrestricted_nodup, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200831_IEDB_peptide_properties_per_Allele/mhcflurry_iedb_mhcrestriction_bestallele_corrected/iedb_peptide_mhcrestriction_parsed_tomhcflurry_bestAllele.csv", sep = ",", dec = ".", quote = F)

## recalculate using mhcflurry
## mhcflurry_iedb_mhcrestriction_bestallele_corrected rocfarriolduran$ mhcflurry-predict iedb_peptide_mhcrestriction_parsed_tomhcflurry_bestAllele.csv --out iedb_peptide_mhcrestriction_parsed_tomhcflurry_bestAllele_out.csv --no-throw --output-delimiter ',' --prediction-column-prefix 'bestAllele_mhcflurry_' --always-include-best-allele

## import recalcualted data
iedb_TAP_mhcflurryrestricted_nodup_recalc <- read.table("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200831_IEDB_peptide_properties_per_Allele/mhcflurry_iedb_mhcrestriction_bestallele_corrected/iedb_peptide_mhcrestriction_parsed_tomhcflurry_bestAllele_out.csv", header = T, sep = ",", dec = ".", stringsAsFactors = T)
dim(iedb_TAP_mhcflurryrestricted_nodup_recalc)

## Join data by peptide without duplicates
iedb_faceted_properties_nodup <- left_join(iedb_faceted,iedb_TAP_mhcflurryrestricted_nodup_recalc, by = "peptide", keep = T)
dim(iedb_faceted_properties_nodup) # 25837

## Ex
write.table(iedb_faceted_properties_nodup, file = "iedb_faceted_properties_nodup.csv", sep = ";", quote = F, dec = ".")

### PLOTING
## Create immunogenicity positive and negative plots
## Filter for positive properties
iedb_faceted_pos <- filter(iedb_faceted_properties_nodup, Immunogenicity == "Positive")
iedb_faceted_neg <- filter(iedb_faceted_properties_nodup, Immunogenicity == "Negative")
## hist(count(as.factor(iedb_faceted_neg$top_allele_facet)))

## mw_peptide
pdf("boxplot_iedb_mw_peptide_perallele.pdf", height = 15, width = 15)
iedb_faceted_mw_peptide <- ggplot (iedb_faceted_properties_nodup, aes (x = Immunogenicity, y = mw_peptide)) + 
  labs(title = "Immunogenicity Vs. Molecular Weight", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_mw_peptide
dev.off()
iedb_faceted_mw_peptide

## mw_tcr_contact
pdf("boxplot_iedb_tcr_contact_perallele.pdf", height = 15, width = 15)
iedb_faceted_mw_tcr_contact <- ggplot (iedb_faceted_properties_nodup, aes (x = Immunogenicity, y = mw_tcr_contact)) + 
  labs(title = "Immunogenicity Vs. Molecular Weight of the TCR contact region", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_mw_tcr_contact
dev.off()
iedb_faceted_mw_tcr_contact

## hydroph_peptide
pdf("boxplot_iedb_hydroph_peptide_perallele.pdf", height = 15, width = 15)
iedb_faceted_hydroph_peptide <- ggplot (iedb_faceted_properties_nodup, aes (x = Immunogenicity, y = hydroph_peptide)) + 
  labs(title = "Immunogenicity Vs. Peptide Hydrophobicity", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_hydroph_peptide
dev.off()
iedb_faceted_hydroph_peptide


## hydroph_tcr_contact
pdf("boxplot_iedb_hydroph_tcr_contact_perallele.pdf", height = 15, width = 15)
iedb_faceted_hydroph_tcr_contact <- ggplot (iedb_faceted_properties_nodup, aes (x = Immunogenicity, y = hydroph_tcr_contact)) + 
  labs(title = "Immunogenicity Vs. Hydrophobicity of the TCR contact region", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_hydroph_tcr_contact
dev.off()
iedb_faceted_hydroph_tcr_contact

## charge_peptide
pdf("boxplot_iedb_charge_peptide_perallele.pdf", height = 15, width = 15)
iedb_faceted_charge_peptide <- ggplot (iedb_faceted_properties_nodup, aes (x = Immunogenicity, y = charge_peptide)) + 
  labs(title = "Immunogenicity Vs. Peptide Charge", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_charge_peptide
dev.off()
iedb_faceted_charge_peptide

## charge_tcr_contact
pdf("boxplot_iedb_charge_tcr_contact_perallele.pdf", height = 15, width = 15)
iedb_faceted_charge_tcr_contact <- ggplot (iedb_faceted_properties_nodup, aes (x = Immunogenicity, y = charge_tcr_contact)) + 
  labs(title = "Immunogenicity Vs. Net charge of the TCR Contact Region", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_charge_tcr_contact
dev.off()
iedb_faceted_charge_tcr_contact

## TAP_efficiency
pdf("boxplot_iedb_TAP_efficiency_perallele.pdf", height = 15, width = 15)
iedb_faceted_TAP_efficiency <- ggplot (iedb_faceted_properties_nodup, aes (x = Immunogenicity, y = TAP_efficiency)) + 
  labs(title = "Immunogenicity Vs. TAP Transport Efficiency", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_TAP_efficiency
dev.off()
iedb_faceted_TAP_efficiency

## Proteasomal_Cleavage
pdf("boxplot_iedb_Proteasomal_Cleavage_perallele.pdf", height = 15, width = 15)
iedb_faceted_Proteasomal_Cleavage <- ggplot (iedb_faceted_properties_nodup, aes (x = Immunogenicity, y = Proteasomal_Cleavage)) + 
  labs(title = "Immunogenicity Vs. Proteasomal Processing", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_Proteasomal_Cleavage
dev.off()
iedb_faceted_Proteasomal_Cleavage

## Affinity scores per Allele are not correct with current data since in peptides with multiple allele 
## only data for the best allele is provided. Recalculate using mhcflurry. ## Processing score is also affected by allele.

## Rename allele to best.allele
colnames(iedb_faceted_properties_nodup)
colnames(iedb_faceted_properties_nodup)[which(names(iedb_faceted_properties_nodup) == "allele")] <- "allele.genotype"
colnames(iedb_faceted_properties_nodup)[which(names(iedb_faceted_properties_nodup) == "peptide.x")] <- "peptide"
iedb_faceted_properties_nodup %>% mutate(allele = top_allele_facet) -> iedb_faceted_properties_nodup

## Export iedb_faceted_properties_nodup to run mhcflurry
write.table(iedb_faceted_properties_nodup, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200831_IEDB_peptide_properties_per_Allele/facet_allele_properties/mhcflurry_faceted_allelespecific/iedb_faceted_properties_nodup_tomhcflurry.csv", sep = ",", quote = F, dec = ".")

## Import mhcflurry calculation
iedb_faceted_properties_nodup_specific_allele <- read.table(file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200831_IEDB_peptide_properties_per_Allele/facet_allele_properties/mhcflurry_faceted_allelespecific/iedb_faceted_properties_nodup_tomhcflurry_out.csv", header = T, sep = ";", dec = ".")
colnames(iedb_faceted_properties_nodup_specific_allele)

"allele.genotype"                               "bestAllele_mhcflurry_affinity"                
[40] "bestAllele_mhcflurry_best_allele"              "bestAllele_mhcflurry_affinity_percentile"      "bestAllele_mhcflurry_processing_score"        
[43] "bestAllele_mhcflurry_presentation_score"       "allele"                                        "specific_allele_mhcflurry_affinity"           
[46] "specific_allele_mhcflurry_best_allele"         "specific_allele_mhcflurry_affinity_percentile" "specific_allele_mhcflurry_processing_score"   
[49] "specific_allele_mhcflurry_presentation_score" 


## specific_allele_mhcflurryaffinity
pdf("boxplot_iedb_mhcflurry_affinity_score_allele_specific_perallele.pdf", height = 15, width = 15)
iedb_faceted_mhcflurry_affinity <- ggplot (iedb_faceted_properties_nodup_specific_allele, aes (x = Immunogenicity, y = specific_allele_mhcflurry_affinity)) + 
  labs(title = "Immunogenicity Vs. HLA Binding Affinity (nM)", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
#iedb_faceted_mhcflurry_affinity_scale <- gg.gap(plot = iedb_faceted_mhcflurry_affinity, ylim = c(0,41000), tick_width = c(100,5000), segments = c(1000,20000), rel_heights = c(0.7, 0, 0.3))
#iedb_faceted_mhcflurry_affinity_scale
iedb_faceted_mhcflurry_affinity
dev.off()
iedb_faceted_mhcflurry_affinity

## specific_allele_mhcflurryaffinity_percentile
pdf("boxplot_iedb_mhcflurry_affinity_percentile_allele_specific_perallele.pdf", height = 15, width = 15)
iedb_faceted_mhcflurry_affinity_percentile <- ggplot (iedb_faceted_properties_nodup_specific_allele, aes (x = Immunogenicity, y = specific_allele_mhcflurry_affinity_percentile)) + 
  labs(title = "Immunogenicity Vs. HLA Binding Affinity (Percentile)", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_mhcflurry_affinity_percentile
dev.off()
iedb_faceted_mhcflurry_affinity_percentile

## specific_allele_mhcflurrypresentation_score
pdf("boxplot_iedb_mhcflurry_presentation_score_allele_specific_perallele.pdf", height = 15, width = 15)
iedb_faceted_mhcflurry_presentation_score <- ggplot (iedb_faceted_properties_nodup_specific_allele, aes (x = Immunogenicity, y = specific_allele_mhcflurry_presentation_score)) + 
  labs(title = "Immunogenicity Vs. MHCflurry presentation score", subtitle = "IEDB immunogenicity dataset - Per Allele") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_mhcflurry_presentation_score
dev.off()
iedb_faceted_mhcflurry_presentation_score

## specific_allele_mhcflurryprocessing_score
pdf("boxplot_iedb_mhcflurry_processing_score_allele_specific_perallele.pdf", height = 15, width = 15)
iedb_faceted_mhcflurry_processing_score <- ggplot (iedb_faceted_properties_nodup_specific_allele, aes (x = Immunogenicity, y = specific_allele_mhcflurry_processing_score)) + 
  labs(title = "Immunogenicity Vs. MHCflurry processing score", subtitle = "IEDB immunogenicity dataset - Allele specific") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  #stat_compare_means(method= "kruskal", label.x.npc = "center", method.args = list(alternative = "two.sided")) + 
  facet_wrap (~top_allele_facet, scales = "free_y")
iedb_faceted_mhcflurry_processing_score
dev.off()
iedb_faceted_mhcflurry_processing_score
