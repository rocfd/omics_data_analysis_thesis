## IEDB Peptide Properties -  FACETED for MHCFlurryscores - NOT_FACETED for other scores

## Libraries
library(ggforce)
library(ggplot2)
library(ggpubr)
library(gg.gap)


## not faceted dataset
## iedb_TAP_mhcflurryrestricted_nodup_recalc

## faceted dataset
## iedb_faceted_properties_nodup_specific_allele

## wd
setwd("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/C_REPORT/figures/plots")

## PLOTTING - PHYSICOCHEMICAL PROPERTIES - NOT FACETED DATASET
## iedb_TAP_mhcflurryrestricted_nodup_recalc


## MW PEPTIDE
pdf("boxplot_iedb_mw_peptide_overall_report.pdf", height = 5, width = 5)
iedb_overall_mw_peptide <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = mw_peptide)) + 
  labs(x = "Immunogenicity", y = "Epitope Molecular Weight", tag = "B") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 1575, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_mw_peptide
dev.off()
iedb_overall_mw_peptide

## MW TCR CONTACT
pdf("boxplot_iedb_mw_tcr_contact_overall_report.pdf", height = 5, width = 5)
iedb_overall_mw_tcr_contact <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = mw_tcr_contact)) + 
  labs(x = "Immunogenicity", y = "TCR contact region Molecular Weight", tag = "B") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 575, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_mw_tcr_contact
dev.off()
iedb_overall_mw_tcr_contact

## HYDROPH PEPTIDE
pdf("boxplot_iedb_hydroph_peptide_overall_report.pdf", height = 5, width = 5)
iedb_overall_hydroph_peptide <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = hydroph_peptide)) + 
  labs(x = "Immunogenicity", y = "Epitope Hydrophobicity", tag = "A") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 5, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_hydroph_peptide
dev.off()
iedb_overall_hydroph_peptide

## HYDROPH TCR CONTACT
pdf("boxplot_iedb_hydroph_tcr_contact_overall_report.pdf", height = 5, width = 5)
iedb_overall_hydroph_tcr_contact <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = hydroph_tcr_contact)) + 
  labs(x = "Immunogenicity", y = "TCR Contact Region Hydrophobicity", tag = "A") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 5.5, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_hydroph_tcr_contact
dev.off()
iedb_overall_hydroph_tcr_contact

## CHARGE PEPTIDE
pdf("boxplot_iedb_charge_peptide_overall_report.pdf", height = 5, width = 5)
iedb_overall_charge_peptide <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = charge_peptide)) + 
  labs(x = "Immunogenicity", y = "Epitope Net Charge", tag = "C") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 7.5, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_charge_peptide
dev.off()
iedb_overall_charge_peptide

## CHARGE TCR CONTACT
pdf("boxplot_iedb_charge_tcr_contact_overall_report.pdf", height = 5, width = 5)
iedb_overall_charge_tcr_contact <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = charge_tcr_contact)) + 
  labs(x = "Immunogenicity", y = "TCR Contact Region Net Charge", tag = "C") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 3.8, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_charge_tcr_contact
dev.off()
iedb_overall_charge_tcr_contact

## STAB PEPTIDE
pdf("boxplot_iedb_stab_peptide_overall_report.pdf", height = 5, width = 5)
iedb_overall_stab_peptide <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = stab_peptide)) + 
  labs(x = "Immunogenicity", y = "Epitope Stability", tag = "D") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 330, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_stab_peptide
dev.off()
iedb_overall_stab_peptide

## STAB PEPTIDE
pdf("boxplot_iedb_stab_tcr_overall_report.pdf", height = 5, width = 5)
iedb_overall_stab_tcr <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = stab_tcr_contact)) + 
  labs(x = "Immunogenicity", y = "TCR Contact Region Stability", tag = "D") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 450, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_stab_tcr
dev.off()
iedb_overall_stab_tcr

## TAP TRANSPORT
pdf("boxplot_iedb_tap_transport_overall_report.pdf", height = 5, width = 5)
iedb_overall_tap <- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = TAP_efficiency)) + 
  labs(x = "Immunogenicity", y = "TAP Transport Efficiency", tag = "B") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 4.5, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_tap
dev.off()
iedb_overall_tap


## PROTEASOMAL CLEAVAGE
pdf("boxplot_iedb_proteasomal_cleavage_overall_report.pdf", height = 5, width = 5)
iedb_overall_prot<- ggplot (iedb_TAP_mhcflurryrestricted_nodup_recalc, aes (x = Immunogenicity, y = Proteasomal_Cleavage)) + 
  labs(x = "Immunogenicity", y = "Proteasomal Processing", tag = "A") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 1.1, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_prot
dev.off()
iedb_overall_prot

##########
##########
## MHCFLURRY PROPERTIES - FACETED DATASET - ALLELE SPECIFIC VS BEST ALLELE

## mhcflurrypresentation_score_overall - ALLELE SPECIFIC
pdf("boxplot_iedb_mhcflurry_presentation_score_Allelespecific_overall_report.pdf", height = 5, width = 5)
iedb_overall_mhcflurry_presentation_score_specific <- ggplot (iedb_faceted_properties_nodup_specific_allele, aes (x = Immunogenicity, y = specific_allele_mhcflurry_presentation_score)) + 
  labs(x = "Immunogenicity", y = "Antigen Presentation Score", tag = "F") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) + 
  stat_compare_means(method= "wilcox.test", label.y = 1.1, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_mhcflurry_presentation_score_specific
dev.off()
iedb_overall_mhcflurry_presentation_score_specific


## mhcflurry_processing_score_perAllele_overall - ALLELE SPECIFIC
pdf("boxplot_iedb_mhcflurry_processing_score_allele_specific_overall_report.pdf", height = 5, width = 5)
iedb_overall_mhcflurry_processing_score_allele_specific <- ggplot (iedb_faceted_properties_nodup_specific_allele, aes (x = Immunogenicity, y = specific_allele_mhcflurry_processing_score)) + 
  labs(x = "Immunogenicity", y = "Antigen Processing Score", tag = "C") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  stat_compare_means(method= "wilcox.test", label.y = 1.1, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_mhcflurry_processing_score_allele_specific
dev.off()
iedb_overall_mhcflurry_processing_score_allele_specific


## AFFINITY

## mhcflurry_affinity_score_AlleleSpecific_overall
pdf("boxplot_iedb_mhcflurry_affinity_allele_specific_overall_report.pdf", height = 5, width = 5)
iedb_overall_mhcflurry_affinity_allele_specific <- ggplot (iedb_faceted_properties_nodup_specific_allele, aes (x = Immunogenicity, y = specific_allele_mhcflurry_affinity)) + 
  labs(x = "Immunogenicity", y = "HLA Binding Affinity (nM)") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  stat_compare_means(method= "wilcox.test", label.y = 40500, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_mhcflurry_affinity_allele_specific_split <- gg.gap(plot = iedb_overall_mhcflurry_affinity_allele_specific, ylim = c(0, 41500), segments = c(1000,35000), tick_width = c(100,1000), rel_heights = c(0.8, 0.0001, 0.2), )
iedb_overall_mhcflurry_affinity_allele_specific_split
dev.off()
iedb_overall_mhcflurry_affinity_allele_specific

## AFFINITY PERCENTILE

## mhcflurry_affinity_percentile_AlleleSpecific_overall
pdf("boxplot_iedb_mhcflurry_affinity_percentile_allele_specific_overall_report.pdf", height = 5, width = 5)
iedb_overall_mhcflurry_affinity_percentile_allele_specific <- ggplot (iedb_faceted_properties_nodup_specific_allele, aes (x = Immunogenicity, y = specific_allele_mhcflurry_affinity_percentile)) + 
  labs(x = "Immunogenicity", y = "HLA Binding Affinity (Percentile)") +
  geom_sina (size = 0.1) + 
  geom_boxplot (alpha = 0.5, outlier.shape = NA) + 
  geom_signif(comparisons = list(c("Negative", "Positive")), map_signif_level=TRUE, textsize = 3) +
  stat_compare_means(method= "wilcox.test", label.y = 90, label.y.npc = "top", label.x.npc = "center", method.args = list(alternative = "two.sided"))
#facet_wrap (~top_allele_facet, scales = "free_y")
iedb_overall_mhcflurry_affinity_allele_percentile_specific_split <- gg.gap(plot = iedb_overall_mhcflurry_affinity_percentile_allele_specific, ylim = c(0, 100), segments = c(4,40), tick_width = c(2,20), rel_heights = c(0.8, 0.0001, 0.2), )
iedb_overall_mhcflurry_affinity_allele_percentile_specific_split
dev.off()
iedb_overall_mhcflurry_affinity_allele_percentile_specific_split
