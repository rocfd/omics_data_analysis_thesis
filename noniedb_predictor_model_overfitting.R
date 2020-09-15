## IMMUNOGENICITY MODELS

## Assess the overfitting effects on the IEDB model
## by using an independently validated dataset of immunogenicity.

## Libraries
library(dplyr)
library(glmnet)

# R Object : iedb_faceted_properties_nodup_specific_allele
#View(iedb_faceted_properties_nodup_specific_allele)

## Export to location
#write.table(x = iedb_faceted_properties_nodup_specific_allele, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/iedb_faceted_dataset.csv", quote = F, sep = ";", dec = ".")

## UPLOAD DATASET MANUALLY!!
#predictor_non_iedb_faceted <- "path/to/file"
#read.table(iedb_file, header = T, sep = ";", dec = ".")
  
  #####  GLM PACKAGE!!!!
  ## NON_IEDB dataset
  ## non_iedb_deconvoluted_all
  
  ## rename
  predictor_noniedb <- non_iedb_deconvoluted_all
  dim(predictor_noniedb) # 8797
  write.table(x = predictor_noniedb, file = "non_iedb_deconvoluted_all.csv", quote = F, sep = ";", dec = ".")
  glimpse(predictor_noniedb)
  table(predictor_noniedb$Immunogenicity)
  
  ## Remove blanks and NAs from MHCflurry data
  predictor_noniedb_filter <- predictor_noniedb %>%
    filter(mhcflurry_best_allele != "")
  dim(predictor_noniedb_filter) #6926
  glimpse(predictor_noniedb_filter)
  View(predictor_noniedb_filter)
  table(predictor_noniedb_filter$mhcflurry_best_allele)

  
  ## Remove NA values
  #predictor_noniedb_nona <- as.data.frame(lapply(predictor_noniedb, na.omit))
  #dim(predictor_noniedb_filter_nona)  
  #new_dataframe = as.data.frame(lapply(df1, na.omit))
  #predictor_noniedb_nona <- lapply(predictor_noniedb, na.omit)
  #class(predictor_noniedb_nona)
  #as.data.frame(predictor_noniedb_nona)
  
  #fn2 <- function (d) {
  #  for (i in 1:nrow(d)) {
  #    if (any(is.na(d[i,]))) {
  #      d <- d[-i, ]
  #    }
  #  }
  #  d # return d
  #}
  
  #predictor_noniedb_nona <- fn2(predictor_noniedb_select)
  #glimpse(predictor_noniedb_nona)
  #View(predictor_noniedb_nona)
  
  ## SELECT VARIABLES
  predictor_noniedb_select <- select(predictor_noniedb_filter, Immunogenicity, mw_peptide, mw_tcr_contact, hydroph_peptide, hydroph_tcr_contact, charge_peptide, charge_tcr_contact, stab_peptide, stab_tcr_contact, TAP, Cle, mhcflurry_affinity, mhcflurry_affinity_percentile, mhcflurry_processing_score, mhcflurry_presentation_score)
  glimpse(predictor_noniedb_select)
  table(predictor_noniedb_select$Immunogenicity)

  ## Standardize continous variables
  ##predictor_noniedb_rescale <- predictor_noniedb %>%
   # mutate_if(is.numeric, funs(as.numeric(scale(.))))
 #glimpse(predictor_noniedb_rescale)
 
 
  ## plot test
  ggplot(predictor_noniedb_select, aes(x = mhcflurry_presentation_score, y = hydroph_peptide)) +
    geom_point(aes(color = Immunogenicity),
               size = 0.5) +
    stat_smooth(method = 'lm',
                formula = y~poly(x, 2),
                se = TRUE,
                aes(color = Immunogenicity)) +
    theme_classic()
  
  library(GGally)
  # Convert data to numeric
  corr <- data.frame(lapply(predictor_noniedb_filter, as.integer))
  # Plot the graph
  ggcorr (corr,
  method = c("pairwise", "spearman"),
  nbreaks = 6,
  hjust = 0.8,
  label = TRUE,
  label_size = 3,
  color = "grey50")
  
  ## Filter to response and predictors only
  
 ## predictor_response_iedb <- select(predictor_iedb_filter, Immunogenicity, mw_peptide, mw_tcr_contact, hydroph_peptide, hydroph_tcr_contact, charge_peptide, charge_tcr_contact, stab_peptide, stab_tcr_contact, TAP_efficiency, Proteasomal_Cleavage, specific_allele_mhcflurry_affinity, specific_allele_mhcflurry_affinity_percentile, specific_allele_mhcflurry_processing_score, specific_allele_mhcflurry_presentation_score)
  
  ### CREATE TEST AND TRAIN SETS
  ## Function
  set.seed(1234)
  create_train_test <- function(data, size = 0.8, train = TRUE) {
    n_row = nrow(data)
    total_row = size * n_row
    train_sample <- 1: total_row
    if (train == TRUE) {
      return (data[train_sample, ])
    } else {
      return (data[-train_sample, ])
    }
  }
  
  ## Training and test
  data_train_noniedb <- create_train_test(predictor_noniedb_select, 0.8, train = TRUE)
  data_test_noniedb <- create_train_test(predictor_noniedb_select, 0.8, train = FALSE)
  dim(data_train_noniedb) #5540
  dim(data_test_noniedb)  #1386
  table(data_train_noniedb$Immunogenicity) ## pos 1213 neg 4327
  table(data_test_noniedb$Immunogenicity) ## pos 323 neg 1063
  View(data_test_noniedb)
  
  ## Remove NA from training.
  
  ## Export sets
  write.table(x = data_train_noniedb, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control/non_iedb_train_test/noniedb_predictor_data_train.csv", sep = ";", dec = ".", quote = F)
  write.table(x = data_test_noniedb, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/overfitting_control/non_iedb_train_test/noniedb_predictor_data_test.csv", sep = ";", dec = ".", quote = F)
  
  ## Build the glm model
  formula_noniedb <- data_train_noniedb$Immunogenicity~.
  logit_noniedb <- glm(formula = formula_noniedb, family = "binomial", data = data_train_noniedb)
  summary(logit_noniedb)
  
  ## Predict on test data
  predict_noniedb <- predict(logit_noniedb, data_test_noniedb, type = 'response')
  table_mat_noniedb <- table(data_test_noniedb$Immunogenicity, predict_noniedb > 0.5)
  class(data_test_noniedb$Immunogenicity)
  data_test_noniedb$Immunogenicity
  class(table_mat_noniedb)
  table_mat_noniedb
  
  ## Accuracy test
  accuracy_Test <- sum(diag(table_mat_noniedb)) / sum(table_mat_noniedb)
  accuracy_Test
  
  ## Precision variable (surrogate for sensibility)
  precision <- function(matrix) {
    # True positive
    tp <- matrix[2, 2]
    # false positive
    fp <- matrix[1, 2]
    return (tp / (tp + fp))
  }
  
  ## Recall function (surrogate for specificity)
  recall <- function(matrix) {
    # true positive
    tp <- matrix[2, 2]# false positive
    fn <- matrix[2, 1]
    return (tp / (tp + fn))
  }
  
  ## Compute precision and recall
  prec_noniedb <- precision(table_mat_noniedb)
  prec_noniedb
  rec_noniedb <- recall(table_mat_noniedb)
  rec_noniedb
  
  ## F score
  f1_noniedb <- 2 * ((prec_noniedb * rec_noniedb) / (prec_noniedb + rec_noniedb))
  f1_noniedb
  
  ## ROC curve
  library(ROCR)
  ROCRpred_noniedb <- prediction(predict_noniedb, data_test_noniedb$Immunogenicity)
  #ROCRpred <- prediction(predict_noniedb, data_test_noniedb$Immunogenicity)
  ROCRpred_noniedb
  ROCRperf_noniedb <- performance(ROCRpred_noniedb, 'tpr', 'fpr')
  ROCRperf_noniedb
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_predictor_model/predictor_ROC_curve_noniedb.pdf")
  ROC_bsc_noniedb <- plot(ROCRperf_noniedb, colorize = TRUE, text.adj = c(-0.2, 1.7),
               main = "ROC Curve \n BSC Predictor Vs. Indep Immunogenicity",
                sub = "Generalized Linear Model (glm) - Logistic Regression")
dev.off()
  
  ### Global F test
  summary(logit)
  
  
  ##### TEST noniedb data on IEDB model
  
  ## Build the glm model
  formula_noniedb <- data_train_noniedb$Immunogenicity~.
  logit_noniedb <- glm(formula = formula_noniedb, family = "binomial", data = data_train_noniedb)
  summary(logit_noniedb)
  
  ## Rename noniedb coluns to match iedbs
  names(predictor_noniedb_select)[names(predictor_noniedb_select) == 'TAP'] <- 'TAP_efficiency'
  names(predictor_noniedb_select)[names(predictor_noniedb_select) == 'Cle'] <- 'Proteasomal_Cleavage'
  names(predictor_noniedb_select)[names(predictor_noniedb_select) == 'mhcflurry_affinity'] <- 'specific_allele_mhcflurry_affinity'
  names(predictor_noniedb_select)[names(predictor_noniedb_select) == 'mhcflurry_affinity_percentile'] <- 'specific_allele_mhcflurry_affinity_percentile'
  names(predictor_noniedb_select)[names(predictor_noniedb_select) == 'mhcflurry_processing_score'] <- 'specific_allele_mhcflurry_processing_score'
  names(predictor_noniedb_select)[names(predictor_noniedb_select) == 'mhcflurry_presentation_score'] <- 'specific_allele_mhcflurry_presentation_score'

  
  ## Predict on test data
  predict_noniedb_ON <- predict(logit, data_test_noniedb, type = 'response')
  ## rename columns
  table_mat_noniedb_ON <- table(data_test_noniedb$Immunogenicity, predict_noniedb_ON > 0.5)
  class(data_test_noniedb$Immunogenicity)
  data_test_noniedb$Immunogenicity
  class(table_mat_noniedb_ON)
  table_mat_noniedb_ON
  
  ## Accuracy test
  accuracy_Test <- sum(diag(table_mat_noniedb_ON)) / sum(table_mat_noniedb_ON)
  accuracy_Test
  
  ## Precision variable (surrogate for sensibility)
  precision <- function(matrix) {
    # True positive
    tp <- matrix[2, 2]
    # false positive
    fp <- matrix[1, 2]
    return (tp / (tp + fp))
  }
  
  ## Recall function (surrogate for specificity)
  recall <- function(matrix) {
    # true positive
    tp <- matrix[2, 2]# false positive
    fn <- matrix[2, 1]
    return (tp / (tp + fn))
  }
  
  ## Compute precision and recall
  prec_noniedb_ON <- precision(table_mat_noniedb_ON)
  prec_noniedb_ON
  rec_noniedb_ON <- recall(table_mat_noniedb_ON)
  rec_noniedb_ON
  
  ## F score
  f1_noniedb_ON <- 2 * ((prec_noniedb_ON * rec_noniedb_ON) / (prec_noniedb_ON + rec_noniedb_ON))
  f1_noniedb_ON
  
  ## ROC curve
  library(ROCR)
  ROCRpred_noniedb_ON <- prediction(predict_noniedb_ON, data_test_noniedb$Immunogenicity)
  #ROCRpred <- prediction(predict_noniedb, data_test_noniedb$Immunogenicity)
  ROCRpred_noniedb_ON
  ROCRperf_noniedb_ON <- performance(ROCRpred_noniedb_ON, 'tpr', 'fpr')
  ROCRperf_noniedb_ON
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_predictor_model/predictor_ROC_curve_overfitting.pdf")
  ROC_bsc_noniedb_ON <- plot(ROCRperf_noniedb_ON, colorize = TRUE, text.adj = c(-0.2, 1.7),
                          main = "ROC Curve (Overfitting test) \n BSC Predictor Vs. Indep Immunogenicity",
                          sub = "Generalized Linear Model (glm) - Logistic Regression")
  dev.off()
  
  ### Global F test
  summary(logit)
  
  ##### TEST noniedb data on IEDB MHCflurry model
  
  ## Predict on test data
  predict_noniedb_mhcflurry <- predict(logit_flurry, data_test_noniedb, type = 'response')
  predict_noniedb_mhcflurry
  ## rename columns
  table_mat_noniedb_mhcflurry <- table(data_test_noniedb$Immunogenicity, predict_noniedb_mhcflurry > 0.5)
  class(data_test_noniedb$Immunogenicity)
  data_test_noniedb$Immunogenicity
  class(table_mat_noniedb_mhcflurry)
  table_mat_noniedb_mhcflurry
  
  ## Accuracy test
  accuracy_Test <- sum(diag(table_mat_noniedb_mhcflurry)) / sum(table_mat_noniedb_mhcflurry)
  accuracy_Test
  
  ## Compute precision and recall
  prec_noniedb_mhcflurry <- precision(table_mat_noniedb_mhcflurry)
  prec_noniedb_mhcflurry
  rec_noniedb_mhcflurry <- recall(table_mat_noniedb_mhcflurry)
  rec_noniedb_mhcflurry
  
  ## F score
  f1_noniedb_mhcflurry <- 2 * ((prec_noniedb_mhcflurry * rec_noniedb_mhcflurry) / (prec_noniedb_mhcflurry + rec_noniedb_mhcflurry))
  f1_noniedb_mhcflurry
  
  ## ROC curve
  library(ROCR)
  ROCRpred_noniedb_mhcflurry <- prediction(predict_noniedb_mhcflurry, data_test_noniedb$Immunogenicity)
  #ROCRpred <- prediction(predict_noniedb, data_test_noniedb$Immunogenicity)
  ROCRpred_noniedb_mhcflurry
  ROCRperf_noniedb_mhcflurry <- performance(ROCRpred_noniedb_mhcflurry, 'tpr', 'fpr')
  ROCRperf_noniedb_mhcflurry
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_predictor_model/predictor_ROC_curve_overfitting_mhcflurry.pdf")
  ROC_bsc_noniedb_mhcflurry <- plot(ROCRperf_noniedb_mhcflurry, colorize = TRUE, text.adj = c(-0.2, 1.7),
                             main = "ROC Curve (Overfitting test) \n BSC Predictor Vs. Indep Immunogenicity",
                             sub = "Generalized Linear Model (glm) - Logistic Regression")
  dev.off()
  
  ### Global F test
  summary(logit)
  
  
