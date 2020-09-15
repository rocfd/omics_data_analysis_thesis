## IMMUNOGENICITY MODEL


## IMMUNOGENICITY PREDICTOR

## DESCRIPTION: logistic regression model based on epitope 
##immunogenicity and antigen processing and presentation plus
## physicochemical properties. The IEDB Immunogenicity dataset has
## been used as training data.

#### PERFORMANCE COMPARISON WITH STATE-OF-THE-ART METRICS
## Compare the performance of modelling MHCflurry Antigen 
## processing score and HLA binding affinity percentile.

## Libraries
library(dplyr)
library(glmnet)

# R Object : iedb_faceted_properties_nodup_specific_allele
View(iedb_faceted_properties_nodup_specific_allele)

## Export to location
#write.table(x = iedb_faceted_properties_nodup_specific_allele, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/iedb_faceted_dataset.csv", quote = F, sep = ";", dec = ".")

## UPLOAD DATASET MANUALLY!!
#iedb_file <- "path/to/file"
#read.table(iedb_file, header = T, sep = ";", dec = ".")


## Rename dataset
predictor_iedb_faceted <- iedb_faceted_properties_nodup_specific_allele
dim(predictor_iedb_faceted)
  colnames(predictor_iedb_faceted)
  ## Remove last register cause it contains blanks. Recheck peptide order to ensure no reordering ocurred. It didnt.
  predictor_iedb_faceted <- slice(predictor_iedb_faceted, 1:25836)
  glimpse(predictor_iedb_filter)
  
  ## Filter response variable (Immunogenicity)  & Predictors & Informative varibles (peptide, top_allele_facet, allele)
  
  predictor_iedb_filter <- select(predictor_iedb_faceted, peptide, top_allele_facet, Immunogenicity, mw_peptide, mw_tcr_contact, hydroph_peptide, hydroph_tcr_contact, charge_peptide, charge_tcr_contact, stab_peptide, stab_tcr_contact, TAP_efficiency, Proteasomal_Cleavage, allele, specific_allele_mhcflurry_affinity, specific_allele_mhcflurry_affinity_percentile, specific_allele_mhcflurry_processing_score, specific_allele_mhcflurry_presentation_score)
  dim(predictor_iedb_filter)
  colnames(predictor_iedb_filter)
  glimpse(predictor_iedb_filter)  
  

### EXTRACT IMMUNOGENICITY and PREDICTORS  

## Pull response variable as immunogenicity_iedb as vector
  immunogenicity_iedb <- pull(predictor_iedb_filter, Immunogenicity)
  class(immunogenicity_iedb) ## factor --> change to vector
  immunogenicity_iedb_vector <- as.vector(immunogenicity_iedb)
  setequal(immunogenicity_iedb, immunogenicity_iedb_vector)
  length(immunogenicity_iedb_vector)
  
  ## Pull predictors as predictors_iedb as matrix
  #predictors_iedb <- pull(predictor_iedb_filter, mw_peptide, mw_tcr_contact, hydroph_peptide, hydroph_tcr_contact, charge_peptide, charge_tcr_contact, stab_peptide, stab_tcr_contact, TAP_efficiency, Proteasomal_Cleavage, specific_allele_mhcflurry_affinity, specific_allele_mhcflurry_affinity_percentile, specific_allele_mhcflurry_processing_score, specific_allele_mhcflurry_presentation_score)
  ## pull only works for a single column
  
## Pull predictors one by one and merge into a predictors matrix
  mw_peptide_iedb <- pull(predictor_iedb_filter, mw_peptide)
  mw_tcr_contact_iedb <- pull(predictor_iedb_filter, mw_tcr_contact)
  hydroph_peptide_iedb <- pull(predictor_iedb_filter, hydroph_peptide) 
  hydroph_tcr_contact_iedb <- pull(predictor_iedb_filter, hydroph_tcr_contact)
  charge_peptide_iedb <- pull(predictor_iedb_filter, charge_peptide)
  charge_tcr_contact_iedb <- pull(predictor_iedb_filter, charge_tcr_contact)
  stab_peptide_iedb <- pull(predictor_iedb_filter, stab_peptide)
  stab_tcr_contact_iedb <- pull(predictor_iedb_filter, stab_tcr_contact)
  TAP_efficiency_iedb <- pull(predictor_iedb_filter, TAP_efficiency)
  Proteasomal_Cleavage_iedb <- pull(predictor_iedb_filter, Proteasomal_Cleavage)
  specific_allele_mhcflurry_affinity_iedb <- pull(predictor_iedb_filter, specific_allele_mhcflurry_affinity)
  specific_allele_mhcflurry_affinity_percentile_iedb <- pull(predictor_iedb_filter, specific_allele_mhcflurry_affinity_percentile)
  specific_allele_mhcflurry_processing_score_iedb <- pull(predictor_iedb_filter, specific_allele_mhcflurry_processing_score)
  specific_allele_mhcflurry_presentation_score_iedb <- pull(predictor_iedb_filter, specific_allele_mhcflurry_presentation_score)
  
## Matrix
  matrix_iedb <- cbind(mw_peptide_iedb, mw_tcr_contact_iedb, hydroph_peptide_iedb, hydroph_tcr_contact_iedb, charge_peptide_iedb, charge_tcr_contact_iedb, stab_peptide_iedb, stab_tcr_contact_iedb, TAP_efficiency_iedb, Proteasomal_Cleavage_iedb, specific_allele_mhcflurry_affinity_iedb, specific_allele_mhcflurry_affinity_percentile_iedb, specific_allele_mhcflurry_processing_score_iedb, specific_allele_mhcflurry_presentation_score_iedb)
  class(matrix_iedb)
  View(matrix_iedb)
  
  
  
  
  
  #####  GLM PACKAGE!!!!
  
  ## Standardize continous variables
  predictor_iedb_rescale <- predictor_iedb_filter %>%
    mutate_if(is.numeric, funs(as.numeric(scale(.))))
  head(predictor_iedb_rescale)
  
  ## plot test
  ggplot(predictor_iedb_rescale, aes(x = specific_allele_mhcflurry_presentation_score, y = hydroph_peptide)) +
    geom_point(aes(color = Immunogenicity),
               size = 0.5) +
    stat_smooth(method = 'lm',
                formula = y~poly(x, 2),
                se = TRUE,
                aes(color = Immunogenicity)) +
    theme_classic()
  
  library(GGally)
  # Convert data to numeric
  corr <- data.frame(lapply(predictor_iedb_filter, as.integer))
  # Plot the graph
  ggcorr (corr,
  method = c("pairwise", "spearman"),
  nbreaks = 6,
  hjust = 0.8,
  label = TRUE,
  label_size = 3,
  color = "grey50")
  
  ## Filter to response and predictors only
  
  predictor_response_iedb <- select(predictor_iedb_filter, Immunogenicity, mw_peptide, mw_tcr_contact, hydroph_peptide, hydroph_tcr_contact, charge_peptide, charge_tcr_contact, stab_peptide, stab_tcr_contact, TAP_efficiency, Proteasomal_Cleavage, specific_allele_mhcflurry_affinity, specific_allele_mhcflurry_affinity_percentile, specific_allele_mhcflurry_processing_score, specific_allele_mhcflurry_presentation_score)
  
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
  data_train <- create_train_test(predictor_response_iedb, 0.8, train = TRUE)
  data_test <- create_train_test(predictor_response_iedb, 0.8, train = FALSE)
  dim(data_train)
  dim(data_test)
  
  ## Export sets
  write.table(x = data_train, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/training_test_data/predictor_data_train.csv", sep = ";", dec = ".", quote = F)
  write.table(x = data_test, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/training_test_data/predictor_data_test.csv", sep = ";", dec = ".", quote = F)
  
  ## Build the glm model
  formula <- data_train$Immunogenicity~.
  logit <- glm(formula = formula, family = "binomial", data = data_train)
  summary(logit)
  
  ## Predict on test data
  predict <- predict(logit, data_test, type = 'response')
  table_mat <- table(data_test$Immunogenicity, predict > 0.5)
  table_mat
  
  ## Accuracy test
  accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
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
  prec <- precision(table_mat)
  prec
  rec <- recall(table_mat)
  rec
  
  ## F score
  f1 <- 2 * ((prec * rec) / (prec + rec))
  f1
  
  ## ROC curve
  library(ROCR)
  ROCRpred <- prediction(predict, data_test$Immunogenicity)
  ROCRpred
  ROCRperf <- performance(ROCRpred, 'tpr', 'fpr')
  ROCRperf
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_predictor_model/predictor_ROC_curve_lab.pdf")
  ROC_bsc <- plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2, 1.7),
               main = "ROC Curve \n PredIG Vs. Immunogenicity",
                sub = "Generalized Linear Model (glm) - Logistic Regression")
dev.off()
  
  ### Global F test
  summary(logit)
  
#### PERFORMANCE COMPARISON VS MHCFLURRY
#### flurry_model --> antigen presentation score
#### flurry_percentile_model --> HLA binding affinity percentile


  ### Compare VS MHCflurry model alone
  
  ## Filter to response and predictors only
  
  #predictor_response_iedb <- select(predictor_iedb_filter, Immunogenicity, mw_peptide, mw_tcr_contact, hydroph_peptide, hydroph_tcr_contact, charge_peptide, charge_tcr_contact, stab_peptide, stab_tcr_contact, TAP_efficiency, Proteasomal_Cleavage, specific_allele_mhcflurry_affinity, specific_allele_mhcflurry_affinity_percentile, specific_allele_mhcflurry_processing_score, specific_allele_mhcflurry_presentation_score)
  mhcflurry_response_iedb <- select(predictor_iedb_filter, Immunogenicity, specific_allele_mhcflurry_presentation_score)

  ## Training and test
  data_train_flurry <- create_train_test(mhcflurry_response_iedb, 0.8, train = TRUE)
  data_test_flurry <- create_train_test(mhcflurry_response_iedb, 0.8, train = FALSE)
  dim(data_train_flurry)
  dim(data_test_flurry)
  
  ## Export sets
  write.table(x = data_train_flurry, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/training_test_data_mhcflurry_model/mhcflurry_data_train.csv", sep = ";", dec = ".", quote = F)
  write.table(x = data_test_flurry, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/training_test_data_mhcflurry_model/mhcflurry_data_test.csv", sep = ";", dec = ".", quote = F)
  
  
  
  ## Build the glm model
  formula_flurry <- data_train_flurry$Immunogenicity~.
  logit_flurry <- glm(formula = formula_flurry, family = "binomial", data = data_train_flurry)
  summary(logit_flurry)
  
  ## Predict on test data
  predict_flurry <- predict(logit_flurry, data_test_flurry, type = 'response')
  table_mat_flurry <- table(data_test_flurry$Immunogenicity, predict > 0.5)
  table_mat_flurry
  
  ## Accuracy test
  accuracy_Test_flurry <- sum(diag(table_mat_flurry)) / sum(table_mat_flurry)
  accuracy_Test_flurry #0.7927632
  

  ## Compute precision and recall
  prec_flurry <- precision(table_mat_flurry)
  prec_flurry
  rec_flurry <- recall(table_mat_flurry)
  rec_flurry
  
  ## F score
  f1_flurry <- 2 * ((prec_flurry * rec_flurry) / (prec_flurry + rec_flurry))
  f1_flurry
  
  ## ROC curve
  library(ROCR)
  ROCRpred_flurry <- prediction(predict_flurry, data_test_flurry$Immunogenicity)
  ROCRperf_flurry <- performance(ROCRpred_flurry, 'tpr', 'fpr')
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_predictor_model/mhcflurry_ROC_curve.pdf")
  ROC_flurry <- plot(ROCRperf_flurry, colorize = TRUE, text.adj = c(-0.2, 1.7),
                  main = "ROC Curve \n MHCflurry Antigen Presentation Score Vs. Immunogenicity",
                  sub = "Generalized Linear Model (glm) - Logistic Regression")
  dev.off()
  
  summary(logit_flurry)
  
  ## Multiplot both ROC curves
  rocs <- grid.arrange(ROC_flurry, ROC_bsc)
  
  
  ### Compare VS MHCflurry BINDING AFFINITY PERCENTILE
  
  ## Filter to response and predictors only
  
  #predictor_response_iedb <- select(predictor_iedb_filter, Immunogenicity, mw_peptide, mw_tcr_contact, hydroph_peptide, hydroph_tcr_contact, charge_peptide, charge_tcr_contact, stab_peptide, stab_tcr_contact, TAP_efficiency, Proteasomal_Cleavage, specific_allele_mhcflurry_affinity, specific_allele_mhcflurry_affinity_percentile, specific_allele_mhcflurry_processing_score, specific_allele_mhcflurry_presentation_score)
  mhcflurry_percentile_response_iedb <- select(predictor_iedb_filter, Immunogenicity, specific_allele_mhcflurry_affinity_percentile)
  
  ## Training and test
  data_train_flurry_percentile <- create_train_test(mhcflurry_percentile_response_iedb, 0.8, train = TRUE)
  data_test_flurry_percentile <- create_train_test(mhcflurry_percentile_response_iedb, 0.8, train = FALSE)
  dim(data_train_flurry_percentile)
  dim(data_test_flurry_percentile)
  
  ## Export sets
  write.table(x = data_train_flurry_percentile, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/training_test_data_mhcflurry_model/mhcflurry_percentile_data_train.csv", sep = ";", dec = ".", quote = F)
  write.table(x = data_test_flurry_percentile, file = "/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/training_test_data_mhcflurry_model/mhcflurry_percentile_data_test.csv", sep = ";", dec = ".", quote = F)
  
  
  
  ## Build the glm model
  formula_flurry_percentile <- data_train_flurry_percentile$Immunogenicity~.
  logit_flurry_percentile <- glm(formula = formula_flurry_percentile, family = "binomial", data = data_train_flurry_percentile)
  summary(logit_flurry_percentile)
  
  ## Predict on test data
  predict_flurry_percentile <- predict(logit_flurry_percentile, data_test_flurry_percentile, type = 'response')
  table_mat_flurry_percentile <- table(data_test_flurry_percentile$Immunogenicity, predict_flurry_percentile > 0.5)
  table_mat_flurry_percentile
  
  ## Accuracy test
  accuracy_Test_flurry_percentile <- sum(diag(table_mat_flurry_percentile)) / sum(table_mat_flurry_percentile)
  accuracy_Test_flurry_percentile #0.7979876
  
  
  ## Compute precision and recall
  prec_flurry_percentile <- precision(table_mat_flurry_percentile)
  prec_flurry_percentile #0.4285714
  rec_flurry_percentile <- recall(table_mat_flurry_percentile)
  rec_flurry_percentile #0.005758157
  
  ## F score
  f1_flurry_percentile <- 2 * ((prec_flurry_percentile * rec_flurry_percentile) / (prec_flurry_percentile + rec_flurry_percentile))
  f1_flurry_percentile #0.01136364
  
  ## ROC curve
  library(ROCR)
  ROCRpred_flurry_percentile <- prediction(predict_flurry_percentile, data_test_flurry_percentile$Immunogenicity)
  ROCRperf_flurry_percentile <- performance(ROCRpred_flurry_percentile, 'tpr', 'fpr')
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_predictor_model/mhcflurry_percentile_ROC_curve.pdf")
  ROC_flurry_percentile <- plot(ROCRperf_flurry_percentile, colorize = TRUE, text.adj = c(-0.2, 1.7),
                     main = "ROC Curve \n MHCflurry HLA Binding Affinity (Percentile) Vs. Immunogenicity",
                     sub = "Generalized Linear Model (glm) - Logistic Regression",
                     tag = "B")
  dev.off()
  
  summary(logit_flurry_percentile)
  
  ## Multiplot both ROC curves
  rocs <- grid.arrange(ROC_flurry_percentile, ROC_bsc)
  
