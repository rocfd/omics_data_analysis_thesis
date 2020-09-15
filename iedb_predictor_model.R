## IMMUNOGENICITY PREDICTOR

## DESCRIPTION: logistic regression model based on epitope 
##immunogenicity and antigen processing and presentation plus
## physicochemical properties. The IEDB Immunogenicity dataset has
## been used as training data.

## Libraries
library(dplyr)
library(glmnet)

# R Object : iedb_faceted_properties_nodup_specific_allele
View(iedb_faceted_properties_nodup_specific_allele)

## Export to location
#write.table(x = iedb_faceted_properties_nodup_specific_allele, file = "path/to/dataset.csv", quote = F, sep = ";", dec = ".")

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
  View(predictor_iedb_filter)  
  

## EXPLORE THE DISTRIBUTION OF CONTINOUS VARIABLES
  continuous <-select_if(predictor_iedb_filter, is.numeric)
  summary(continuous)
  
  # Histogram with kernel density curve
  library(ggplot2)
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_mw_peptide.pdf")
  ggplot(continuous, aes(x = mw_peptide)) +
    labs(title = "Peptide Molecular Weight Distribution", subtitle = "IEDB Immunogenicity dataset") +
    geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_mw_tcr_contact.pdf")
  ggplot(continuous, aes(x = mw_tcr_contact)) +
    labs(title = "TCR Contact Region Molecular Weight Distribution", subtitle = "IEDB Immunogenicity dataset") +
    geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_hydroph_peptide.pdf")
  ggplot(continuous, aes(x = hydroph_peptide)) +
    labs(title = "Peptide Hydrophobicity Distribution", subtitle = "IEDB Immunogenicity dataset") +
    geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_hydroph_tcr_contact.pdf")  
  ggplot(continuous, aes(x = hydroph_tcr_contact)) +
    labs(title = "TCR Contact Region Hydrophobicity Distribution", subtitle = "IEDB Immunogenicity dataset") +
      geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_charge_peptide.pdf")  
  ggplot(continuous, aes(x = charge_peptide)) +
    labs(title = "Peptide Net Charge Distribution", subtitle = "IEDB Immunogenicity dataset") +
    geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_charge_tcr_contact.pdf")  
  ggplot(continuous, aes(x = charge_tcr_contact)) +
    labs(title = "TCR Contact Region Net Charge Distribution", subtitle = "IEDB Immunogenicity dataset") +
    geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_stability_peptide.pdf")  
  ggplot(continuous, aes(x = stab_peptide)) +
    labs(title = "Peptide Stability Distribution", subtitle = "IEDB Immunogenicity dataset") +
    geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_stability_tcr_contact.pdf")  
  ggplot(continuous, aes(x = stab_tcr_contact)) +
    labs(title = "TCR contact region Stability Distribution", subtitle = "IEDB Immunogenicity dataset") +
      geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_TAP_efficiency.pdf")  
  ggplot(continuous, aes(x = TAP_efficiency)) +
    labs(title = "TAP Transport Efficiency Distribution", subtitle = "IEDB Immunogenicity dataset") +
      geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_proteasomal_cleavage.pdf")  
  ggplot(continuous, aes(x = Proteasomal_Cleavage)) +
    labs(title = "Proteasomal Cleavage Distribution", subtitle = "IEDB Immunogenicity dataset") +
      geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_affinity.pdf")  
  ggplot(continuous, aes(x = specific_allele_mhcflurry_affinity)) +
    labs(title = "Binding Afffinity (nM) Distribution", subtitle = "IEDB Immunogenicity dataset") +
      geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_affinity_percentile.pdf")  
  ggplot(continuous, aes(x = specific_allele_mhcflurry_affinity_percentile)) +
    labs(title = "Binding Afffinity (Percentile) Distribution", subtitle = "IEDB Immunogenicity dataset") +
      geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_processing.pdf")  
  ggplot(continuous, aes(x = specific_allele_mhcflurry_processing_score)) +
    labs(title = "Processing Efficiency (MHCflurry) Distribution", subtitle = "IEDB Immunogenicity dataset") +
      geom_density(alpha = .2, fill = "#FF6666")
  dev.off()
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_feature_distributions/density_presentation.pdf")  
  ggplot(continuous, aes(x = specific_allele_mhcflurry_presentation_score)) +
    labs(title = "Presentation Efficiency (MHCflurry) Distribution", subtitle = "IEDB Immunogenicity dataset") +
      geom_density(alpha = .2, fill = "#FF6666")
  dev.off()

###########################
  
## EXTRACT IMMUNOGENICITY (RESPONSE VARIABLE)
## Pull response variable as immunogenicity_iedb as vector
 

 immunogenicity_iedb <- pull(predictor_iedb_filter, Immunogenicity)
  class(immunogenicity_iedb) ## factor --> change to vector
  immunogenicity_iedb_vector <- as.vector(immunogenicity_iedb)
  setequal(immunogenicity_iedb, immunogenicity_iedb_vector)
  length(immunogenicity_iedb_vector)
  
## EXTRACT PREDICTORS OR COVARIABLES
  ## Pull predictors as predictors_iedb as matrix
  #predictors_iedb <- pull(predictor_iedb_filter, mw_peptide, mw_tcr_contact, hydroph_peptide, hydroph_tcr_contact, charge_peptide, charge_tcr_contact, stab_peptide, stab_tcr_contact, TAP_efficiency, Proteasomal_Cleavage, specific_allele_mhcflurry_affinity, specific_allele_mhcflurry_affinity_percentile, specific_allele_mhcflurry_processing_score, specific_allele_mhcflurry_presentation_score)
  ## pull only works for a single column
  
## MERGE EXTRACTED PREDICTORS INTO MATRIX
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
  #####  GLM PACKAGE
  
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
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_predictor_model/predictor_ROC_curve.pdf")
  ROC_bsc <- plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2, 1.7),
               main = "ROC Curve \n BSC Predictor Vs. Immunogenicity",
                sub = "Generalized Linear Model (glm) - Logistic Regression")
dev.off()
  
  ### Global F test
  summary(logit)
  
  
######## IMMUNOGENICITY PREDICTOR MODEL (glmnet)#######
## ADVANCED MODEL WITH FEATURE SELECTION
  
  ## fit_lasso = glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial")
  ## does not work due to blanks at last row
  ## fixed

  
  #### LASSO PENALIZATION (alpha = 1)
  ## Compute the model with all the variables  
  fit_lasso = glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", alpha = 1)
  fit_lasso
  
  ##Plot the initial model
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_lasso_model/lasso_coefficients_plot.pdf")
  plot(fit_lasso, label = TRUE)
  dev.off()
  plot(fit_lasso, xvar = "dev", label = TRUE)
  
  ## Model summary
  print(fit_lasso)
  
  ## Check coefficients
  coef(fit_lasso,s=0.1) ## do not understand the output
  
  ## Cross-validation
  cvfit_lasso = cv.glmnet(matrix_iedb, immunogenicity_iedb_vector, family = "binomial", standardize = TRUE)
  ##cvfit_lasso = cv.glmnet(matrix_iedb, immunogenicity_iedb_vector, family = "binomial") equal results

  # Plot cross-validated data
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_lasso_model/lasso_crossvalidation.pdf")
  plot(cvfit_lasso)  
  dev.off()
  
  ## Minimum crossvalidated error
  cvfit_lasso$lambda.min
  
  ## lambda.1se, which gives the most regularized model such that error is within one standard error of the minimum
  cvfit_lasso$lambda.1se
  
  ## Coefficients with minimum crossvalidated error
  coef(cvfit_lasso, s = "lambda.min")

  ## Coefficients with 1se
  coef(cvfit_lasso, s = "lambda.1se")

  ## LOGISTIC REGRESSION SPECIFIC INSTRUCTIONS
  
  ## Types of measue
  # Class
  cvfit_lasso_class = cv.glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", type.measure = "class")
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_lasso_model/lasso_class.pdf")
  plot(cvfit_lasso_class)
  dev.off()
  
  # deviance
  cvfit_lasso_deviance = cv.glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", type.measure = "deviance")
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_lasso_model/lasso_deviance.pdf")
  plot(cvfit_lasso_deviance)
  dev.off()

  # mean absolute error  
  cvfit_lasso_mae = cv.glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", type.measure = "mae")
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_lasso_model/lasso_mae.pdf")
  plot(cvfit_lasso_mae)  
  dev.off()
  

  # AUC score
  cvfit_lasso_auc = cv.glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", type.measure = "auc")
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_lasso_model/lasso_auc.pdf")
  plot(cvfit_lasso_auc)  
  dev.off()
  
  
  ## COEFFICIENTS
  cvfit_lasso$lambda.min #0.0003917452
  cvfit_lasso$lambda.1se # 0.004400567
  
  # coef in lambda min
  coef(cvfit_lasso, s = "lambda.min")
  
  # coef in lambda 1se
  coef(cvfit_lasso, s = "lambda.1se")
  
  
  ################
  
  ## ELASTIC NET PENALIZATION (alpha = 0.5)
  
  ## Compute the model with all the variables  
  fit_net = glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", alpha = 0.5)
  
  ##Plot the initial model
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_net_model/net_")
  plot(fit_net, label = TRUE)
  dev.off()
  plot(fit_net, xvar = "dev", label = TRUE)
  
  ## Model summary
  print(fit_net)
  
  ## Check coefficients
  coef(fit_net,s=0.1) ## do not understand the output
  
  ## Cross-validation
  cvfit_net_crossV = cv.glmnet(matrix_iedb, immunogenicity_iedb_vector, family = "binomial", alpha = 0.5, standardize = TRUE)
  ##cvfit_net = cv.glmnet(matrix_iedb, immunogenicity_iedb_vector, family = "binomial") equal results
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_net_model/net_crossV.pdf")
  plot(cvfit_net_crossV)  
  dev.off()
  
  ## Minimum crossvalidated error
  cvfit_net_crossV$lambda.min #0.0003917452
  
  ## lambda.1se, which gives the most regularized model such that error is within one standard error of the minimum
  cvfit_net_crossV$lambda.1se #0.008439894
  
  ## Coefficients with minimum crossvalidated error
  coef(cvfit_net_crossV, s = "lambda.min") 
  
  ## Coefficients with 1se
  coef(cvfit_net_crossV, s = "lambda.1se")
  
  ## LOGISTIC REGRESSION SPECIFIC INSTRUCTIONS
  
  ## Types of measue
  # Class
  cvfit_net_crossV_class = cv.glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", alpha = 0.5, standardize = TRUE, type.measure = "class")
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_net_model/net_crossV_class.pdf")
  plot(cvfit_net_crossV_class)
  dev.off()
  
  # deviance
  cvfit_net_crossV_deviance = cv.glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", alpha = 0.5, standardize = TRUE, type.measure = "deviance")
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_net_model/net_crossV_deviance.pdf")
  plot(cvfit_net_crossV_deviance)  
  dev.off()
  
  # mean absolute error  
  cvfit_net_crossV_mae = cv.glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", alpha = 0.5, standardize = TRUE, type.measure = "mae")
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_net_model/net_crossV_mae.pdf")
  plot(cvfit_net_crossV_mae)  
  dev.off()
  
  # AUC score
  cvfit_net_crossV_auc = cv.glmnet(x = matrix_iedb, y = immunogenicity_iedb_vector, family = "binomial", alpha = 0.5, standardize = TRUE, type.measure = "auc")
  
  pdf("/Users/rocfarriolduran/Desktop/BSC/A_NEOANTIGENS/A_DEVELOPMENT/20200908_immunogenicity_predictor_model/figs_net_model/net_crossV_auc.pdf")
  plot(cvfit_net_crossV_auc)
  dev.off()
  
  ## COEFFICIENTS
  cvfit_net_crossV$lambda.min #0.0134387
  cvfit_net_crossV$lambda.1se #0.07171829
  
  # coef in lambda min
  coef(cvfit_net_crossV, s = "lambda.min")
  coef(cvfit_net_crossV, s = "lambda.1se")  
  
  
  ##### GLOBAL F test, T-test, SUMMARY
  
  ## models fit_lasso // fit_net
  
  summary(cvfit_lasso)
  summary(cvfit_net_crossV)

  ## PERFORMANCE LASSO
  ## PERFORMANCE ELASTIC NET
