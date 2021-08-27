if(!require(e1071)) install.packages('e1071', dependencies = TRUE)
if(!require(caret)) install.packages('caret', repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages('dplyr', repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages('tidyverse', repos = "http://cran.us.r-project.org")
if(!require(readxl)) install.packages('readxl', repos = "http://cran.us.r-project.org")
if(!require(xlsx)) install.packages('xlsx', repos = "http://cran.us.r-project.org")
if(!require(glmnet)) install.packages('glmnet', repos = "http://cran.us.r-project.org")
if(!require(dummies)) install.packages('dummies', repos = "http://cran.us.r-project.org")
if(!require(mRMRe)) install.packages('mRMRe', repos = "http://cran.us.r-project.org")
if(!require(geometry)) install.packages('geometry', repos = "http://cran.us.r-project.org")
if(!require(survival)) install.packages('survival', repos = "http://cran.us.r-project.org")
if(!require(survminer)) install.packages('survminer', repos = "http://cran.us.r-project.org")
if(!require(ggfortify)) install.packages('ggfortify', repos = "http://cran.us.r-project.org")
if(!require(ggpubr)) install.packages('ggpubr', repos = "http://cran.us.r-project.org")
if(!require(MatchIt)) install.packages('MatchIt', repos = "http://cran.us.r-project.org")
if(!require(cobalt)) install.packages('cobalt', repos = "http://cran.us.r-project.org")
# if(!require(optmatch)) install.packages('optmatch', repos = "http://cran.us.r-project.org")
# if(!require(CBPS)) install.packages('CBPS', repos = "http://cran.us.r-project.org")


if(!require(writexl)) install.packages('writexl', repos = "http://cran.us.r-project.org")

set.seed(123)

num_mrmr_features = 100

#================================================
# Propensity score 1-1 matching for LUAD and LUSC#
#================================================

local_path_luad_lusc = '/Users/rinading/Desktop/UCLA/WQE/data/'
# Read in the files
data_for_matching_dev = read_excel(paste0(local_path_luad_lusc, "for_propensity_score_matching_dev.xlsx"), col_names = TRUE)

# Check initial imbalanced data before matching
unmatched_results_dev = matchit(subtype ~ stage + gender, data = data_for_matching_dev,
                  method = NULL, distance = "glm")
# 1-1 nearest neighbor logistic regression matching
matched_results_dev = matchit(subtype ~ stage + gender, data = data_for_matching_dev,
                                method = "nearest", distance = "glm")

plot(matched_results_dev, type = "jitter", interactive = FALSE)
love.plot(bal.tab(matched_results_dev), stat = "mean.diffs", threshold = .1, 
          var.order = "unadjusted", drop.distance = TRUE)

matched_indices_dev = matched_results_dev$match.matrix

luad_radiomic_dev = as.data.frame(read_excel(paste0(local_path_luad_lusc, "luad_radiomic_dev.xlsx"), col_names = TRUE))
luad_clinical_dev = as.data.frame(read_excel(paste0(local_path_luad_lusc, "luad_clinical_dev.xlsx"), col_names = TRUE))
luad_genomic_dev = as.data.frame(read_excel(paste0(local_path_luad_lusc, "luad_genomic_dev.xlsx"), col_names = TRUE))

luad_genomic_dev$...1 = NULL
luad_radiomic_test = as.data.frame(read_excel(paste0(local_path_luad_lusc, "luad_radiomic_test.xlsx"), col_names = TRUE))
luad_clinical_test = as.data.frame(read_excel(paste0(local_path_luad_lusc, "luad_clinical_test.xlsx"), col_names = TRUE))
luad_genomic_test = as.data.frame(read_excel(paste0(local_path_luad_lusc, "luad_genomic_test.xlsx"), col_names = TRUE))
luad_genomic_test$...1 = NULL

lusc_radiomic_dev = as.data.frame(read_excel(paste0(local_path_luad_lusc, "lusc_radiomic_dev.xlsx"), col_names = TRUE))
lusc_clinical_dev = as.data.frame(read_excel(paste0(local_path_luad_lusc, "lusc_clinical_dev.xlsx"), col_names = TRUE))
lusc_genomic_dev = as.data.frame(read_excel(paste0(local_path_luad_lusc, "lusc_genomic_dev.xlsx"), col_names = TRUE))
lusc_genomic_dev$...1 = NULL

lusc_radiomic_test = as.data.frame(read_excel(paste0(local_path_luad_lusc, "lusc_radiomic_test.xlsx"), col_names = TRUE))
lusc_clinical_test = as.data.frame(read_excel(paste0(local_path_luad_lusc, "lusc_clinical_test.xlsx"), col_names = TRUE))
lusc_genomic_test = as.data.frame(read_excel(paste0(local_path_luad_lusc, "lusc_genomic_test.xlsx"), col_names = TRUE))
lusc_genomic_test$...1 = NULL


# Only select the matched LUAD cases
luad_radiomic_dev_matched = luad_radiomic_dev[matched_indices_dev, ]
luad_clinical_dev_matched = luad_clinical_dev[matched_indices_dev, ]
luad_genomic_dev_matched = luad_genomic_dev[, as.integer(matched_indices_dev)]

df_radiomic_dev = luad_radiomic_dev_matched
df_clinical_dev = luad_clinical_dev_matched
df_radiomic_test = luad_radiomic_test
df_clinical_test = luad_clinical_test

# Simplify the radiomic feature names
for (col in 1:ncol(df_radiomic_dev)){
  colnames(df_radiomic_dev)[col] = sub("imaging.radiomics.", "", colnames(df_radiomic_dev)[col])
}

for (col in 1:ncol(df_radiomic_test)){
  colnames(df_radiomic_test)[col] = sub("imaging.radiomics.", "", colnames(df_radiomic_test)[col])
}

luad_genomic_all = cbind(luad_genomic_dev_matched, luad_genomic_test)
lusc_genomic_all = cbind(lusc_genomic_dev, lusc_genomic_test)
#write_xlsx(luad_genomic_all_matched, '/Users/rinading/Desktop/UCLA/WQE/data/GSEA/matched_luad_genomic_all.xlsx')
#write_xlsx(lusc_genomic_all, '/Users/rinading/Desktop/UCLA/WQE/data/GSEA/lusc_genomic_all.xlsx')
#write_xlsx(luad_clinical_dev_matched, '/Users/rinading/Desktop/UCLA/WQE/data/matched_luad_clinical_dev.xlsx')
#write_xlsx(luad_clinical_test_matched, '/Users/rinading/Desktop/UCLA/WQE/data/matched_luad_clinical_test.xlsx')

# Function for plotting the KM curves and calculate p-value from log-rank test
KM_curves_from_risk = function(risk_scores, risk_threshold, time, status, coxph_data){
  risk_groups = as.factor(case_when(
    risk_scores >= risk_threshold ~ "High-risk",
    risk_scores < risk_threshold ~ "Low-risk",
  ))
  
  model = coxph(Surv(time, status) ~ risk_scores, data = coxph_data)
  summary(model)
  
  # Plot the KM curves for the 2 risk groups
  coxph_data_2_groups = as.data.frame(cbind(time, status, risk_groups))
  fit = survfit(Surv(time, status) ~ risk_groups, data = coxph_data_2_groups)
  ggsurvplot(
    fit,  # survfit object
    data = coxph_data_2_groups, # data used to fit survival curves
    risk.table = TRUE,  # show risk table.
    pval = TRUE,      # show p-value of log-rank test.
    conf.int = TRUE,         # show confidence intervals for
    # point estimates of survival curves.
    xlab = "Time in months",   # customize X axis label.
    break.time.by = 20,     # break X axis in time intervals by xxx.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE, # show bars instead of names in text annotations
    # in legend of risk table
    #risk.table.col = "strata",# Risk table color by groups
    legend.labs = c("High-risk", "Low-risk")    # Change legend labels
  )
  return (model)
}
#================================================
# Pre-selection of features #
#================================================
# Radiomic data
radiomic_with_events = cbind(df_radiomic_dev, df_clinical_dev$surv.clinical.Event)
radiomic_mrmr_type = mRMR.data(data = radiomic_with_events)
radiomic_selected = mRMR.classic(data = radiomic_mrmr_type, target_indices = ncol(radiomic_with_events), feature_count = num_mrmr_features)
radiomic_selected_indices = unlist(radiomic_selected@filters)
df_radiomic_dev_mrmr_selected = df_radiomic_dev[, radiomic_selected_indices]
# Use the same features as in the development set for the testing set.
df_radiomic_test_mrmr_selected = df_radiomic_test[, radiomic_selected_indices]

#================================================
# Training the cox regression model #
#================================================

Y_dev = as.matrix(cbind(df_clinical_dev$surv.clinical.Time.Months, df_clinical_dev$surv.clinical.Event))
colnames(Y_dev) = c('time', 'status')
X_dev = as.matrix(df_radiomic_dev_mrmr_selected)
# Fit the Cox regression model
cvfit = cv.glmnet(X_dev, Y_dev, family = 'cox', type.measure = 'C', alpha = 1, nfolds = 5)
plot(cvfit)
# See what features are selected and their weights
beta = coef(cvfit, s = "lambda.min")
#beta = coef(cvfit)
lasso_selected_feature_indices = beta@i
lasso_selected_feature_beta = beta@x

# Get the risk scores for each patient
risk_scores_dev = predict(cvfit, newx = X_dev, s = "lambda.min")
c_index_dev = apply(risk_scores_dev, 2, Cindex, y = Y_dev) # Calculate concordance index using the risk scores
# Use the median risk score as the risk threshold
risk_threshold = median(risk_scores_dev)
risk_groups_dev = as.factor(case_when(
  risk_scores_dev >= risk_threshold ~ "High-risk",
  risk_scores_dev < risk_threshold ~ "Low-risk",
))

time_dev = df_clinical_dev$surv.clinical.Time.Months
status_dev = df_clinical_dev$surv.clinical.Event
coxph_data_dev = as.data.frame(cbind(time_dev, status_dev, risk_scores_dev))

KM_curves_from_risk(risk_scores_dev, risk_threshold, time_dev, status_dev, coxph_data_dev)
#================================================
# Testing the cox regression model #
#================================================
Y_test = as.matrix(cbind(df_clinical_test$surv.clinical.Time.Months, df_clinical_test$surv.clinical.Event))
colnames(Y_test) = c('time', 'status')
X_test = as.matrix(df_radiomic_test_mrmr_selected)
risk_scores_test = predict(cvfit, newx = X_test, s = "lambda.min")
c_index_test = apply(risk_scores_test, 2, Cindex, y=Y_test)

risk_groups_test = as.factor(case_when(
  risk_scores_test >= risk_threshold ~ "High-risk",
  risk_scores_test < risk_threshold ~ "Low-risk",
))

time_test = df_clinical_test$surv.clinical.Time.Months
status_test = df_clinical_test$surv.clinical.Event
coxph_data_test = as.data.frame(cbind(time_test, status_test, risk_scores_test))

test = KM_curves_from_risk(risk_scores_test, risk_threshold, time_test, status_test, coxph_data_test)

# Save predicted risk scores for dev and test sets
risk_scores_all = as.data.frame(rbind(risk_scores_dev, risk_scores_test))
write_xlsx(as.data.frame(risk_scores_dev), '/Users/rinading/Desktop/UCLA/WQE/data/GSEA/luad_risk_scores_all.xlsx')
#================================================
# Visualize top feature values #
#================================================
# Sort features indicated by LASSO feature weight absolute value
sorted_features_indices = sort.int(abs(lasso_selected_feature_beta), decreasing = TRUE, index.return = TRUE)$ix
sorted_features_weights = round(sort.int(abs(lasso_selected_feature_beta), decreasing = TRUE, index.return = TRUE)$x, digits = 2)
selected_features = colnames(X_dev[, lasso_selected_feature_indices[sorted_features_indices]])
selected_features_and_weights = cbind(selected_features, sorted_features_weights)

# Save the selected features to a matrix 
selected_features_dev = X_dev[, lasso_selected_feature_indices[sorted_features_indices]]
selected_features_test = X_test[, lasso_selected_feature_indices[sorted_features_indices]]
selected_features_all = as.data.frame(rbind(selected_features_dev, selected_features_test))
write_xlsx(selected_features_all, '/Users/rinading/Desktop/UCLA/WQE/data/enrichment_radiomic_corr/selected_radiomic_luad.xlsx')

# Compare top features in train vs test cohorts
train_labels = rep("development", nrow(X_dev))
test_labels = rep("test", nrow(X_test))
top_feature = 5
train_labels_top_feature = cbind(X_dev[, lasso_selected_feature_indices[sorted_features_indices[top_feature]]], train_labels)
test_labels_top_feature = cbind(X_test[, lasso_selected_feature_indices[sorted_features_indices[top_feature]]], test_labels)

train_test_labels_top_feature = rbind(train_labels_top_feature, test_labels_top_feature)
df_train_test_labels_top_feature = as.data.frame(train_test_labels_top_feature)
colnames(df_train_test_labels_top_feature) = c("top_feature_value", "cohort")
df_train_test_labels_top_feature = transform(df_train_test_labels_top_feature, top_feature_value = as.double(top_feature_value))

ggboxplot(df_train_test_labels_top_feature, x = "cohort", y = "top_feature_value",
               color = "cohort", add = "jitter") + stat_compare_means(method = "wilcox.test")

# Check the Wilxocon p value for all selected features
wilcoxon_results = c()
for (i in 1:length(sorted_features_indices)){
  x = X_dev[, lasso_selected_feature_indices[sorted_features_indices[i]]]
  y = X_test[, lasso_selected_feature_indices[sorted_features_indices[i]]]
  wilcoxon_results = c(wilcoxon_results, wilcox.test(x, y, p.adjust.method = "fdr")$p.value)
}
