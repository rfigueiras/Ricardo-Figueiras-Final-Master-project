
# Testing proteomic signature


# Data loading
vsn_knn_combat<-read.csv("vsn_knn_combat.csv", row.names = 1, check.names = FALSE)       
metadata<-read.csv("metadata_merged.csv", row.names = 1, check.names = FALSE)
selected_features<-read.csv("rfe_selected_features_vsn_knn_combat.csv", row.names = 1, check.names = FALSE)


# library
library(caret)
library(glmnet)
library(caretEnsemble)
library(e1071)
library(doParallel)
library(caret)
library(glmnet)
library(caretEnsemble)
library(e1071)
library(doParallel)
library(ggplot2)
library(reshape2)
library(pROC)
library(dbplyr)
library(tidyverse)
library(caret)
library(caretEnsemble)
library(pROC)
library(ggplot2)


# Assigning training data and testing data
train_data <- vsn_knn_combat[, metadata$MS_batch %in% c("Batch1", "Batch2")]
test_data <-vsn_knn_combat[, metadata$MS_batch == "Batch3"]
  
train_labels <- metadata$Group[metadata$MS_batch %in% c("Batch1", "Batch2")]
test_labels <- metadata$Group[metadata$MS_batch == "Batch3"]
  
train_labels <- factor(train_labels)
test_labels <- factor(test_labels)
  
train_df <- as.data.frame(t(train_data))
train_df$Group <- as.factor(train_labels)
test_df <- as.data.frame(t(test_data))
test_df$Group <- as.factor(test_labels)
  



# Subset the training and test data based on selected features, selected
# previously with RFE
train_df_selected <- train_df[, c(selected_features$x, "Group")]
test_df_selected <- test_df[, c(selected_features$x, "Group")]


dim(train_df_selected)
dim(test_df_selected)

#set HCC as reference level
train_df_selected$Group <- relevel(train_df_selected$Group, ref = "HCC")
test_df_selected$Group <- relevel(test_df_selected$Group, ref = "HCC")


table(train_df_selected$Group)
table(test_df_selected$Group)




# Train models using selected features by RFE
set.seed(12355)
# Define train control with repeated cross-validation
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 100,
                              classProbs = TRUE,
                              savePredictions = TRUE,
                              summaryFunction = defaultSummary)


# Set up parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


models_pipeline3 <- caretList(
  Group ~ ., data = train_df_selected,
  trControl = train_control,
  methodList = c("rf", "svmRadial", "nnet", "glmnet", "knn", "lda")# algorithms to use
)

# Stop parallel processing
stopCluster(cl)


# Collect resampling results
models_pipeline3_results <- resamples(models_pipeline3)

# Print the summary of results
summary(models_pipeline3_results)


# Prediction of batch3 samples using random forest model (which performed the best
#during cross-validation on the trainning data)


rf_model<-models_pipeline3$rf



# remove the target column from the test data for predictions
test_data_for_predictions <- test_df_selected[, -ncol(test_df_selected)]

# predictions on the test set 
pred_rf <- predict(rf_model, newdata = test_data_for_predictions, type = "prob")
pred_rf_class <- predict(rf_model, newdata = test_data_for_predictions)


# Calculate the ROC curve
roc_curve_rf <- roc(test_df_selected$Group, pred_rf[, 2], levels = levels(test_df_selected$Group), ci = TRUE)
auc_value_rf <- auc(roc_curve_rf)
ci_value_rf <- ci.auc(roc_curve_rf)


roc_data <- data.frame(
  fpr = 1-roc_curve_rf$specificities, #get False positive rate
  tpr = roc_curve_rf$sensitivities #get True positive rate
)

auc_value_rf <- auc(roc_curve_rf)
ci_values <- ci(roc_curve_rf)#cofidence interval



# confusion matrix metrics
conf_matrix_rf <- confusionMatrix(pred_rf_class, test_df_selected$Group)
print(conf_matrix_rf)


# function to calculate peformance metrics from confusion matrix:

# initialize data frame
results <- data.frame(
  Model = character(),
  Accuracy = numeric(),
  P_Value = numeric(),
  AUC = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  F1_Score = numeric(),
  stringsAsFactors = FALSE
)


calculate_metrics <- function(pred_class, true_labels, roc_obj) {
  auc_value <- auc(roc_obj)
  ci_value <- ci.auc(roc_obj)
  conf_matrix <- confusionMatrix(pred_class, true_labels)
  
  accuracy <- conf_matrix$overall['Accuracy']
  sensitivity <- conf_matrix$byClass['Sensitivity']
  specificity <- conf_matrix$byClass['Specificity']
  f1_score <- conf_matrix$byClass['F1']
  p_value_Acc_NIR<- conf_matrix$overall['AccuracyPValue']
  
  return(list(Accuracy = accuracy, P_Value_Acc_NIR = p_value_Acc_NIR, AUC = auc_value, 
              Sensitivity = sensitivity, Specificity = specificity, F1_Score = f1_score))
}

# apply function
metrics_testing_rf <- calculate_metrics(pred_rf_class, test_df_selected$Group, roc_curve_rf)

write.csv(as.data.frame(metrics_testing_rf), "performance_metrics_testing_MS.csv")
write.csv(conf_matrix_rf$table, "confusion_matrix_table_testing_MS.csv")


# PLOT ROC CURVE

p <- ggroc(roc_curve_rf, color = "indianred", size = 2) +
  scale_x_reverse() +  # Reverse the x-axis
  ggtitle("Mass-spectometry testing \nrandom forest model ROC Curve") +
  xlab("Specificity (1 - FPR)") +  
  ylab("Sensitivity (TPR)") +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black") +  # Adjust the diagonal line to match the reversed scale
  theme_bw() +
  annotate("text", x = 0.3, y = 0.25, 
           label = paste0("AUC: ", round(auc_value_rf, 2), "\n95% CI: (", 
                          round(ci_values[1], 2), "-", round(ci_values[3], 2), ")"), 
           size = 14, color = "black") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30), 
        plot.title = element_text(size = 32),panel.grid.major.y = element_blank(), panel.grid.major = element_line(size = 1),  
        panel.grid.minor = element_line(size = 1), 
        panel.border = element_rect(size = 1)       
  ) 

# Display the plot
print(p)



ggsave("testing MS.png", plot = p, width = 12, height = 12, dpi=600)



# Get variable importance ranking of the random-forest model used

importance_rf <- varImp(rf_model, scale = FALSE)
print(importance_rf$importance)

write(as.matrix(importance_rf$importance), "variable_importance_rf.csv")
plot(importance_rf)



write.csv(roc_data, "ROC curve validation batch3.csv")


protein_name_map<-read.csv("protein_names.csv")



importance_rf_anonimised <- importance_rf
name_mapping <- setNames(protein_name_map$Name, protein_name_map$Real_name)
# Replace the row names with their corresponding anonymous name
for (real_name in rownames(importance_rf_anonimised$importance)) {
  if (real_name %in% names(name_mapping)) {
    rownames(importance_rf_anonimised$importance)[rownames(importance_rf_anonimised$importance) == real_name] <- name_mapping[real_name]
  }
}


importance_rf_anonimised_df <- data.frame(Protein = rownames(importance_rf_anonimised$importance), 
                                       Importance = importance_rf_anonimised$importance[1])




# PLOT VARIABLE IMPORTANCE RANKINGS

library(ggplot2)

feat_imp<-ggplot(importance_rf_anonimised_df, aes(x = reorder(Protein, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "dodgerblue4") +
  coord_flip() +  # Flip the coordinates for better readability
  labs(title = "Feature Importance \ncross-validation with selected features",
       x = "Proteins",
       y = "Importance") +
  theme_bw() +  # Apply a theme with a white background
  theme(
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32, face="bold"),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.grid.major.y = element_blank(), panel.grid.major = element_line(size = 1),  
    panel.grid.minor = element_line(size = 1), 
    panel.border = element_rect(size = 1)       
  ) 


ggsave("Variable importance RFE pipeline3_v2.png", plot = feat_imp, width = 12, height = 18, dpi=600)


