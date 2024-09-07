
#libraries
library(glmnet)
library(caret)
library(doParallel)
library(pROC)
library(caretEnsemble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(gridExtra)
library(grid)


# loading data
metadata_RNA<-read.csv("metadata_filtered_RNA.csv", row.names = 1, check.names = FALSE)
metadata_protein<-read.csv( "metadata_filtered_proteins.csv", row.names = 1, check.names = FALSE)

data_miRNA<-read.csv("Seurat_counts_filt_miRNA.csv", row.names = 1, check.names = FALSE)
data_proteins<-read.csv("vsn_knn_combat.csv", row.names = 1, check.names = FALSE)

AFP_data<-read.csv("AFP.csv", row.names = 1, check.names = FALSE)

keep_miRNAs <- rownames(data_miRNA)
keep_proteins <- rownames(data_proteins)


data_miRNA_selected<-data_miRNA[keep_miRNAs,]
data_proteins_selected<-data_proteins[keep_proteins,]
metadata_protein_filt<- metadata_protein[metadata_protein$MS_batch=="Batch1"|metadata_protein$MS_batch=="Batch2",]
data_proteins_selected_filt<- data_proteins_selected[,colnames(data_proteins_selected) %in% metadata_protein_filt$Sample_ID]
colnames(data_miRNA_selected) <- gsub("MR_([0-9]+)A?", "\\1", colnames(data_miRNA_selected))
colnames(data_miRNA_selected) <- gsub("A", "", colnames(data_miRNA_selected))


#intersect for common samples on the different omics
common_samples <- intersect(intersect(colnames(data_proteins_selected_filt), colnames(data_miRNA_selected)), metadata_protein_filt$Sample_ID)

data_proteins_matched <- data_proteins_selected_filt[,common_samples]
data_miRNA_matched <- data_miRNA_selected[,common_samples]
metadata_matched <- metadata_protein_filt[common_samples, ]

data_AFP_matched<- AFP_data[, colnames(AFP_data) %in% metadata_matched$Sample_ID]
data_AFP_matched <- as.data.frame(lapply(data_AFP_matched, as.numeric), check.names = FALSE)



dim(data_proteins_matched)
dim(data_miRNA_matched)
dim(metadata_matched)
dim(data_AFP_matched)
rownames(data_AFP_matched)<-"AFP"




#combine omics data
combined_data <- rbind(data_proteins_matched, data_miRNA_matched, data_AFP_matched)
combined_data_t <- t(combined_data)
combined_data_df <- as.data.frame(combined_data_t)


# Outcome variable
group<- as.factor(metadata_matched$Group)
combined_data_df$Group <- group



#------------------------ random splits function--------------------------------------------------------------------

# Function for random splits
split_data <- function(train_df, splits_number) {
  training_datasets <- list()
  testing_datasets <- list()
  
  for (i in 1:splits_number) {
    set.seed(1234 + i)  #use different seeds for each split
    train_df_HCC <- train_df[train_df$Group == "HCC", ]
    train_df_Cirrhosis <- train_df[train_df$Group == "Cirrhosis", ]
    # 75% of HCC samples for training
    indices_HCC <- sample(1:nrow(train_df_HCC), size = 0.75 * nrow(train_df_HCC))
    train_data_HCC <- train_df_HCC[indices_HCC, ]
    # remaining HCC samples for testing
    test_data_HCC <- train_df_HCC[-indices_HCC, ]
    # 75% of cirrhosis samples for training
    indices_Cirrhosis <- sample(1:nrow(train_df_Cirrhosis), size = 0.75 * nrow(train_df_Cirrhosis))
    train_data_Cirrhosis <- train_df_Cirrhosis[indices_Cirrhosis, ]
    # remaining cirrhosis samples for testing
    test_data_Cirrhosis <- train_df_Cirrhosis[-indices_Cirrhosis, ]
    
    train_data_ns <- rbind(train_data_HCC, train_data_Cirrhosis)
    test_data_ns <- rbind(test_data_HCC, test_data_Cirrhosis)
    # scale training data
    scaler <- preProcess(train_data_ns[, -ncol(train_data_ns)], method = c("center", "scale"))
    train_data <- predict(scaler, train_data_ns[, -ncol(train_data_ns)])
    # apply scaling parameter for testing data
    test_data <- predict(scaler, test_data_ns[, -ncol(test_data_ns)])
    
    train_data <- as.data.frame(train_data)
    test_data <- as.data.frame(test_data)
    # add group labels
    train_data$Group <- train_data_ns$Group
    test_data$Group <- test_data_ns$Group
    # set HCC as reference level
    train_data$Group <- relevel(train_data$Group, ref = "HCC")
    test_data$Group <- relevel(test_data$Group, ref = "HCC")
    
    training_datasets[[paste0("train_data_", i)]] <- train_data
    testing_datasets[[paste0("test_data_", i)]] <- test_data
  }
  
  return(list(training_datasets = training_datasets, testing_datasets = testing_datasets))
}

#------------------------ random splits function for early-stage testing----------------


# Function for random splits with early-HCC on the testing subset

split_data_early <- function(train_df, splits_number) {
  training_datasets <- list()
  testing_datasets <- list()
  
  for (i in 1:splits_number) {
    set.seed(123 + i)  # use different seeds for each split
    
    # split the dataset into HCC and Cirrhosis
    hcc <- train_df[train_df$Group == "HCC", ]
    cirrhosis <- train_df[train_df$Group == "Cirrhosis", ]
    
    # filter early stage HCC samples based for testing set
    hcc_bclc_0A <- hcc[metadata_matched_HCC$`BCLC stage` %in% c("0", "A"), ]
    hcc_others <- hcc[!(metadata_matched_HCC$`BCLC stage` %in% c("0", "A")), ]
    
    # select early stage 12 HCC samples for test set
    test_hcc_0A_indices <- sample(1:nrow(hcc_bclc_0A), size = 12)
    
    test_hcc <- hcc_bclc_0A[test_hcc_0A_indices, ]
    
    #remaining HCC samples for training set
    train_hcc <- rbind(hcc_others, hcc_bclc_0A[-test_hcc_0A_indices, ])
    
    #split Cirrhosis samples into 75% training and 25% testing sets
    indices_cirrhosis <- sample(1:nrow(cirrhosis), size = 0.75 * nrow(cirrhosis))
    train_cirrhosis <- cirrhosis[indices_cirrhosis, ]
    test_cirrhosis <- cirrhosis[-indices_cirrhosis, ]
    
    # combine the training sets
    train_data_ns <- rbind(train_hcc, train_cirrhosis)
    
    # combine the testing sets
    test_data_ns <- rbind(test_hcc, test_cirrhosis)
    print(dim(test_hcc))
    print(dim(test_data_ns))
    print(table(test_data_ns$Group))
    # apply scaling
    scaler <- preProcess(train_data_ns[, -ncol(train_data_ns)], method = c("center", "scale"))
    train_data <- predict(scaler, train_data_ns[, -ncol(train_data_ns)])
    test_data <- predict(scaler, test_data_ns[, -ncol(test_data_ns)])
    
    train_data <- as.data.frame(train_data)
    test_data <- as.data.frame(test_data)
    
    train_data$Group <- train_data_ns$Group
    test_data$Group <- test_data_ns$Group
    
    # Relevel HCC as the reference levl
    train_data$Group <- relevel(train_data$Group, ref = "HCC")
    test_data$Group <- relevel(test_data$Group, ref = "HCC")
    print(table(train_data$Group))
    print(table(test_data$Group))
    
    training_datasets[[paste0("train_data_", i)]] <- train_data
    testing_datasets[[paste0("test_data_", i)]] <- test_data
  }
  
  return(list(training_datasets = training_datasets, testing_datasets = testing_datasets))
}

#--------------------------------- Training and testing---------------------------

# use the function to obtain 10 random splits
datasets <- split_data(combined_data_df, 10)

training_datasets<-datasets$training_datasets
testing_datasets<-datasets$testing_datasets


training_datasets$train_data_1


# RFE control with repeated cross-validation
rfe_control_rf <- caret::rfeControl(functions = rfFuncs, method = "repeatedcv", number = 10, repeats = 100, returnResamp = "final", allowParallel = TRUE, verbose = TRUE)



cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


# Recursive feature elimination on the 10 splits 
set.seed(11034)
rfe_results_1 <- caret::rfe(training_datasets$train_data_1[, -which(colnames(training_datasets$train_data_1) == "Group")], training_datasets$train_data_1$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_2 <- caret::rfe(training_datasets$train_data_2[, -which(colnames(training_datasets$train_data_2) == "Group")], training_datasets$train_data_2$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_3 <- caret::rfe(training_datasets$train_data_3[, -which(colnames(training_datasets$train_data_3) == "Group")], training_datasets$train_data_3$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_4 <- caret::rfe(training_datasets$train_data_4[, -which(colnames(training_datasets$train_data_4) == "Group")], training_datasets$train_data_4$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_5 <- caret::rfe(training_datasets$train_data_5[, -which(colnames(training_datasets$train_data_5) == "Group")], training_datasets$train_data_5$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_6 <- caret::rfe(training_datasets$train_data_6[, -which(colnames(training_datasets$train_data_6) == "Group")], training_datasets$train_data_6$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_7 <- caret::rfe(training_datasets$train_data_7[, -which(colnames(training_datasets$train_data_7) == "Group")], training_datasets$train_data_7$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_8 <- caret::rfe(training_datasets$train_data_8[, -which(colnames(training_datasets$train_data_8) == "Group")], training_datasets$train_data_8$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_9 <- caret::rfe(training_datasets$train_data_9[, -which(colnames(training_datasets$train_data_9) == "Group")], training_datasets$train_data_9$Group, sizes = c(2:30), rfeControl = rfe_control_rf)

rfe_results_10 <- caret::rfe(training_datasets$train_data_10[, -which(colnames(training_datasets$train_data_10) == "Group")], training_datasets$train_data_10$Group, sizes = c(2:30), rfeControl = rfe_control_rf)
stopCluster(cl)

#obtain and subset selected features
selected_features_rf_1 <- predictors(rfe_results_1)
selected_features_rf_2 <- predictors(rfe_results_2)
selected_features_rf_3 <- predictors(rfe_results_3)
selected_features_rf_4 <- predictors(rfe_results_4)
selected_features_rf_5 <- predictors(rfe_results_5)
selected_features_rf_6 <- predictors(rfe_results_6)
selected_features_rf_7 <- predictors(rfe_results_7)
selected_features_rf_8 <- predictors(rfe_results_8)
selected_features_rf_9 <- predictors(rfe_results_9)
selected_features_rf_10 <- predictors(rfe_results_10)



train_data_1_selected<-as.data.frame(training_datasets$train_data_1[,c(selected_features_rf_1, "Group")])
train_data_2_selected<-as.data.frame(training_datasets$train_data_2[,c(selected_features_rf_2, "Group")])
train_data_3_selected<-as.data.frame(training_datasets$train_data_3[,c(selected_features_rf_3, "Group")])
train_data_4_selected<-as.data.frame(training_datasets$train_data_4[,c(selected_features_rf_4, "Group")])
train_data_5_selected<-as.data.frame(training_datasets$train_data_5[,c(selected_features_rf_5, "Group")])
train_data_6_selected<-as.data.frame(training_datasets$train_data_6[,c(selected_features_rf_6, "Group")])
train_data_7_selected<-as.data.frame(training_datasets$train_data_7[,c(selected_features_rf_7, "Group")])
train_data_8_selected<-as.data.frame(training_datasets$train_data_8[,c(selected_features_rf_8, "Group")])
train_data_9_selected<-as.data.frame(training_datasets$train_data_9[,c(selected_features_rf_9, "Group")])
train_data_10_selected<-as.data.frame(training_datasets$train_data_10[,c(selected_features_rf_10, "Group")])




# use the selected features to train the model with the training data

# Set up parallel processing
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)


# Define the training control
train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 100, savePredictions = "final", classProbs = TRUE)


models_featsel_1 <- caretList(
  Group ~ ., data = train_data_1_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)


models_featsel_2 <- caretList(
  Group ~ ., data = train_data_2_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)



models_featsel_3 <- caretList(
  Group ~ ., data = train_data_3_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)



models_featsel_4 <- caretList(
  Group ~ ., data = train_data_4_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)

models_featsel_5 <- caretList(
  Group ~ ., data = train_data_5_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)

models_featsel_6 <- caretList(
  Group ~ ., data = train_data_6_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)

models_featsel_7 <- caretList(
  Group ~ ., data = train_data_7_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)


models_featsel_8 <- caretList(
  Group ~ ., data = train_data_8_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)

models_featsel_9 <- caretList(
  Group ~ ., data = train_data_9_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)

models_featsel_10 <- caretList(
  Group ~ ., data = train_data_10_selected,
  trControl = train_control,
  methodList = c("rf", "glmnet")
)

# Stop parallel processing
stopCluster(cl)


# Get resampling results
results_models_featsel_1 <- resamples(models_featsel_1)
# Print the summary of results
summary(results_models_featsel_1)




# subset selected features on the testing data 
test_data_1_selected<-as.data.frame(testing_datasets$test_data_1[,c(selected_features_rf_1, "Group")])
test_data_2_selected<-as.data.frame(testing_datasets$test_data_2[,c(selected_features_rf_2, "Group")])
test_data_3_selected<-as.data.frame(testing_datasets$test_data_3[,c(selected_features_rf_3, "Group")])
test_data_4_selected<-as.data.frame(testing_datasets$test_data_4[,c(selected_features_rf_4, "Group")])
test_data_5_selected<-as.data.frame(testing_datasets$test_data_5[,c(selected_features_rf_5, "Group")])
test_data_6_selected<-as.data.frame(testing_datasets$test_data_4[,c(selected_features_rf_6, "Group")])
test_data_7_selected<-as.data.frame(testing_datasets$test_data_4[,c(selected_features_rf_7, "Group")])
test_data_8_selected<-as.data.frame(testing_datasets$test_data_4[,c(selected_features_rf_8, "Group")])
test_data_9_selected<-as.data.frame(testing_datasets$test_data_4[,c(selected_features_rf_9, "Group")])
test_data_10_selected<-as.data.frame(testing_datasets$test_data_10[,c(selected_features_rf_10, "Group")])


#Predict 10 dataframes of test data using the respective trained model

# Predict rf_1 
pred_rf_1 <- predict(models_featsel_1$rf, newdata = test_data_1_selected[, -ncol(test_data_1_selected)], type = "prob")
pred_rf_class_1 <- predict(models_featsel_1$rf, newdata = test_data_1_selected[, -ncol(test_data_1_selected)])

# Calculate the ROC curve
roc_curve_rf_1 <- roc(test_data_1_selected$Group, pred_rf_1[, 2], levels = levels(test_data_1_selected$Group), ci = TRUE)
auc_value_rf_1 <- auc(roc_curve_rf_1)
ci_value_rf_1<- ci.auc(roc_curve_rf_1)

# Plot the ROC curve
plot(roc_curve_rf_1, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_1, 3), ")"))
ci(roc_curve_rf_1)


print(paste("AUC:", round(auc_value_rf_1, 3)))
print(paste("910% CI:", round(ci_value_rf_1[1], 3), "-", round(ci_value_rf_1[3], 3)))

# confusion matrix
conf_matrix_rf_1 <- confusionMatrix(pred_rf_class_1 , testing_datasets$test_data_1$Group)
print(conf_matrix_rf_1)




rf_model2 <- models_featsel_2$rf

# Predict rf_2 
pred_rf_2 <- predict(models_featsel_2$rf, newdata = test_data_2_selected[, -ncol(test_data_2_selected)], type = "prob")
pred_rf_class_2 <- predict(models_featsel_2$rf, newdata = test_data_2_selected[, -ncol(test_data_2_selected)])

# Calculate the ROC curve
roc_curve_rf_2 <- roc(test_data_2_selected$Group, pred_rf_2[, 2], levels = levels(test_data_2_selected$Group), ci = TRUE)
auc_value_rf_2 <- auc(roc_curve_rf_2)
ci_value_rf_2<- ci.auc(roc_curve_rf_2)

# Plot the ROC curve
plot(roc_curve_rf_2, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_2, 3), ")"))
ci(roc_curve_rf_2)


print(paste("AUC:", round(auc_value_rf_2, 3)))
print(paste("910% CI:", round(ci_value_rf_2[1], 3), "-", round(ci_value_rf_2[3], 3)))

# confusion matrix
conf_matrix_rf_2 <- confusionMatrix(pred_rf_class_2 , test_data_2_selected$Group)
print(conf_matrix_rf_2)






# Predict rf_3 
pred_rf_3 <- predict(models_featsel_3$rf, newdata = test_data_3_selected[, -ncol(test_data_3_selected)], type = "prob")
pred_rf_class_3 <- predict(models_featsel_3$rf, newdata = test_data_3_selected[, -ncol(test_data_3_selected)])

# Calculate the ROC curve
roc_curve_rf_3 <- roc(test_data_3_selected$Group, pred_rf_3[, 2], levels = levels(test_data_3_selected$Group), ci = TRUE)
auc_value_rf_3 <- auc(roc_curve_rf_3)
ci_value_rf_3<- ci.auc(roc_curve_rf_3)

# Plot the ROC curve
plot(roc_curve_rf_3, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_3, 3), ")"))
ci(roc_curve_rf_3)


print(paste("AUC:", round(auc_value_rf_3, 3)))
print(paste("910% CI:", round(ci_value_rf_3[1], 3), "-", round(ci_value_rf_3[3], 3)))

# confusion matrix
conf_matrix_rf_3 <- confusionMatrix(pred_rf_class_3 , test_data_3_selected$Group)
print(conf_matrix_rf_3)




# Predict rf_4 
pred_rf_4 <- predict(models_featsel_4$rf, newdata = test_data_4_selected[, -ncol(test_data_4_selected)], type = "prob")
pred_rf_class_4 <- predict(models_featsel_4$rf, newdata = test_data_4_selected[, -ncol(test_data_4_selected)])

# Calculate the ROC curve
roc_curve_rf_4 <- roc(test_data_4_selected$Group, pred_rf_4[, 2], levels = levels(test_data_4_selected$Group), ci = TRUE)
auc_value_rf_4 <- auc(roc_curve_rf_4)
ci_value_rf_4<- ci.auc(roc_curve_rf_4)

# Plot the ROC curve
plot(roc_curve_rf_4, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_4, 3), ")"))
ci(roc_curve_rf_4)


print(paste("AUC:", round(auc_value_rf_4, 3)))
print(paste("910% CI:", round(ci_value_rf_4[1], 3), "-", round(ci_value_rf_4[3], 3)))

# confusion matrix
conf_matrix_rf_4 <- confusionMatrix(pred_rf_class_4 , test_data_4_selected$Group)
print(conf_matrix_rf_4)

# Predict on the test set using the best model (rf_5 in this case)
pred_rf_5 <- predict(models_featsel_5$rf, newdata = test_data_5_selected[, -ncol(test_data_5_selected)], type = "prob")
pred_rf_class_5 <- predict(models_featsel_5$rf, newdata = test_data_5_selected[, -ncol(test_data_5_selected)])

# Calculate the ROC curve
roc_curve_rf_5 <- roc(test_data_5_selected$Group, pred_rf_5[, 2], levels = levels(test_data_5_selected$Group), ci = TRUE)
auc_value_rf_5 <- auc(roc_curve_rf_5)
ci_value_rf_5<- ci.auc(roc_curve_rf_5)

# Plot the ROC curve
plot(roc_curve_rf_5, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_5, 3), ")"))
ci(roc_curve_rf_5)


print(paste("AUC:", round(auc_value_rf_5, 3)))
print(paste("910% CI:", round(ci_value_rf_5[1], 3), "-", round(ci_value_rf_5[3], 3)))

# confusion matrix
conf_matrix_rf_5 <- confusionMatrix(pred_rf_class_5 , test_data_5_selected$Group)
print(conf_matrix_rf_5)



# Predict rf_6 
pred_rf_6 <- predict(models_featsel_6$rf, newdata = test_data_6_selected[, -ncol(test_data_6_selected)], type = "prob")
pred_rf_class_6 <- predict(models_featsel_6$rf, newdata = test_data_6_selected[, -ncol(test_data_6_selected)])

# Calculate the ROC curve
roc_curve_rf_6 <- roc(test_data_6_selected$Group, pred_rf_6[, 2], levels = levels(test_data_6_selected$Group), ci = TRUE)
auc_value_rf_6 <- auc(roc_curve_rf_6)
ci_value_rf_6<- ci.auc(roc_curve_rf_6)

# Plot the ROC curve
plot(roc_curve_rf_6, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_6, 3), ")"))
ci(roc_curve_rf_6)


print(paste("AUC:", round(auc_value_rf_6, 3)))
print(paste("910% CI:", round(ci_value_rf_6[1], 3), "-", round(ci_value_rf_6[3], 3)))

# confusion matrix
conf_matrix_rf_6 <- confusionMatrix(pred_rf_class_6 , test_data_6_selected$Group)
print(conf_matrix_rf_6)


# Predict rf_7 
pred_rf_7 <- predict(models_featsel_7$rf, newdata = test_data_7_selected[, -ncol(test_data_7_selected)], type = "prob")
pred_rf_class_7 <- predict(models_featsel_7$rf, newdata = test_data_7_selected[, -ncol(test_data_7_selected)])

# Calculate the ROC curve
roc_curve_rf_7 <- roc(test_data_7_selected$Group, pred_rf_7[, 2], levels = levels(test_data_7_selected$Group), ci = TRUE)
auc_value_rf_7 <- auc(roc_curve_rf_7)
ci_value_rf_7<- ci.auc(roc_curve_rf_7)

# Plot the ROC curve
plot(roc_curve_rf_7, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_7, 3), ")"))
ci(roc_curve_rf_7)


print(paste("AUC:", round(auc_value_rf_7, 3)))
print(paste("910% CI:", round(ci_value_rf_7[1], 3), "-", round(ci_value_rf_7[3], 3)))

# confusion matrix
conf_matrix_rf_7 <- confusionMatrix(pred_rf_class_7 , test_data_7_selected$Group)
print(conf_matrix_rf_7)



# Predict rf_8
pred_rf_8 <- predict(models_featsel_8$rf, newdata = test_data_8_selected[, -ncol(test_data_8_selected)], type = "prob")
pred_rf_class_8 <- predict(models_featsel_8$rf, newdata = test_data_8_selected[, -ncol(test_data_8_selected)])

# Calculate the ROC curve
roc_curve_rf_8 <- roc(test_data_8_selected$Group, pred_rf_8[, 2], levels = levels(test_data_8_selected$Group), ci = TRUE)
auc_value_rf_8 <- auc(roc_curve_rf_8)
ci_value_rf_8<- ci.auc(roc_curve_rf_8)

# Plot the ROC curve
plot(roc_curve_rf_8, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_8, 3), ")"))
ci(roc_curve_rf_8)


print(paste("AUC:", round(auc_value_rf_8, 3)))
print(paste("910% CI:", round(ci_value_rf_8[1], 3), "-", round(ci_value_rf_8[3], 3)))

# confusion matrix
conf_matrix_rf_8 <- confusionMatrix(pred_rf_class_8 , test_data_8_selected$Group)
print(conf_matrix_rf_8)


# Predict rf_9 
pred_rf_9 <- predict(models_featsel_9$rf, newdata = test_data_9_selected[, -ncol(test_data_9_selected)], type = "prob")
pred_rf_class_9 <- predict(models_featsel_9$rf, newdata = test_data_9_selected[, -ncol(test_data_9_selected)])

# Calculate the ROC curve
roc_curve_rf_9 <- roc(test_data_9_selected$Group, pred_rf_9[, 2], levels = levels(test_data_9_selected$Group), ci = TRUE)
auc_value_rf_9 <- auc(roc_curve_rf_9)
ci_value_rf_9<- ci.auc(roc_curve_rf_9)

# Plot the ROC curve
plot(roc_curve_rf_9, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_9, 3), ")"))
ci(roc_curve_rf_9)


print(paste("AUC:", round(auc_value_rf_9, 3)))
print(paste("910% CI:", round(ci_value_rf_9[1], 3), "-", round(ci_value_rf_9[3], 3)))

# confusion matrix
conf_matrix_rf_9 <- confusionMatrix(pred_rf_class_9 , test_data_9_selected$Group)
print(conf_matrix_rf_9)






# Predict rf_10 
pred_rf_10 <- predict(models_featsel_10$rf, newdata = test_data_10_selected[, -ncol(test_data_10_selected)], type = "prob")
pred_rf_class_10 <- predict(models_featsel_10$rf, newdata = test_data_10_selected[, -ncol(test_data_10_selected)])

# Calculate the ROC curve
roc_curve_rf_10 <- roc(test_data_10_selected$Group, pred_rf_10[, 2], levels = levels(test_data_10_selected$Group), ci = TRUE)
auc_value_rf_10 <- auc(roc_curve_rf_10)
ci_value_rf_10<- ci.auc(roc_curve_rf_10)

# Plot the ROC curve
plot(roc_curve_rf_10, col = "blue", main = paste0("ROC Curve (AUC = ", round(auc_value_rf_10, 3), ")"))
ci(roc_curve_rf_10)
print(paste("AUC:", round(auc_value_rf_10, 3)))
print(paste("910% CI:", round(ci_value_rf_10[1], 3), "-", round(ci_value_rf_10[3], 3)))
# confusion matrix
conf_matrix_rf_10 <- confusionMatrix(pred_rf_class_10 , test_data_10_selected$Group)
print(conf_matrix_rf_10)



#-------------------- Aggregate results--------------------------------------------------------------------------------------


# Initialize a data frame 
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

# Define a function to calculate performance metrics
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

#Calculate metrics for each model
metrics_rf_1 <- calculate_metrics(pred_rf_class_1, test_data_1_selected$Group, roc_curve_rf_1)
metrics_rf_2 <- calculate_metrics(pred_rf_class_2, test_data_2_selected$Group, roc_curve_rf_2)
metrics_rf_3 <- calculate_metrics(pred_rf_class_3, test_data_3_selected$Group, roc_curve_rf_3)
metrics_rf_4 <- calculate_metrics(pred_rf_class_4, test_data_4_selected$Group, roc_curve_rf_4)
metrics_rf_5 <- calculate_metrics(pred_rf_class_5, test_data_5_selected$Group, roc_curve_rf_5)
metrics_rf_6 <- calculate_metrics(pred_rf_class_6, test_data_6_selected$Group, roc_curve_rf_6)
metrics_rf_7 <- calculate_metrics(pred_rf_class_7, test_data_7_selected$Group, roc_curve_rf_7)
metrics_rf_8 <- calculate_metrics(pred_rf_class_8, test_data_8_selected$Group, roc_curve_rf_8)
metrics_rf_9 <- calculate_metrics(pred_rf_class_9, test_data_9_selected$Group, roc_curve_rf_9)
metrics_rf_10 <- calculate_metrics(pred_rf_class_10, test_data_10_selected$Group, roc_curve_rf_10)


# combine metrics for the 10 splits

results_pred_rf <- rbind(data.frame(Model = "rf_1", t(metrics_rf_1)),
                         data.frame(Model = "rf_2", t(metrics_rf_2)),
                         data.frame(Model = "rf_3", t(metrics_rf_3)),
                         data.frame(Model = "rf_4", t(metrics_rf_4)),
                         data.frame(Model = "rf_5", t(metrics_rf_5)),
                         data.frame(Model = "rf_6", t(metrics_rf_6)),
                         data.frame(Model = "rf_7", t(metrics_rf_7)),
                         data.frame(Model = "rf_8", t(metrics_rf_8)),
                         data.frame(Model = "rf_9", t(metrics_rf_9)),
                         data.frame(Model = "rf_10", t(metrics_rf_10))
)

summary(as.numeric(results_pred_rf$Accuracy))
print(results_pred_rf)
write.csv(results_pred_rf, "results_pred_rf_10x.csv")

results_pred_rf <- data.frame(
  Model = as.character(results_pred_rf$Model),
  Accuracy = as.numeric(results_pred_rf$Accuracy),
  P_Value_Acc_NIR = as.numeric(results_pred_rf$P_Value_Acc_NIR),
  AUC = as.numeric(results_pred_rf$AUC),
  Sensitivity = as.numeric(results_pred_rf$Sensitivity),
  Specificity = as.numeric(results_pred_rf$Specificity),
  F1_Score = as.numeric(results_pred_rf$F1_Score)
)

results_pred_rf_table <- data.frame(
  Model = as.character(results_pred_rf$Model),
  Accuracy = round(as.numeric(results_pred_rf$Accuracy),3),
  P_Value_Acc_NIR = formatC(as.numeric(results_pred_rf$P_Value_Acc_NIR), format = "e", digits = 2),
  AUC = round(as.numeric(results_pred_rf$AUC),3),
  Sensitivity = round(as.numeric(results_pred_rf$Sensitivity),3),
  Specificity = round(as.numeric(results_pred_rf$Specificity),3),
  F1_Score = round(as.numeric(results_pred_rf$F1_Score),3)
)


# Combining all roc curves into one list
roc_list_rf <- list(
  rf_1 = roc_curve_rf_1,
  rf_2 = roc_curve_rf_2,
  rf_3 = roc_curve_rf_3,
  rf_4 = roc_curve_rf_4,
  rf_5 = roc_curve_rf_5,
  rf_6 = roc_curve_rf_6,
  rf_7 = roc_curve_rf_7,
  rf_8 = roc_curve_rf_8,
  rf_9 = roc_curve_rf_9,
  rf_10 = roc_curve_rf_10
)



roc_df_rf <- data.frame()
#loop through list to obtain FPR and TPR
for (name in names(roc_list_rf)) {
  roc_data <- data.frame(
    fpr = 1 - roc_list_rf[[name]]$specificities,  # Calculate FPR correctly
    tpr = roc_list_rf[[name]]$sensitivities,
    model = name
  )
  roc_df_rf <-  rbind(roc_df_rf, roc_data)
}


#Calculate mean and confidence intervals
mean_roc_rf <- roc_df_rf %>%
  group_by(fpr) %>%
  summarise(
    mean_tpr = mean(tpr),
    ci_lower = quantile(tpr, probs = 0.025),
    ci_upper = quantile(tpr, probs = 0.975),
  )


mean_AUC_rf = mean(results_pred_rf_table$AUC)
ci_lower_AUC_rf = quantile(results_pred_rf_table$AUC, probs = 0.025)
ci_upper_AUC_rf = quantile(results_pred_rf_table$AUC, probs = 0.975)


# Combining all selected features into one list
all_selected_features <- c(
  selected_features_rf_1,
  selected_features_rf_2,
  selected_features_rf_3,
  selected_features_rf_4,
  selected_features_rf_5,
  selected_features_rf_6,
  selected_features_rf_7,
  selected_features_rf_8,
  selected_features_rf_9,
  selected_features_rf_10
)


selected_features_count <- table(all_selected_features)

#order the data frame by the Count column in descending order
selected_features_df <- selected_features_df[order(-selected_features_df$Count), ]


print(selected_features_df)


