 

# Example of code used for Cross-validation and Recursive feature elimination 
# here using microRNA data, similar code was used for the proteomics and multiomics data




# Data loading

Seurat_counts<- read.csv("Seurat_counts_filt_miRNA.csv", row.names = 1, check.names = FALSE)
metadata<-read.csv("metadata_filtered_RNA.csv", row.names = 1, check.names = FALSE)
metadata_filt<-metadata[metadata$Group!="QC",]
Seurat_counts_filt<-Seurat_counts[,colnames(Seurat_counts) %in% metadata_filt$Sample]


#Library  
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
library(dplyr)
library(tidyr)
library(patchwork)
library(gridExtra)
library(grid)





#-------------------------------------------- Cross-validation with all features ---------------------------------------------------------------------------

# Get labels
train_labels <- metadata_filt$Group
train_df <- as.data.frame(t(Seurat_counts_filt))
# Add labels to the training data
train_df$Group <- as.factor(train_labels)
# Use HCC as reference level
train_df$Group <- relevel(train_df$Group, ref = "HCC")

# Core parallelization 
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Train models cross-validation
set.seed(1534)
train_control <- trainControl(method = "repeatedcv", number = 5, repeats = 100, savePredictions = "final", classProbs = FALSE)

models <- caretList(
  Group ~ ., data = train_df,
  trControl = train_control,
  methodList = c("rf", "svmRadial", "nnet", "glmnet", "knn", "lda"))#algorithms to test

results <- resamples(models)
summary_results <- summary(results)


#get summary results for Accuracy and Kappa
acc_mirnas<-summary_results$statistics$Accuracy
kappa_mirnas<-summary_results$statistics$Kappa

acc_mirnas<-as.data.frame(acc_mirnas)
kappa_mirnas<-as.data.frame(kappa_mirnas)

acc_mirnas<-acc_mirnas[,-7]
kappa_mirnas<-kappa_mirnas[,-7]

write.csv(acc_mirnas, "miRNAs_cross_validation_acc.csv")

write.csv(kappa_mirnas, "miRNAs_cross_validation_kappa.csv")


acc_mirnas$models<-rownames(acc_mirnas)
kappa_mirnas$models<-rownames(kappa_mirnas)



# Plot summary results of cross-validation
mirnas_long_acc <- acc_mirnas %>%
  pivot_longer(cols = Min.:Max., names_to = "Statistic", values_to = "Accuracy")

mirnas_long_Kappa <- kappa_mirnas %>%
  pivot_longer(cols = Min.:Max., names_to = "Statistic", values_to = "Kappa")


plot_accuracy_mirnas <- ggplot(mirnas_long_acc, aes(x = Accuracy, y = models)) +
  geom_boxplot(fill = "beige") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Accuracy",
       y = "Model Type")+ theme(axis.title = element_text(size = 30),
                                axis.text.y = element_text(size = 30, face="bold"),
                                axis.text.x = element_text(size = 30))

plot_kappa_mirnas <- ggplot(mirnas_long_Kappa, aes(x = Kappa, y = models, fill = pipeline)) +
  geom_boxplot(fill = "beige") +
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  labs(x = "Kappa",
       fill = "Pipeline")+ theme(axis.title = element_text(size = 30),
                                 axis.text = element_text(size = 30)) 

# Combine the plots
combined_plot_mirna <- plot_accuracy_mirnas + plot_kappa_mirnas + plot_layout(nrow = 1) +
  plot_annotation(title = "MicroRNA cross-validation training dataset \nafter RFE (k= 5 folds, 100 repeats)",
                  theme = theme(plot.title = element_text(size = 32, hjust = 0.5), panel.grid.major = element_line(size = 1),  
                                panel.grid.minor = element_line(size = 1), 
                                panel.border = element_rect(size = 1)       
                  ))

# Display the combined plot
print(combined_plot_mirna)

ggsave("Accuracies and kappa cross-validation training miRNAs.png", plot = combined_plot_mirna, width = 12, height = 12, dpi=600)



#-----------------------------------------------------Feature selection---------------------------------------------------------------------#

# Core parallelization 
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
# RFE feature selection
rfe_control_rf <- caret::rfeControl(functions = rfFuncs, method = "repeatedcv", number = 5, repeats = 100, returnResamp = "final", allowParallel = TRUE, verbose = TRUE)
#option "sizes = c(1:20, 25, 35:40)" was used for proteomics and multiomics datasets

set.seed(1534)
rfe_results <- caret::rfe(train_df[, -ncol(train_df)], train_df$Group, rfeControl = rfe_control_rf)

png("Micro RNA Recursive Feature Elimination.png", width = 1000, height = 1200, res=200)
plot(rfe_results)
dev.off()



# Extract the RFE results and create a data frame for ggplot
rfe_df <- data.frame(Variables = rfe_results$results$Variables, 
                     Accuracy = rfe_results$results$Accuracy)

# Plotting results
plot1 <- ggplot(rfe_df, aes(x = Variables, y = Accuracy)) +
  geom_point(size = 5, color = "dodgerblue4") +
  geom_vline(xintercept = 16, linetype = "dashed", color = "indianred", size = 1) +
  annotate("text", x = 16+2, y = 0.625, label = "Selected", color = "indianred", angle = 90, vjust = 0.5, hjust = 1, size=9) +
  scale_x_continuous(breaks = seq(0, 58, by = 5), limits = c(0, 60)) +
  labs(x = "Number of variables", y = "Accuracy", title = NULL) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 30),
        axis.text.y = element_text(size = 30),
        panel.grid.major.y = element_blank(), 
        panel.grid.major = element_line(size = 1),  
        panel.grid.minor = element_line(size = 1), 
        panel.border = element_rect(size = 1)       
  )



plot1 <- plot1 + ggtitle("MicroRNA \nRecursive Feature Elimination (RFE) Result") + theme(plot.title = element_text(size = 32, hjust = 0.5))

ggsave("MicroRNA Recursive Feature Elimination.png", plot = plot1, width = 12, height = 12, dpi=400) 





write.csv(rfe_results$results, "RF_variable_acc_microRNAs.csv")



selected_features <- predictors(rfe_results)
write.csv(selected_features, "rfe_selected_features_microRNAs.csv")


#---------------------------- Plot rankings of variable importance during RFE------------------------------------------------------#

microrna_name_map<-read.csv("microRNA_counts_anonimized.csv")

#get the importance scores for the selected variables
importance_rfe <- varImp(rfe_results)
importance_rfe$microrna<-rownames(importance_rfe)
importance_rfe <- importance_rfe[rownames(importance_rfe) %in% selected_features,]

#Replace the row names with their corresponding anonymous name
importance_rfe_anonimised <- importance_rfe
name_mapping <- setNames(microrna_name_map$Name, microrna_name_map$Real_name)
for (real_name in rownames(importance_rfe_anonimised)) {
  if (real_name %in% names(name_mapping)) {
    rownames(importance_rfe_anonimised)[rownames(importance_rfe_anonimised) == real_name] <- name_mapping[real_name]
  }
}


# Data frame for plotting
importance_rfe_anonimised_df <- data.frame(microrna = rownames(importance_rfe_anonimised), 
                                          Importance = importance_rfe_anonimised[1])  
  
# Plot ranking  
feat_imp<-ggplot(importance_rfe_anonimised_df, aes(x = reorder(microrna, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "dodgerblue4") +
  coord_flip() +  # Flip the coordinates for better readability
  labs(title = "Feature Importance RFE",
       x = "miRNAs",
       y = "Importance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 32),
    axis.text.y = element_text(size = 32, face="bold"),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    plot.title = element_text(size = 30, hjust = 0.5)
    , panel.grid.major = element_line(size = 1),  
    panel.grid.minor = element_line(size = 1), 
    panel.border = element_rect(size = 1)       
  )

  
ggsave("Variable importance RFE microRNA.png", plot = feat_imp, width = 12, height = 18, dpi=600)





