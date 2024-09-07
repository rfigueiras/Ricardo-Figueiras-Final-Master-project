
#BiocManager::install('mixOmics')
library(mixOmics)
library(ggplot2)
library(dplyr)
library(tidyr)

# Code used for MIXOMICs integration and plots using features previously selected by
# recursive feature elimination.

# Define variables categories and perform N-integration using Partial Least Squares 
#Discriminant Analysis (PLS-DA)

data_proteinsAFP_matched_sel<-rbind(data_proteins_matched_sel_anonimised, data_AFP_matched_sel_anonimised)

data_list_sel_AFP <- list(protein_AFP= t(data_proteinsAFP_matched_sel), miRNA = t(data_miRNA_matched_sel_anonimised))
result.diablo_sel_AFP <- block.plsda(data_list_sel_AFP, group, scale = TRUE)

data_list_sel <- list(protein= t(data_proteins_matched_sel_anonimised), miRNA = t(data_miRNA_matched_sel_anonimised))
result.diablo_sel <- block.plsda(data_list_sel, group, scale = TRUE)


# Plot variable contributions for protein_AFP block
plotVar(result.diablo_sel_AFP)


png("Mixomics correlation circle plot all_labels.png", width = 3000, height = 2500, res=400)

plotVar(result.diablo_sel_AFP)

dev.off()


# Plot variable representation plots
plotVar(result.diablo_sel), block = 'miRNA')
plotVar(result.diablo_sel, block = 'protein')

png("Mixomics plotVar all_labels.png", width = 1200, height = 1000, res=200)

# Plot patients scatter plot
plotIndiv(result.diablo_sel, group = group, ind.names = FALSE, legend = TRUE)

dev.off()




png("Circos plot selected variables_labels.png", width = 4000, height = 3500, res=600)

# Circular plot
circosPlot(result.diablo_sel_AFP, cutoff = 0.70, line = TRUE, 
           color.blocks = c('salmon3', 'dodgerblue2'),
           color.cor = c("firebrick4", "bisque3"), size.labels = 0.9, size.variables = 0.,
           showIntraLinks = FALSE, ncol.legend = 2)

dev.off()



# Get DIABLO correlations
corMat <- circosPlot(result.diablo_sel_AFP, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1)
corMat



# Extract data for miRNA and protein blocks to perform Pearson correlation
miRNA_data <- as.data.frame(t(data_miRNA_matched_sel_anonimised))  # miRNA block data
protein_data <- as.data.frame(t(data_proteins_matched_sel_anonimised))  # protein block data


# Calculate correlation matrix
cor_matrix <- cor(result.diablo_sel_AFP$X$protein_AFP, result.diablo_sel_AFP$X$miRNA, use = "complete.ob", method = "pearson")
cor_matrix <- cor(t(data_proteinsAFP_matched_sel), t(data_miRNA_matched_sel_anonimised), use = "complete.ob", method = "pearson")

# View the correlation matrix
print(cor_matrix)


write.csv(cor_matrix,"correlation_matrix.csv")

write.csv(corMat, "diablo_correlations.csv")



