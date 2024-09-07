
#library(remotes)
#remotes::install_github("dleelab/pvca")
library(pvca)
library(Matrix)
library(lme4)
library(ggplot2)

#Data loading

merged_combat_knn_vsn<- read.csv("batch_merged_combat_knn_vsn.csv", row.names = 1, check.names = FALSE)#pipeline 1
knn_vsn_combat<-read.csv("knn_vsn_combat.csv", row.names = 1, check.names = FALSE)#pipeline 2      
vsn_knn_combat<-read.csv("vsn_knn_combat.csv", row.names = 1, check.names = FALSE)#pipeline 3   
minprob_qsmooth_combat<- read.csv("minprob_qsmooth_combat.csv", row.names = 1, check.names = FALSE)#pipeline 4
minprob_vsn_combat<- read.csv("minprob_vsn_combat.csv", row.names = 1, check.names = FALSE)#pipeline 5
metadata_merged<-read.csv("metadata_merged.csv", row.names = 1)#metadata

# selecting variable to be assessed 
valid_technical_factors <- c("MS_batch", "Blood.fractionation.protocol", "Prep_day", "Rounded.volume")
valid_biological_factors <- c("Group", "Etiology2", "Cirrhosis..compensation.status", "X.Portal.hypertension" )


all_factors <- c(valid_technical_factors, valid_biological_factors)

#subset the metadata for the variables of interest
PVCA_metadata <- metadata_merged[,all_factors]
PVCA_metadata$Prep_day<-as.factor(PVCA_metadata$Prep_day)
colnames(PVCA_metadata) <- make.names(colnames(PVCA_metadata))




#custom_PVCA function adapted from https://github.com/dleelab/pvca
custom_PVCA <- function (counts, meta, threshold, inter, already_normalized = TRUE) 
{  # if the counts matrix is not normalized, center the data
  if (!already_normalized) {
    counts.center <- t(apply(counts, 1, scale, center = TRUE, scale = FALSE))
  } else {
    counts.center <- counts # else proceed
  }
  #obtain correlation matrix
  cor.counts <- cor(counts.center)
  print("Correlation matrix calculated.")
  #eigen decomposition on the correlation matrix
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  
  print("Eigen decomposition done.")
  #proportion of variance explained by each principal component.
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val / eigen.val.sum
  
  meta <- as.data.frame(meta)#convert metadata into data.frame
  all <- 0
  npc.in <- 0
  #number of principal components (PCs) needed to explain at least the threshold of variance
  for (i in 1:n.eigen) {
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if (all > threshold) {
      break
    }
  }
  if (npc.in < 3) {
    npc.in <- 3# minumum number of PCs
  }
  pred.list <- colnames(meta)
  meta <- droplevels(meta)
  n.preds <- ncol(meta) + 1
  if (inter) {
    n.preds <- n.preds + choose(ncol(meta), 2)
  }
  ran.pred.list <- c()
  for (i in 1:ncol(meta)) {
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ")"))
  }
  if (inter) {
    for (i in 1:(ncol(meta) - 1)) {
      for (j in (i + 1):ncol(meta)) {
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep = " ~ ")
  ran.var.mat <- NULL
  
  print("Starting mixed model fitting.")#  fitting mixed-effects model for each selected PC.
  
  for (i in 1:npc.in) {
    dat <- data.frame(pc = eigen.mat[, i], meta)
    
    print(paste("Processing PC", i, "with data dimension:", dim(dat)))
    
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit)
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  
  print("Mixed model fitting completed.")
  #calculate the proportion of variance explained by each component
  ran.var.mat.std <- ran.var.mat / rowSums(ran.var.mat)
  wgt.vec <- eigen.val / eigen.val.sum
  prop.var <- colSums(ran.var.mat.std * wgt.vec[1:npc.in])
  std.prop.var <- prop.var / sum(prop.var)
  
  print("Variance proportion calculation completed.")
  
  return(std.prop.var)
}


#function to perform PVCA (using the previous function) on random subsets of the 
#data with resanmpling done with replacement:
random_subset_PVCA <- function(counts, meta, threshold, inter, already_normalized = TRUE, n_subsets = 100, subset_size = 0.8) {
  results <- list()
  
  for (i in 1:n_subsets) {
    #generate a random subset of indices
    subset_indices <- sample(1:ncol(counts), size = floor(subset_size * ncol(counts)), replace = TRUE)#rounding down number of samples, resampling with replacement
    
    #subset the counts and metadata with the selected indices
    counts_subset <- counts[,subset_indices ]
    meta_subset <- meta[subset_indices, ]
    
    #apply custom_PVCA function to the subset
    pvca_result <- custom_PVCA(counts_subset, meta_subset, threshold, inter, already_normalized)
    
    # Store the result
    results[[i]] <- pvca_result
  }
  
  #aggregate results
  aggregate_results <- do.call(rbind, results)
  mean_results <- colMeans(aggregate_results)
  sd_results <- apply(aggregate_results, 2, sd)
  
  return(list(mean = mean_results, sd = sd_results, all_results = results))
}




#perform PVCA on 100 random subsets of 90% of the samples using principal components that explain of 75% of variance 


#pipeline 1
PVCA_merged_combat_knn_vsn<-random_subset_PVCA(as.data.frame(merged_combat_knn_vsn), PVCA_metadata, threshold = 0.75, inter = FALSE, already_normalized = FALSE, n_subsets = 100, subset_size = 0.9)
#pipeline 2
PVCA_knn_vsn_combat<-random_subset_PVCA(as.data.frame(knn_vsn_combat), PVCA_metadata, threshold = 0.75, inter = FALSE, already_normalized = FALSE, n_subsets =100, subset_size = 0.9)
#pipeline 3
PVCA_vsn_knn_combat<-random_subset_PVCA(as.data.frame(vsn_knn_combat), PVCA_metadata, threshold = 0.75, inter = FALSE, already_normalized = FALSE, n_subsets = 100, subset_size = 0.9)
#pipeline 4
PVCA_minprob_vsn_combat<-random_subset_PVCA(as.data.frame(minprob_vsn_combat), PVCA_metadata, threshold = 0.75, inter = FALSE, already_normalized = FALSE, n_subsets = 100, subset_size = 0.9)
#pipeline 5
PVCA_minprob_qsmooth_combat<-random_subset_PVCA(as.data.frame(minprob_qsmooth_combat), PVCA_metadata, threshold = 0.75, inter = FALSE, already_normalized = FALSE, n_subsets = 100, subset_size = 0.9)



#get results for biological variable group
group_values_merged_combat_knn_vsn <- sapply(PVCA_merged_combat_knn_vsn$all_results, function(result) result["Group"])
group_values_knn_vsn_combat <- sapply(PVCA_knn_vsn_combat$all_results, function(result) result["Group"])
group_values_vsn_knn_combat<- sapply(PVCA_vsn_knn_combat$all_results, function(result) result["Group"])
group_values_minprob_vsn_combat <- sapply(PVCA_minprob_vsn_combat$all_results, function(result) result["Group"])
group_values_minprob_qsmooth_combat <- sapply(PVCA_minprob_qsmooth_combat$all_results, function(result) result["Group"])

group_values_list <- list(
  pipeline1 = group_values_merged_combat_knn_vsn,
  pipeline2 = group_values_knn_vsn_combat,
  pipeline3 = group_values_vsn_knn_combat,
  pipeline4 = group_values_minprob_vsn_combat,
  pipeline5 = group_values_minprob_qsmooth_combat
)

group_values_df <- do.call(rbind, lapply(names(group_values_list), function(name) {
  data.frame(VarianceExplained_group = group_values_list[[name]], Dataset = name)
}))



# Create boxplots using ggplot2
library(ggplot2)

pvca_comp<-ggplot(group_values_df, aes(x = Dataset, y = VarianceExplained_group, fill= Dataset)) +
  geom_boxplot() +
  labs(title = "PVCA subsets test \nvariance Explained by 'Group' Variable",
       x = "Dataset",
       y = "Variance Explained")+theme_bw() +
  theme(
    axis.title.y = element_text(size = 26, hjust = 0.5), 
    axis.title.x = element_text(size = 26, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 26),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    plot.title = element_text(size = 28, hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    plot.caption = element_text(size = 26, hjust = 0.5)
  )

#perform Kruskal-Wallis test
kruskal_test <- kruskal.test(VarianceExplained_group ~ Dataset, data = group_values_df)
print(kruskal_test)

#perform pairwise Wilcoxon test (post-hoc analysis)
pairwise_wilcox <- pairwise.wilcox.test(group_values_df$VarianceExplained_group, group_values_df$Dataset, p.adjust.method = "BH")
print(pairwise_wilcox)


as.data.frame(pairwise_wilcox$p.value)

ggsave("PVCA subsets Group explained.png", plot = pvca_comp, width = 8, height = 8, dpi=300)