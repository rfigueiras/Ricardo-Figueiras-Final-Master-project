# library

remotes::install_github("ruma1974/EVqualityMS@master")
library("ruma1974/EVqualityMS@master")
library(EVqualityMS)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOstats)
library(GO.db)
library(RColorBrewer)
library(enrichR)
library(dplyr)
library(pheatmap)
library(enrichR)

#--------------------------------EV_quality_MS_tool---------------------------------------------------------------------------

# Data input

MS_batch1<-read.csv("MS_batch1_fp.csv", check.names = FALSE)
MS_batch2<-read.csv("MS_batch2_fp.csv", check.names = FALSE)
MS_batch3<-read.csv("MS_batch3_fp.csv", check.names = FALSE)

metadata_batch1<-read.csv("metadata_filt_batch1.csv", check.names = FALSE)
metadata_batch2<-read.csv("metadata_filt_batch2.csv", check.names = FALSE)
metadata_batch3<-read.csv("metadata_filt_batch3.csv", check.names = FALSE)


# Change protein group name to GN to be compatible with EVqualityMS tool
colnames(MS_batch1)<-c("GN", colnames(MS_batch1[,-1]))
colnames(MS_batch2)<-c("GN", colnames(MS_batch2[,-1]))
colnames(MS_batch3)<-c("GN", colnames(MS_batch3[,-1]))

merged_data <- merge(MS_batch1, MS_batch2, by = "GN", all = FALSE)
merged_data <- merge(merged_data, MS_batch3, by = "GN", all = FALSE)

# View the merged data
head(merged_data)
dim(merged_data)


#merge batches by protein group names
merged_data <- merge(MS_batch1, MS_batch2, by = "GN", all = FALSE)
merged_data <- merge(merged_data, MS_batch3, by = "GN", all = FALSE)

# View the merged data
head(merged_data)
dim(merged_data)

#heatmap with subpopulation EV markers from Kowal et al 2016 paper for each 
#batch

png("EVsub_b1.png", height=1200, width=1600, res=300)
EVqualityMS::heatmapEVsub(MS_batch1)
dev.off()

png("EVsub_b2.png", height=1200, width=1600, res=300)
EVqualityMS::heatmapEVsub(MS_batch2)
dev.off()

png("EVsub_b3.png", height=1200, width=1600, res=300)
EVqualityMS::heatmapEVsub(MS_batch3)
dev.off()

#heatmap with Protein markers co-isolated with CD63, CD81, CD9 in Kowal et al 2016 
#paper using protein groups common in 3 batches
png("EVsub_CD_all.png", height=1200, width=1600, res=300)
EVqualityMS::heatmapEVsubCD(merged_data)
dev.off()

#--------------------KEGG pathway analysis-------------------------------------------------------------------------------------------#
  
# Data frame with Protein groups that intersect in the 3 batches  
Protein_Groups_int<-intersect(intersect(MS_batch1$GN, MS_batch2$GN), MS_batch3$GN)


Protein_Groups_int<-gsub("\\(.*?\\)","", Protein_Groups_int)


# Convert for EntrezID
eg_protein_group<-bitr(Protein_Groups_int, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
OrgDb <- org.Hs.eg.db
universe <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# enrichKEGG adjusted p-val 0.05
ego_MS_KEGG <-enrichKEGG(
  eg_protein_group$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  #universe,
  minGSSize = 10,
  maxGSSize = 200,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)


# Dotplot with results
p_KEGG <- dotplot(ego_MS_KEGG, title = "KEGG pathway analysis", showCategory = 25, font.size = 15)

p_KEGG <-p_KEGG + theme(
  axis.text.y = element_text(size = 24),        
  axis.text.x = element_text(size = 28),        
  plot.title = element_text(size = 30),  
  legend.text = element_text(size = 25),       
  legend.title = element_text(size = 24),       
  axis.title.x = element_text(size = 28),       
  axis.title.y = element_text(size = 28),
  plot.margin = margin(t = 15, r = 15, b = 15, l = 140)
)

ggsave("KEGG_pathway.png", plot = p_KEGG, width = 12, height = 20, dpi=600)

KEGG_MS<-ego_MS_KEGG@result
write.csv(KEGG_MS, "KEGG_MS.csv")

write.csv(Protein_Groups_int, "Proteins_present_in_all_batches.csv")

#----------------ENRICHR GTEX tissue upregulated markers-----------------------------------------------------------------------------------------

# avaialable databases 
dbs <- listEnrichrDbs()
head(dbs)


# Apply EnrichR using GTEx_Tissue_Expression_Up database
results_NEXPLOR_GTEx_Tissue_Expression_Up <- enrichr(eg_protein_group$SYMBOL, c("GTEx_Tissue_Expression_Up"))



p_gtex2<-plotEnrich(results_NEXPLOR_GTEx_Tissue_Expression_Up$GTEx_Tissue_Expression_Up, showTerms = 12, numChar = 45, y = "Count", orderBy = "P.value",
                   title = "GTEx_Tissue_Expression_Up")

p_gtex2 <- p_gtex2 + theme(
  axis.text.y = element_text(size = 30),      
  axis.text.x = element_text(size = 30),      
  plot.title = element_text(size = 30),  
  legend.text = element_text(size = 26),      
  legend.title = element_text(size = 26),       
  axis.title.x = element_text(size = 28),       
  axis.title.y = element_text(size = 28)      
)


ggsave("Enrichr GTEx_Tissue_exp_UP.png", plot = p_gtex2, width = 22, height = 12, dpi=400)

#--------------------Heatmaps of MISEV2023 markers-------------------------------#
MS_data_batch1<- as.data.frame(MS_batch1)
MS_data_batch2<- as.data.frame(MS_batch2)
MS_data_batch3<- as.data.frame(MS_batch3)

rownames(MS_data_batch1)<-MS_data_batch1$GN
rownames(MS_data_batch2)<-MS_data_batch2$GN
rownames(MS_data_batch3)<-MS_data_batch3$GN

MS_data_batch1<-MS_data_batch1[,-1]
MS_data_batch2<-MS_data_batch2[,-1]
MS_data_batch3<-MS_data_batch3[,-1]

# Perform outer join on the row names
merged_data <- merge(merge(MS_data_batch1, MS_data_batch2, by = "row.names", all = TRUE), MS_data_batch3, by.x = "Row.names", by.y = "row.names", all = TRUE)
rownames(merged_data) <- merged_data$Row.names
merged_data <- merged_data[, -which(names(merged_data) == "Row.names")]


metadata_merged<- rbind(metadata_batch1, metadata_batch2, metadata_batch3)

# Protein names for each category of MISEV 2023 markers
proteins_1a <- c("\\bCD9\\b", "CD63", "CD81", "CD82", "CD47", "GNA.*", "TSAP6")
proteins_1b <- c("ITGA.*", "ITGB.*", "TFR2", "LAMP1", "LAMP2", "\\bSDC.*\\b", "BSG", "ADAM10")
proteins_1c <- c("GPC1", "NT5E", "CD59")
proteins_2a <- c("\\bTSG101\\b", "\\bCHMP.*\\b", "\\bPDCD6IP\\b", "\\bVPS4A\\b", "\\bVPS4B\\b", "\\bFLOT1\\b", "\\bFLOT2\\b", "\\bCAV.*\\b", "\\bSDCBP\\b")
proteins_2b <- c("\\bHSPA8\\b", "\\bHSP90AB1\\b", "\\bACT(?!N)\\w*\\b", "\\bTUB.*\\b", "\\bGAPDH\\b")
proteins_3a <- c("\\bAPO.*\\b")
proteins_3b <- c("\\bIG.*\\b", "\\bUMOD\\b", "\\bALB\\b", "\\b14-3-3.*\\b", "\\bAGO.*\\b")
proteins_3c <- c("\\bHSP90AA\\b", "\\bHSP90AB\\b", "\\bTGFBI\\b", "\\bHSPA13\\b", "\\bLDHA\\b", "\\bLDHB\\b")
proteins_4a <- c("\\bHIST1H.*\\b", "\\bLMNA\\b", "\\bLMNC\\b")
proteins_4b <- c("\\bVDAC.*\\b", "\\bCYC1\\b", "\\bTOMM20\\b")
proteins_4c <- c("\\bCANX\\b", "\\bHSP90B1\\b", "\\bHSPA5\\b", "\\bGOLGA2\\b")
proteins_4d <- c("\\bMAP1LC3A\\b", "\\bACTN1\\b", "\\bACTN4\\b")
proteins_5a <- c("\\bAPO\\b.*", "\\bC[1-9].\\b", "\\bFGB\\b", "\\bFGA\\b", "\\bFGG\\b")
proteins_5b <- c("\\bTGFB1\\b", "\\bTGFB2\\b", "\\bIFNG\\b", "\\bVEGFA\\b", "\\bFGF1\\b", "\\bFGF2\\b", "\\bPDGF.*\\b", "\\bEGF\\b", "\\bIL.*\\b")
proteins_5c <- c("\\bFN1\\b", "\\bCOL.*\\b", "\\bMFGE8\\b", "\\bLGALS3BP\\b", "\\bCD5L\\b", "\\bAHSG\\b")


# Function using regex to filter data based on the lists of protein groups names
filter_data_regex <- function(data, proteins) {
  pattern <- paste(proteins, collapse = "|")
  filtered_data <- data[grepl(pattern, rownames(data), ignore.case = TRUE, perl = TRUE), ]
  return(filtered_data)
}


merged_data_null_imp<-merged_data
#change missing values to 0 for vizualisation purposes
merged_data_null_imp[is.na(merged_data_null_imp)]<-0

# Create subsets for each category
data_1a <- filter_data_regex(merged_data_null_imp, proteins_1a)
data_1b <- filter_data_regex(merged_data_null_imp, proteins_1b)
data_1c <- filter_data_regex(merged_data_null_imp, proteins_1c)
data_2a <- filter_data_regex(merged_data_null_imp, proteins_2a)
data_2b <- filter_data_regex(merged_data_null_imp, proteins_2b)
data_3a <- filter_data_regex(merged_data_null_imp, proteins_3a)
data_3b <- filter_data_regex(merged_data_null_imp, proteins_3b)
data_3c <- filter_data_regex(merged_data_null_imp, proteins_3c)
data_4a <- filter_data_regex(merged_data_null_imp, proteins_4a)
data_4b <- filter_data_regex(merged_data_null_imp, proteins_4b)
data_4c <- filter_data_regex(merged_data_null_imp, proteins_4c)
data_4d <- filter_data_regex(merged_data_null_imp, proteins_4d)
data_5a <- filter_data_regex(merged_data_null_imp, proteins_5a)
data_5b <- filter_data_regex(merged_data_null_imp, proteins_5b)
data_5c <- filter_data_regex(merged_data_null_imp, proteins_5c)

length(metadata_merged$Sample_ID)
length(colnames(merged_data_null_imp))



#heatmaps with pheatmap
annotation<-data.frame(group=metadata_merged$Group,
                       batch=metadata_merged$MS_batch)

rownames(annotation)<-colnames(merged_data)

colnames(data_1a)

annotation$group <- factor(annotation$group, levels = c("HCC", "Cirrhosis"))

annotation$batch <- factor(annotation$batch, levels = c("Batch1", "Batch2", "Batch3"))

color_group<-c("orange","cyan")

color_phenotype<-c(brewer.pal(length(unique(metadata_merged$MS_batch)), "Set1"))

annotation_color<-list(group = setNames(color_group, levels(annotation$group)),batch= setNames(color_phenotype, levels(annotation$batch)))

# Transmembrane EV proteins
data_1subset<-rbind(data_1a[c("CD81", "GNAI2"),],data_1b[c("ADAM10", "LAMP1","LAMP2"),])

p1s<-pheatmap(data_1subset,
         annotation = annotation,
         annotation_color = annotation_color,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE, main = "MISEV2023 Markers 1  \n      Transmembrane proteins in EVs",
         fontsize=28, show_colnames = FALSE)

ggsave("MISEV2023 Markers 1 subset.png", plot = p1s, width = 12, height = 12, dpi=600)


# Citosolic EV proteins
rownames(data_2a)<- gsub("\\(.*?\\)", "", rownames(data_2a))
p2a<-pheatmap(data_2a,
         annotation = annotation,
         annotation_color = annotation_color,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE, main = "MISEV2023 Markers 2 \nCytosolic proteins in EVs",
         fontsize=28, show_colnames = FALSE)


ggsave("MISEV2023 Markers 2.png", plot = p2a, width = 12, height = 12, dpi=600)



