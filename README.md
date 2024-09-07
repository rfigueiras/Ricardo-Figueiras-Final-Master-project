# Ricardo-Figueiras-Final-Master-project  
## Master in Omics Data Analysis, Universitat de Vic (UVIC)  
Some of the code used for the final master project: "Multi center Targeted Multiomic Profiling of Plasma Circulating Extracellular Vesicles for Hepatocellular Carcinoma Detection in Liver Cirrhosis" including:  

### 1- Proteomics analysis:  
1- Proteomics_functional_analysis>R: includes analysis with “EVqualityMS” (https://github.com/ruma1974/EVqualityMS), KEGG pathway enrichment analysis, liver markers enrichment according to GTEX database and heatmaps with MIVEV2023 markers.  
2- Proteomic_pipeline3.R: Preprocessing with the selected pipeline from the different 5 pipelines tested:  
<img src="https://github.com/rfigueiras/Ricardo-Figueiras-Final-Master-project/blob/main/1-Proteomics_analysis/Pipelines%20test%20for%20proteomics%20preprocessing.png?raw=true" alt="Proteomics preprocessing pipelines tested" width="50%">  

3- pvca_subset_test.R: Principal variance component analysis with resampling to evaluate the different pipelines.  
4- Testing_proteomic_signature_on_batch_3_samples: Prediction of batch 3 samples with the proteomic signature developed by using batch 1 and batch 2 samples as the training dataset.   
### 2- exceRpt microRNA analysis:  
1- Sequencing alignment with exceRpt Docker.txt: Code used to apply the exceRpt Docker to several samples in paralalel.    
2- analysis_microRNA_exceRpt.R: Preprocessing of microRNA data.  
### 3- Multiomics integration and testing:  
1- mixOmics integration and plots.R: Code used for N-integration with mixOmics and circular plot.  
2- Training and testing with random data partitioning.R: Random data partitioning without replacement. Scaling, RFE, and random forest model training for the respective optimal subset of selected features was applied to the training subset and then used to predict the classes in the testing subset.   
### 4- Cross-validation and Recursive Feature Elimination:    
1- Cross-validation and Recursive Feature elimination.R: Code used for testing with cross-validation with all features and for recursive feature elimination. This approach was used for the proteomics, microRNA and multiomics data.  
