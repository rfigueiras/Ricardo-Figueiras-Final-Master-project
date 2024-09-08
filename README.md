# Ricardo-Figueiras-Final-Master-Project  
## Master in Omics Data Analysis, Universitat de Vic (UVIC)  
Some of the code used for the final master project: "Multiomic Profiling of Circulating Hepatocyte-Derived Extracellular Vesicles for Hepatocellular Carcinoma Detection in Liver Cirrhosis."

### 1- Proteomics Analysis:  
1- `Proteomics_functional_analysis.R`: Includes analysis with “EVqualityMS” (https://github.com/ruma1974/EVqualityMS), KEGG pathway enrichment analysis, liver marker enrichment according to the GTEX database, and heatmaps with MISEV2023 markers.  
2- `Proteomic_pipeline3.R`: Preprocessing with the selected pipeline 3. Similar analyses and code were applied to the other 4 pipelines tested, with some differences described here:  
<img src="https://github.com/rfigueiras/Ricardo-Figueiras-Final-Master-project/blob/main/1-Proteomics_analysis/Pipelines%20test%20for%20proteomics%20preprocessing.png?raw=true" alt="Proteomics preprocessing pipelines tested" width="50%">  

3- `pvca_subset_test.R`: Principal variance component analysis with resampling to evaluate the different pipelines.  
4- `Testing_proteomic_signature_on_batch_3_samples.R`: Prediction of batch 3 samples using the proteomic signature developed by training on batch 1 and batch 2 samples.  

### 2- exceRpt microRNA Analysis:  
1- `Sequencing_alignment_with_exceRpt_Docker.txt`: Code used to apply the exceRpt Docker to several samples in parallel.    
2- `analysis_microRNA_exceRpt.R`: Preprocessing of exceRpt microRNA data. The same script was also applied to the COMPSRA microRNA data.

### 3- Multiomics Integration and Testing:  
1- `mixOmics_integration_and_plots.R`: Code used for N-integration with mixOmics and circular plots.  
2- `Training_and_testing_with_random_data_partitioning.R`: Random data partitioning without replacement. Scaling, recursive feature elimination, and random forest model training were applied to the training subset using the optimal subset of selected features, which were then used to predict classes in the testing subset.

### 4- Cross-validation and Recursive Feature Elimination:  
1- `Cross_validation_and_recursive_feature_elimination.R`: Code used for testing with cross-validation using all features, and for feature selection with recursive feature elimination. This approach was used for the proteomics, microRNA, and multiomics data.

