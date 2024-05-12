## ecDNA gene expression

Lin, M.S.; Jo, S., Luebeck, J.; Chang, H.Y.; Wu, S.; Mischel, P.S.; Bafna, V. “Transcriptional immune suppression and upregulation of double stranded DNA damage and repair repertoires in ecDNA-containing tumors.” eLife. 12:RP88895 (Reviewed preprint). February 27, 2024. https://doi.org/10.7554/eLife.88895.2.

**Impute missing values using K-nearest neighbors (KNN)**
```
python -u ./src/cBioPortal_data_preprocessing.py --cBioPortal_data ./data/reference/TCGA/cBioportal --output ./data/KNN_imputed_data
```

**Generate gene expression matrices**
```
python -u ./src/Generate_TCGA_KNN_matrix.py --minimum_ecdna 3 --protein_coding 1 --KNN_imputed_data_dir ./data/KNN_imputed_data --output ./data/RNAseq_expression_matrix
```
```
python -u ./src/Generate_TCGA_KNN_matrix.py --minimum_ecdna 0 --protein_coding 0 --KNN_imputed_data_dir ./data/KNN_imputed_data --output ./data/RNAseq_expression_matrix 
```

**Generate training/testing datasets**
```
python -u ./src/Generate_testing_training_datasets.py --matrix ./data/RNAseq_expression_matrix/AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt
```

**Compute cliff delta for all TCGA genes in matrix**
```
python -u ./src/Compute_cliff_delta_all_TCGA_genes.py --matrix ./data/RNAseq_expression_matrix/AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt
```
```
python -u ./src/Compute_cliff_delta_all_TCGA_genes.py --matrix ./data/RNAseq_expression_matrix/AC049_TCGA_KNN_imputed_matrix.txt
```

**Generate read counts matrix**
```
python -u ./src/Generate_TCGA_level3_RSEM_raw_count_matrix.py --minimum_ecdna 3 --protein_coding 1 --output ./data/RNAseq_expression_matrix
```
```
python -u ./src/Generate_TCGA_level3_RSEM_raw_count_matrix.py --minimum_ecdna 0 --protein_coding 0 --output ./data/RNAseq_expression_matrix
```

**Compute DESeq2 LFC for all TCGA genes in matrix**
```
./src/Run_deseq2.R
```

**Determine gene directionality based on cliff's delta and DESeq2 LFC**
```
python -u ./src/Gene_directionality_cliffd_LFC.py --matrix ./data/RNAseq_expression_matrix/AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt
```
```
python -u ./src/Gene_directionality_cliffd_LFC.py --matrix ./data/RNAseq_expression_matrix/AC049_TCGA_KNN_imputed_matrix.txt
```

**Boruta analysis**
```
python -u ./src/Perform_Boruta_Trials.py --matrix ./data/RNAseq_expression_matrix/AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt
```

**Co-expression**
```
./src/Parse_pvclust_object.R
```
```
python -u ./src/Coexpressed_genes.py --boruta_dir ./data/AC049_TCGA_KNN_imputed_min_ecdna_3_matrix/Boruta_Trials
```

**Precision/recall CorEx genes** 
```
python -u ./src/CorEx_genes_precision_recall.py --matrix ./data/RNAseq_expression_matrix/AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt
```

**Geneset enrichment & clustering of genesets** 
```
python -u ./src/Perform_geneset_enrichment.py --matrix ./data/RNAseq_expression_matrix/AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt
```