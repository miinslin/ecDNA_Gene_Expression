suppressPackageStartupMessages(library("DESeq2"))

#matrixType = "AC049_TCGA_KNN_imputed_matrix" 
matrixType = "AC049_TCGA_KNN_imputed_min_ecdna_3_matrix" 

basedir <- dirname(sys.frame(1)$ofile)
basedir = dirname(basedir)
print(basedir)

matrix_dir = file.path(basedir, "data", "RNAseq_expression_matrix")


if (matrixType=="AC049_TCGA_KNN_imputed_matrix") {
  AC_0_4_9_TCGA_matrix_fn = file.path(matrix_dir, "AC049_TCGA_RSEM_estimated_counts_matrix.txt")
  AC_0_4_9_TCGA_matrix_samples_fn = file.path(matrix_dir, "AC049_TCGA_RSEM_estimated_counts_matrix_samples.txt")
  DeseqPrefix = "DESeq_DE"
} else if (matrixType=="AC049_TCGA_KNN_imputed_min_ecdna_3_matrix") {
  AC_0_4_9_TCGA_matrix_fn = file.path(matrix_dir, "AC049_TCGA_RSEM_estimated_counts_min_ecdna_3_matrix.txt")
  AC_0_4_9_TCGA_matrix_samples_fn = file.path(matrix_dir, "AC049_TCGA_RSEM_estimated_counts_min_ecdna_3_matrix_samples.txt")
  DeseqPrefix = "DESeq_DE_ecdna3"
}


rsem_counts_dir = file.path(basedir, "data", "RSEM_estimated_counts", DeseqPrefix, "TCGA_all")
dir.create(file.path(rsem_counts_dir), recursive = TRUE)


deseq_batch_corr_all_genes = file.path(rsem_counts_dir, "condition_ecDNA_vs_Non_ecDNA_results_batch_corr_ALL.tsv") 

########## READ matrix ########## 
AC_0_4_9_TCGA_matrix <- read.csv(AC_0_4_9_TCGA_matrix_fn, header= TRUE, sep='\t',check.names=FALSE)
print(dim(AC_0_4_9_TCGA_matrix))

AC_0_4_9_TCGA_matrix_samples <- read.csv(AC_0_4_9_TCGA_matrix_samples_fn, header= TRUE, sep='\t',check.names=FALSE)
AC_0_4_9_TCGA_matrix_samples = AC_0_4_9_TCGA_matrix_samples[AC_0_4_9_TCGA_matrix_samples$SampleType != "Metastatic", ]  

ecdna_status = AC_0_4_9_TCGA_matrix_samples$"#ecDNA_status"
tumor = AC_0_4_9_TCGA_matrix_samples$Tumor
samples = AC_0_4_9_TCGA_matrix_samples$Sample_barcode
non_metastatic_samples = c()
for (i in 1:length(samples)) {
  non_metastatic_samples[i] = paste(ecdna_status[i],tumor[i],samples[i],sep='|')
}

rownames(AC_0_4_9_TCGA_matrix_samples) <- AC_0_4_9_TCGA_matrix_samples$Sample_barcode
rownames(AC_0_4_9_TCGA_matrix) <- AC_0_4_9_TCGA_matrix$"Gene_name|Gene_id"

AC_0_4_9_TCGA_matrix = AC_0_4_9_TCGA_matrix[,(names(AC_0_4_9_TCGA_matrix) %in% non_metastatic_samples)]

current_column_names = colnames(AC_0_4_9_TCGA_matrix)
ecDNA_status_list = c()
tumor_list = c()
samplebarcode_list = c()
new_column_names = c()
batch_list = c()
center_list = c()
platform_list = c()
for (i in 1:length(current_column_names)) {
  entry = current_column_names[i]
  #"Non_ecDNA|HNSC|TCGA-CR-7391-01"
  ecDNA_status = strsplit(entry, split = "\\|")[[1]][1]
  ecDNA_status_list[i] = ecDNA_status
  tumor = strsplit(entry, split = "\\|")[[1]][2]
  tumor_list[i] = tumor
  
  samplebarcode = strsplit(entry, split = "\\|")[[1]][3]
  samplebarcode_list[i] = samplebarcode
  new_column_names[i] = samplebarcode
  
  center_list[i] = AC_0_4_9_TCGA_matrix_samples[samplebarcode,]$Center
  platform_list[i] = AC_0_4_9_TCGA_matrix_samples[samplebarcode,]$platform 
  batch_list[i] = paste(AC_0_4_9_TCGA_matrix_samples[samplebarcode,]$Center, AC_0_4_9_TCGA_matrix_samples[samplebarcode,]$platform, sep="_") 
}

colnames(AC_0_4_9_TCGA_matrix) <- new_column_names
AC_0_4_9_TCGA_matrix <- as.matrix(AC_0_4_9_TCGA_matrix)

coldata <-data.frame("condition" = ecDNA_status_list,
                     "tumor" = tumor_list,
                     "batch" = batch_list,
                     "center" = center_list,
                     "platform" = platform_list,
                     "sample_barcode" = samplebarcode_list)
rownames(coldata) = coldata$sample_barcode 

coldata$condition <- factor(coldata$condition)
coldata$tumor <- factor(coldata$tumor)
coldata$batch <- factor(coldata$batch)

all(rownames(coldata) %in% colnames(AC_0_4_9_TCGA_matrix))
all(rownames(coldata) == colnames(AC_0_4_9_TCGA_matrix))

dds_batch_corr <- DESeqDataSetFromMatrix(countData = AC_0_4_9_TCGA_matrix,
                                         colData = coldata,
                                         design = ~batch+condition)
keep <- rowSums(counts(dds_batch_corr)) >= 10
dds_batch_corr <- dds_batch_corr[keep,]

dds_batch_corr$condition <- relevel(dds_batch_corr$condition, ref = "Non_ecDNA")

dds_batch_corr <- DESeq(dds_batch_corr)
resultsNames(dds_batch_corr)

lfc.cutoff <- log2(1.1)


## results ##
res_batch_corr <- results(dds_batch_corr, lfcThreshold=lfc.cutoff, contrast=c("condition","ecDNA","Non_ecDNA"))
summary(res_batch_corr)

## Log fold change shrinkage results ## 
resLFC_batch_corr <- lfcShrink(dds_batch_corr, lfcThreshold=lfc.cutoff, coef="condition_ecDNA_vs_Non_ecDNA", type="apeglm")
summary(resLFC_batch_corr)

#add shrinkage post lfc
res_batch_corr$posteriorLFC <- resLFC_batch_corr$log2FoldChange
res_batch_corr$svalue <- resLFC_batch_corr$svalue
write.table(as.data.frame(res_batch_corr), sep = "\t", file=deseq_batch_corr_all_genes, quote=FALSE)