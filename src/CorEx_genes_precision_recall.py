# -*- coding: utf-8 -*-
import os
import argparse
import pandas as pd
from TCGA_sample_info import read_expression_matrix
from gene_validation_precision_recall import validate_genes

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--matrix", dest = "matrix",\
                    default = os.path.join(base_dir, 'data', \
                    'RNAseq_expression_matrix',\
                    'AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt'),\
                    help="Gene expression matrix path.")
args = parser.parse_args()

TCGA_matrix_fn = os.path.normpath(args.matrix)
matrix_type = os.path.basename(TCGA_matrix_fn).replace('.txt','')
training_testing_dir = os.path.join(base_dir, 'data', matrix_type, 'Training_Testing_datasets')

num_runs = 200
percentage_of_runs = 0.05
min_runs = int(num_runs*percentage_of_runs)

au_threshold = 95
coexpression_dir = os.path.join(base_dir, 'data', matrix_type, 'co_expression')
corex_gene_summary = os.path.join(coexpression_dir, 'au'+str(au_threshold)+'_cluster_boruta_genes.txt')

corex_gene_summary_df = pd.read_csv(corex_gene_summary, sep='\t', header=0, dtype='str')
corex_gene_summary_df = corex_gene_summary_df.astype({'#Trials(200)_by_cluster':float})
CorEx_genes = corex_gene_summary_df.loc[corex_gene_summary_df['#Trials(200)_by_cluster']>=min_runs]['Gene|GeneId'].tolist()
CorEx_genes = sorted(CorEx_genes)
print('CorEx_genes: ', len(CorEx_genes))

model_performance_dir = os.path.join(base_dir, 'data', matrix_type,\
                              'Model_performance',\
                              'CorEx_'+str(len(CorEx_genes))+'_genes')
if not os.path.isdir(model_performance_dir):
    os.makedirs(model_performance_dir)
        
parse_matrix = read_expression_matrix(TCGA_matrix_fn)
TCGA_combined_exp_matrix = parse_matrix['TCGA_combined_exp_matrix']
TCGA_all_tumor_samples = parse_matrix['TCGA_all_tumor_samples']

TCGA_all_tumor_samples = sorted(TCGA_all_tumor_samples)
print('#samples: ', len(TCGA_all_tumor_samples))

TCGA_data_sub = {}
for sample in TCGA_all_tumor_samples:
    values = [TCGA_combined_exp_matrix[g][sample] for g in CorEx_genes]
    TCGA_data_sub[sample] = values

validate_genes(training_testing_dir, model_performance_dir, TCGA_data_sub, CorEx_genes, num_runs, '', '', '', '')