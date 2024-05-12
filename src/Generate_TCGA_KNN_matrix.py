# -*- coding: utf-8 -*-
import os
import argparse
from map_geneid_symbol import map_gene_ids
from TCGA_sample_info import retrieve_sample_info
import pandas as pd
import numpy as np

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-D", "--KNN_imputed_data_dir", dest = "data", default = os.path.join(base_dir, 'data', 'KNN_imputed_data'), \
                    help="Output directory for processed cbioportal data generated via cBioPortal_data_preprocessing.py.")
parser.add_argument("-e", "--minimum_ecdna", dest = "min_ecdna", default = 0,\
                    help="Minimum number of ecDNA(+) samples required per tumor type.")
parser.add_argument("-p", "--protein_coding", dest = "prot_coding", default = 0,\
                    help="Include only protein-coding genes in matrix. 0:False; 1:True.")
parser.add_argument("-O", "--output", dest = "output", default = os.path.join(base_dir, 'data', 'RNAseq_expression_matrix'), \
                    help="Output directory for gene expression matrices.")    
args = parser.parse_args()

min_ecdna_n = int(args.min_ecdna)
protein_coding = False
if int(args.prot_coding)==1:
    protein_coding = True

print('min_ecdna_n: ', min_ecdna_n)
print('protein_coding: ', protein_coding)

KNN_imputed_data_dir = os.path.normpath(args.data)

RNAseq_expression_dir = os.path.normpath(args.output)
if not os.path.isdir(RNAseq_expression_dir):
    print('Creating ', RNAseq_expression_dir)
    os.makedirs(RNAseq_expression_dir)

TCGA_knn_matrix_fn = os.path.join(RNAseq_expression_dir,\
                                'AC049_TCGA_KNN_imputed_matrix.txt')
TCGA_knn_matrix_samples_fn = os.path.join(RNAseq_expression_dir,\
                                'AC049_TCGA_KNN_imputed_matrix_samples.txt')
if min_ecdna_n!=0:
    TCGA_knn_matrix_fn = os.path.join(RNAseq_expression_dir,\
                                    'AC049_TCGA_KNN_imputed_min_ecdna_'+str(min_ecdna_n)+'_matrix.txt')
    TCGA_knn_matrix_samples_fn = os.path.join(RNAseq_expression_dir,\
                                    'AC049_TCGA_KNN_imputed_min_ecdna_'+str(min_ecdna_n)+'_matrix_samples.txt')

TCGA_info = retrieve_sample_info(min_ecdna = min_ecdna_n, remove_metastatic = True) 
parse_hgnc = map_gene_ids(protein_coding_bool = protein_coding)
hgnc_entrez = parse_hgnc['hgnc_entrez']
hgnc_symbol = parse_hgnc['hgnc_symbol']

sample_info = TCGA_info['sample_info']
all_samples = set(sample_info.keys())
all_ecDNA_samples = TCGA_info['class_counts']['ecDNA']
all_nonecDNA_samples = TCGA_info['class_counts']['Non_ecDNA']
TCGA_tumors = TCGA_info['TCGA_tumors']

#samples: ecDNA_status|Tumor|Sample_barcode
ecDNA_samples = sorted(['ecDNA'+'|'+sample_info[x]['Tumor']+'|'+x for x in all_ecDNA_samples])
non_ecDNA_samples = sorted(['Non_ecDNA'+'|'+sample_info[x]['Tumor']+'|'+x for x in all_nonecDNA_samples])

#write TCGA_knn_matrix_samples_fn
with open(TCGA_knn_matrix_samples_fn,'w') as sampleinfo:
    sampleinfo.write('\t'.join(['#ecDNA_status','Tumor', 'Sample_barcode',\
                                'SampleType', 'full_barcode',\
                                'Center','platform'])+'\n')
    for sample in ecDNA_samples:
        ecDNA_status, tumor, samplebarcode = sample.split('|')
        sample_type = sample_info[samplebarcode]['SampleType']
        full_barcode = sample_info[samplebarcode]['full_barcode']
        Center = sample_info[samplebarcode]['Center']
        platform = sample_info[samplebarcode]['platform']
        sampleinfo.write('\t'.join([ecDNA_status, tumor, samplebarcode,\
                                    sample_type, full_barcode,\
                                    Center, platform])+'\n')
    for sample in non_ecDNA_samples:
        ecDNA_status, tumor, samplebarcode = sample.split('|')
        sample_type = sample_info[samplebarcode]['SampleType']
        full_barcode = sample_info[samplebarcode]['full_barcode']
        Center = sample_info[samplebarcode]['Center']
        platform = sample_info[samplebarcode]['platform']
        sampleinfo.write('\t'.join([ecDNA_status, tumor, samplebarcode,\
                                    sample_type, full_barcode,\
                                    Center, platform])+'\n')

#read KNN imputed matrices
sharedGenes = []
for dataset in os.listdir(KNN_imputed_data_dir):
    TCGA_code = dataset.split('_')[0].upper()
    if TCGA_code not in TCGA_tumors:
        continue
    KNN_RNA_Seq_fn = os.path.join(KNN_imputed_data_dir, dataset, \
                                  'data_RNA_Seq_v2_expression_median_KNN_imputed.txt')
    KNN_RNA_Seq_data = pd.read_csv(KNN_RNA_Seq_fn, sep='\t', header=0, index_col=0, dtype='str')
    sharedGenes.append(set(KNN_RNA_Seq_data.index))
sharedGenes = set.intersection(*map(set,sharedGenes))
#print('sharedGenes: ', len(sharedGenes))

#
matrix_genes = set()
if protein_coding == False:
    matrix_genes.update(sharedGenes)
else:
    for genename_geneid in sharedGenes:
        genename, gene_id = genename_geneid.split('|')
        skip = False
        if gene_id not in hgnc_entrez:
            skip = True
            if genename in hgnc_symbol:
                skip = False
            if skip==True:
                continue
        matrix_genes.add(genename_geneid)
print('matrix_genes: ', len(matrix_genes))


TCGA_gene_values = {g:{s:'' for s in all_samples} for g in matrix_genes}
for dataset in os.listdir(KNN_imputed_data_dir):
    TCGA_code = dataset.split('_')[0].upper()
    if TCGA_code not in TCGA_tumors:
        continue
    KNN_RNA_Seq_fn = os.path.join(KNN_imputed_data_dir, dataset,\
                                    'data_RNA_Seq_v2_expression_median_KNN_imputed.txt')
    KNN_RNA_Seq_data = pd.read_csv(KNN_RNA_Seq_fn, sep='\t', header=0, index_col=0, dtype='str')
    KNN_RNA_Seq_data = KNN_RNA_Seq_data[KNN_RNA_Seq_data.columns.intersection(all_samples)]
    #print(TCGA_code, KNN_RNA_Seq_data.shape)
    for gene in matrix_genes:
        for sample in KNN_RNA_Seq_data.columns:
            value = float(KNN_RNA_Seq_data[sample][gene])
            if np.isnan(value):
                raise ValueError('np.nan value')
            TCGA_gene_values[gene][sample] = value

#write TCGA_knn_matrix_fn
with open(TCGA_knn_matrix_fn,'w') as outfile:
    outfile.write('Gene_name|Gene_id'+'\t'+'\t'.join(ecDNA_samples)+'\t'+'\t'.join(non_ecDNA_samples)+'\n')
    for genename_geneid in sorted(matrix_genes):
        values = []
        for entry in ecDNA_samples:
            sample = entry.split('|')[2]
            values.append(TCGA_gene_values[genename_geneid][sample])
        for entry in non_ecDNA_samples:
            sample = entry.split('|')[2]
            values.append(TCGA_gene_values[genename_geneid][sample])
        outfile.write(genename_geneid+'\t'+'\t'.join([str(x) for x in values])+'\n')