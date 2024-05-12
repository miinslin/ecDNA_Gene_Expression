# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from sklearn.impute import KNNImputer
import pandas as pd
from map_geneid_symbol import map_gene_ids
 
'''
cBioPortal "data_RNA_Seq_v2_expression_median.txt" data contains batch corrected
values of the upper-quartile (UQ) normalized RSEM estimated counts data 
from Broad firehose. Due to the batch effect correction process, there are
missing values in the matrices that require preprocessing. This script uses the
K-nearest neighbors (KNN) method to replace missing values as described in
Troyanskaya et al., 2001 and implemented in scikit-learn's KNNImputer class.
'''

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-D", "--cBioPortal_data",
                    dest = "cbioportal",
                    default = os.path.join(base_dir, 'data', 'reference', 'TCGA', 'cBioportal'),
                    help="Directory containing cbioportal data where subdirectories are individual studies (e.g., brca_tcga_pan_can_atlas_2018).")
parser.add_argument("-O", "--output",
                    dest = "output",
                    default = os.path.join(base_dir, 'data', 'KNN_imputed_data'),
                    help="Output directory for processed cbioportal data.")
args = parser.parse_args()

cBioPortal_TCGA_data = os.path.normpath(args.cbioportal)
KNN_imputed_data_dir = os.path.normpath(args.output)
if not os.path.isdir(KNN_imputed_data_dir):
    print(KNN_imputed_data_dir, ' does not exist. Creating directory.')
    os.makedirs(KNN_imputed_data_dir)

# hugo name, gene id mappings
parse_hgnc = map_gene_ids(protein_coding_bool=True)
hgnc_entrez = parse_hgnc['hgnc_entrez']
print('hgnc_entrez: ', len(hgnc_entrez)) #19,218
hgnc_entrez['340602']['symbol'] = 'EZHIP'
hgnc_entrez['100134869'] = {'symbol':'UBE2Q2P2'}
hgnc_entrez['390284'] = {'symbol':'SRP14P1'}
hgnc_entrez['10357'] = {'symbol':'HMGB1P1'}
hgnc_entrez['653553']= {'symbol':'HSPB1P1'}

def updateHugoSymbol(geneName, geneId):
    if str(geneName)=='nan':
        if geneId in hgnc_entrez:
            HugoSymbol = hgnc_entrez[geneId]['symbol']
            return HugoSymbol
        else:
            return ''
    else:
        return geneName

def logTransform(exp_value):
    return np.log2(float(exp_value) + 1)

def backTransform(exp_value):
    return 2**exp_value

def createGeneID(geneName, geneId):
    return '|'.join([geneName, geneId])

def checkNanCells(matrix):
    #print('Documenting NA cells...')
    totalCells = matrix.shape[0]*matrix.shape[1]
    nonNanCells = matrix.stack().index.tolist()
    nanCells = np.argwhere(matrix.isnull().values).tolist()
    if len(nonNanCells)+len(nanCells)!=totalCells:
        raise ValueError('incorrect Nan cell determination...')
    return nanCells

def updateValuesifNan(origmatrix, imputedmatrix, nanPairs):
    print('updating ', len(nanPairs), ' cells with imputed values...')
    genes = origmatrix.index.tolist()
    samples = origmatrix.columns.tolist()
    for gene_i, sample_i in nanPairs:
        sample = samples[sample_i]
        gene = genes[gene_i]
        imputedValue = imputedmatrix[sample][gene]
        origmatrix.loc[gene, sample] = backTransform(imputedValue)
    return origmatrix

# Perform KNN imputation of missing values separately for each tumor type
def impute_missing_values(input_fn, output_fn):

    cBioPortal_data = pd.read_csv(input_fn, sep='\t', header=0, dtype='str')
    tumor_sampleIDs = list(cBioPortal_data.columns)[2:]
    tumor_samples_n = float(len(tumor_sampleIDs))
    cBioPortal_data[tumor_sampleIDs] = cBioPortal_data[tumor_sampleIDs].astype(float)

    #Replace SEPT15 typo with SEP15
    SEPT15_row = cBioPortal_data.index[cBioPortal_data['Hugo_Symbol']=='SEPT15'].tolist()
    if SEPT15_row:
        cBioPortal_data.replace('SEPT15', 'SEP15', inplace=True)

    #Update Hugo Symbol if not provided in cBioPortal file
    updated_Hugo_Symbols = cBioPortal_data.apply(lambda row: updateHugoSymbol(row.Hugo_Symbol, row.Entrez_Gene_Id), axis=1)
    cBioPortal_data['Hugo_Symbol'] = updated_Hugo_Symbols
    geneID_values = cBioPortal_data.apply(lambda row: createGeneID(row.Hugo_Symbol, row.Entrez_Gene_Id), axis=1)
    cBioPortal_data = cBioPortal_data.copy()
    cBioPortal_data.insert(0, "Gene_name|Gene_id", geneID_values, False)
    cBioPortal_data.drop(columns=['Hugo_Symbol', 'Entrez_Gene_Id'], inplace=True)

    #Remove genes where >60% of samples have missing values
    cBioPortal_data = cBioPortal_data.set_index('Gene_name|Gene_id')
    tumor_geneIDs = list(cBioPortal_data.index)
    genes_60p_nan = [ci for ci,c in enumerate(cBioPortal_data.isnull().sum(axis=1).tolist()) if float(c)/tumor_samples_n>0.6]
    tumor_genes_60_NA = [tumor_geneIDs[i] for i in genes_60p_nan]
    cBioPortal_data.drop(tumor_genes_60_NA, inplace=True)

    gene_list = list(cBioPortal_data.index)
    print(len(tumor_genes_60_NA), ' genes with >60% NA values')
    print('#samples: ', len(tumor_sampleIDs))
    print('#genes: ', len(gene_list))

    #Perform KNN imputation of missing values
    nanCells = checkNanCells(cBioPortal_data)
    if nanCells:
        #print('running KNNImputer...')
        lt_cBioPortal_data = cBioPortal_data.copy()
        lt_cBioPortal_data[tumor_sampleIDs] = lt_cBioPortal_data[tumor_sampleIDs].map(logTransform)
        imputer = KNNImputer(n_neighbors=5, weights="uniform", missing_values=np.nan)
        imputed_matrix = imputer.fit_transform(lt_cBioPortal_data)
        imputed_matrix = pd.DataFrame(imputed_matrix, columns=lt_cBioPortal_data.columns, index=lt_cBioPortal_data.index)
        # only use imputed values for nan values in original matrix
        cBioPortal_data = updateValuesifNan(cBioPortal_data, imputed_matrix, nanCells)
    cBioPortal_data.to_csv(output_fn, sep='\t', header=True, index=True)

    return tumor_sampleIDs

cBioportal_sample_to_tumor_type = {}
for dataset in os.listdir(cBioPortal_TCGA_data):
    print('\n',dataset)
    TCGA_code = dataset.split('_')[0].upper()
    data_RNA_Seq_v2 = os.path.join(cBioPortal_TCGA_data, dataset, 'data_RNA_Seq_v2_expression_median.txt')
    tumor_output_dir = os.path.join(KNN_imputed_data_dir, dataset)
    if not os.path.isdir(tumor_output_dir):
        os.mkdir(tumor_output_dir)
    tumor_knn_matrix_fn = os.path.join(tumor_output_dir, 'data_RNA_Seq_v2_expression_median_KNN_imputed.txt')
    cBioportal_sample_to_tumor_type[TCGA_code] = impute_missing_values(data_RNA_Seq_v2, tumor_knn_matrix_fn)


with open(os.path.join(base_dir, 'data', 'reference', 'TCGA','cBioportal_TCGA_sample_to_tumor_type.txt'),'w') as outfile:
    outfile.write('\t'.join(['tumor', 'sample'])+'\n')
    for tumor in cBioportal_sample_to_tumor_type:
        for sample in cBioportal_sample_to_tumor_type[tumor]:
            outfile.write('\t'.join([tumor, sample])+'\n')
    