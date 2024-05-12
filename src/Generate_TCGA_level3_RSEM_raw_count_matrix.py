# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from TCGA_sample_info import read_expression_matrix
from TCGA_sample_info import retrieve_sample_info

base_dir = os.path.abspath(os.path.join(__file__, "../.."))

parser = argparse.ArgumentParser()
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

TCGA_level3_preprocessed_dir = os.path.join(base_dir, 'data', 'reference', 'TCGA', 'Broad_Firehose',\
                                            'mRNAseq_Preprocess.Level_3', 'data')
RNAseq_expression_dir = os.path.normpath(args.output)

TCGA_knn_matrix_fn = os.path.join(RNAseq_expression_dir,\
                                'AC049_TCGA_KNN_imputed_matrix.txt')
TCGA_knn_matrix_samples_fn = os.path.join(RNAseq_expression_dir,\
                                'AC049_TCGA_KNN_imputed_matrix_samples.txt')
TCGA_matrix_fn = os.path.join(RNAseq_expression_dir,\
                                'AC049_TCGA_RSEM_estimated_counts_matrix.txt')
TCGA_matrix_samples_fn = os.path.join(RNAseq_expression_dir,\
                                'AC049_TCGA_RSEM_estimated_counts_matrix_samples.txt')
if min_ecdna_n!=0:
    TCGA_knn_matrix_fn = os.path.join(RNAseq_expression_dir,\
                                    'AC049_TCGA_KNN_imputed_min_ecdna_'+str(min_ecdna_n)+'_matrix.txt')
    TCGA_knn_matrix_samples_fn = os.path.join(RNAseq_expression_dir,\
                                    'AC049_TCGA_KNN_imputed_min_ecdna_'+str(min_ecdna_n)+'_matrix_samples.txt')
    TCGA_matrix_fn = os.path.join(RNAseq_expression_dir,\
                                    'AC049_TCGA_RSEM_estimated_counts_min_ecdna_'+str(min_ecdna_n)+'_matrix.txt')
    TCGA_matrix_samples_fn = os.path.join(RNAseq_expression_dir,\
                                    'AC049_TCGA_RSEM_estimated_counts_min_ecdna_'+str(min_ecdna_n)+'_matrix_samples.txt')

#####################################
#read gene expression matrix sample info
TCGA_info = retrieve_sample_info(min_ecdna = min_ecdna_n, remove_metastatic = True)
cbioportal_sample_info = TCGA_info['sample_info']
TCGA_tumors = TCGA_info['TCGA_tumors']
class_counts = TCGA_info['class_counts']

#read gene expression matrix
parse_matrix = read_expression_matrix(TCGA_knn_matrix_fn)
TCGA_combined_exp_matrix = parse_matrix['TCGA_combined_exp_matrix']
TCGA_KNN_matrix_samples = parse_matrix['TCGA_all_tumor_samples']

TCGA_KNN_matrix_genes = set()
TCGA_KNN_matrix_id_to_symbol = {}
for gene in TCGA_combined_exp_matrix.keys():
    genename, geneid = gene.split('|')
    if not genename:
        genename = '?'
    TCGA_KNN_matrix_id_to_symbol[geneid] = genename
    TCGA_KNN_matrix_genes.add(genename+'|'+geneid)
    
print('TCGA_KNN_matrix_genes: ', len(TCGA_KNN_matrix_genes))
print('TCGA_KNN_matrix_samples: ', len(TCGA_KNN_matrix_samples))

########### TCGA PanCancer Atlas Dataset  ###########
TCGA_all_tumor_samples = set()
TCGA_sample_to_tumor = {}
TCGA_gene_values = {}

all_ecDNA_samples = set()
all_nonecDNA_samples = set()
 
tumor_to_class = {}

'gdac.broadinstitute.org_BLCA.mRNAseq_Preprocess.Level_3.2016012800.0.0'
for dataset in os.listdir(TCGA_level3_preprocessed_dir):
    TCGA_code = dataset.split('_')[1].split('.')[0]
    
    if TCGA_code not in TCGA_tumors:
        print('skipping ', TCGA_code)
        continue
    
    #.uncv2.mRNAseq_raw_counts.txt
    rna_file = os.path.join(TCGA_level3_preprocessed_dir, dataset,\
                            TCGA_code+'.uncv2.mRNAseq_raw_counts.txt')
    tumor_to_class[TCGA_code] = {'ecDNA':set(), 'Non_ecDNA':set()}
    
    sampleIDs = []        
    header = ''
    colidx = ''
    with open(rna_file,'r') as infile:    
        for i,line in enumerate(infile):
            sp = line.rstrip().split('\t')
            if i ==0:
                header = line.rstrip().split('\t')
                #'HYBRIDIZATION R' 'TCGA-2F-A9KO-01'
                colidx = {x:xi for xi,x in enumerate(sp)}
                sampleIDs = [x for x in sp if x in TCGA_KNN_matrix_samples]                    
                if not sampleIDs:
                    raise ValueError('sampleIDs')
                    
                for x in sampleIDs:
                    TCGA_sample_to_tumor[x] = TCGA_code
                    AC_class = cbioportal_sample_info[x]['ecDNA_status']
                    if AC_class=='ecDNA':
                        all_ecDNA_samples.add(x)
                    if AC_class=='Non_ecDNA':
                        all_nonecDNA_samples.add(x)
                    tumor_to_class[TCGA_code][AC_class].add(x)
                TCGA_all_tumor_samples.update(sampleIDs)
            else:
                genename_geneid = sp[colidx['HYBRIDIZATION R']]
                genename, geneid = genename_geneid.split('|')
                if genename == '?':
                    if geneid in TCGA_KNN_matrix_id_to_symbol:
                        genename = TCGA_KNN_matrix_id_to_symbol[geneid]
                        if genename:
                            genename_geneid = genename+'|'+geneid
                
                if genename_geneid not in TCGA_KNN_matrix_genes:
                    continue                    
                
                if genename_geneid not in TCGA_gene_values:
                    TCGA_gene_values[genename_geneid] = {}
                    
                for sample in sampleIDs:
                    exp_value = sp[colidx[sample]].replace('"','').replace("'",'')
                    if exp_value == 'NA':
                        TCGA_gene_values[genename_geneid][sample] = {'float':np.nan, 'string':exp_value}
                    elif exp_value!='NA':
                        TCGA_gene_values[genename_geneid][sample] = {'float':float(exp_value), 'string':str(round(float(exp_value)))}

TCGA_all_tumor_samples = list(TCGA_all_tumor_samples) 
print('\nSamples with rnaseq data: ', len(TCGA_all_tumor_samples))
print('ecDNA+: ', len(all_ecDNA_samples))
print('ecDNA-: ', len(all_nonecDNA_samples))

TCGA_all_gene = sorted(TCGA_gene_values.keys())
print('\nGenes in matrix: ', str(len(TCGA_all_gene))) #16309

#only include genes with no NA values
TCGA_all_gene = sorted(TCGA_all_gene)
print(str(len(TCGA_all_gene))+' TCGA_all_genes')

############## write to file  ##############     
ecDNA_samples = sorted(['ecDNA'+'|'+TCGA_sample_to_tumor[x]+'|'+x for x in all_ecDNA_samples]) #
non_ecDNA_samples = sorted(['Non_ecDNA'+'|'+TCGA_sample_to_tumor[x]+'|'+x for x in all_nonecDNA_samples]) #
 
with open(TCGA_matrix_samples_fn,'w') as sampleinfo:
    sampleinfo.write('\t'.join(['#ecDNA_status','Tumor', 'Sample_barcode', 'SampleType',\
                                'full_barcode', 'Center','platform'])+'\n')
    for sample in ecDNA_samples:
        ecDNA_status, tumor, samplebarcode = sample.split('|')
        sample_type = cbioportal_sample_info[samplebarcode]['SampleType']
        full_barcode = cbioportal_sample_info[samplebarcode]['full_barcode']
        Center = cbioportal_sample_info[samplebarcode]['Center']
        platform = cbioportal_sample_info[samplebarcode]['platform']
        sampleinfo.write('\t'.join([ecDNA_status, tumor, samplebarcode,\
                                    sample_type, full_barcode,\
                                    Center, platform])+'\n')
    for sample in non_ecDNA_samples:
        ecDNA_status, tumor, samplebarcode = sample.split('|')
        sample_type = cbioportal_sample_info[samplebarcode]['SampleType']
        full_barcode = cbioportal_sample_info[samplebarcode]['full_barcode']
        Center = cbioportal_sample_info[samplebarcode]['Center']
        platform = cbioportal_sample_info[samplebarcode]['platform']
        sampleinfo.write('\t'.join([ecDNA_status, tumor, samplebarcode,\
                                    sample_type, full_barcode,\
                                    Center, platform])+'\n')

with open(TCGA_matrix_fn,'w') as outfile:
    outfile.write('Gene_name|Gene_id'+'\t'+'\t'.join(ecDNA_samples)+'\t'+'\t'.join(non_ecDNA_samples)+'\n')
    for genename_geneid in TCGA_all_gene:
        values = []
        for entry in ecDNA_samples:
            sample = entry.split('|')[2]
            values.append(TCGA_gene_values[genename_geneid][sample]['string'])
        for entry in non_ecDNA_samples:
            sample = entry.split('|')[2]
            values.append(TCGA_gene_values[genename_geneid][sample]['string'])
        outfile.write(genename_geneid+'\t'+'\t'.join(values)+'\n')

