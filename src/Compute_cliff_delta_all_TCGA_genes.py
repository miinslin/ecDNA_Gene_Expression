# -*- coding: utf-8 -*-
import os
from scipy.stats import mannwhitneyu
import pickle
from cliffs_delta import cliffs_delta
import argparse
from TCGA_sample_info import read_expression_matrix

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

alpha = 0.05
 
#read gene expression matrix
parse_matrix = read_expression_matrix(TCGA_matrix_fn)
TCGA_combined_exp_matrix = parse_matrix['TCGA_combined_exp_matrix']
TCGA_all_tumor_samples = parse_matrix['TCGA_all_tumor_samples']
sample_classification = parse_matrix['sample_classification']
TCGA_sample_to_tumor = parse_matrix['TCGA_sample_to_tumor']
TCGA_sample_to_class = parse_matrix['TCGA_sample_to_class']

TCGA_all_gene = sorted(TCGA_combined_exp_matrix.keys())
print('#genes: ', len(TCGA_all_gene))
TCGA_all_tumor_samples = sorted(TCGA_all_tumor_samples)
print('#samples: ', len(TCGA_all_tumor_samples)) 
TCGA_tumors = sorted(set(TCGA_sample_to_tumor.values()))
print('#tumors: ', len(TCGA_tumors))

all_ecDNA_samples = sorted(sample_classification['ecDNA'])
all_nonecDNA_samples = sorted(sample_classification['Non_ecDNA'])
print('ecDNA(+): ', len(all_ecDNA_samples))
print('ecDNA(-): ', len(all_nonecDNA_samples))

TCGA_tumor_samples = {tumor:{'ecDNA':set(),'Non_ecDNA':set()} for tumor in TCGA_tumors}
for sample_barcode in TCGA_all_tumor_samples:
    TCGA_tumor_samples[TCGA_sample_to_tumor[sample_barcode]][TCGA_sample_to_class[sample_barcode]].add(sample_barcode)

#output files
Cliff_delta_dir = os.path.join(base_dir, 'data', matrix_type, 'Cliff_Delta')
if not os.path.isdir(Cliff_delta_dir):
    os.makedirs(Cliff_delta_dir)
    
MWU_pkl = os.path.join(Cliff_delta_dir, 'All_TCGA_genes_cliffd.pkl')
results_fn = os.path.join(Cliff_delta_dir, 'All_TCGA_genes_cliffd.txt')

Tumor_MWU_pkl = os.path.join(Cliff_delta_dir, 'TCGA_tumor_genes_cliffd.pkl')
Tumor_results_fn = os.path.join(Cliff_delta_dir, 'TCGA_tumor_genes_cliffd.txt')

 
'''
d : -1 ≤ δ ≤ 1. 
measures how often values in one distribution are larger than values in
a second distribution. non-parametric.
Values near ±1 indicate the absence of overlap between the two samples,
while values near zero indicate a lot of overlap between the two samples.

The magnitude is assessed using the thresholds provided in (Romano 2006),
|d|<0.147 "negligible" 
|d|<0.33 "small"
|d|<0.474 "medium" 
otherwise "large"
'''

def sign_num(num):
    x = ''
    if abs(num)/num==1:
        x = 'positive'
    elif abs(num)/num==-1:
        x = 'negative'
    else:
        raise ValueError()
    return x

    
results_ecDNA_up = {}
results_ecDNA_down = {}
mannwhitneyu_list = set()
MWU_all = {}
for gg in TCGA_all_gene:
    
    ecDNA_gene_values = [TCGA_combined_exp_matrix[gg][sample] for sample in all_ecDNA_samples]
    noecDNA_gene_values = [TCGA_combined_exp_matrix[gg][sample] for sample in all_nonecDNA_samples]
    
    nx, ny = len(ecDNA_gene_values), len(noecDNA_gene_values)
    
    cliff_d, res = cliffs_delta(ecDNA_gene_values, noecDNA_gene_values)
    
    MWU_direction = ''
    
    #ecDNA up, noecDNA down
    U1_up, p_up = mannwhitneyu(ecDNA_gene_values, noecDNA_gene_values, alternative='greater')
    U2_up = nx*ny - U1_up
    if p_up<=alpha:
        results_ecDNA_up[gg] = {'U1':U1_up, 'U2':U2_up, 'p':p_up, 'cliff_d':cliff_d} 
        mannwhitneyu_list.add(gg)
        MWU_all[gg] = {'U1':U1_up, 'U2':U2_up, 'p':p_up, 'direction':'UP', 'cliff_d':cliff_d} 
        MWU_direction = 'UP'

    #ecDNA down, noecDNA up
    U1_down, p_down = mannwhitneyu(ecDNA_gene_values, noecDNA_gene_values, alternative='less')
    U2_down = nx*ny - U1_down    
    if p_down<=alpha:
        results_ecDNA_down[gg] = {'U1':U1_down, 'U2':U2_down, 'p':p_down, 'cliff_d':cliff_d} 
        mannwhitneyu_list.add(gg)
        if gg in MWU_all:
            print(gg, 'UP & DOWN', p_up, p_down)
            raise ValueError()
        MWU_all[gg] = {'U1':U1_down, 'U2':U2_down, 'p':p_down, 'direction':'DOWN', 'cliff_d':cliff_d} 
        MWU_direction = 'DOWN'

    #if both MWU tests no-significant, direction of non-sig is smaller p-value
    if gg not in MWU_all:
        #if p_down smaller
        if p_down<p_up:
            MWU_all[gg] = {'U1':U1_down, 'U2':U2_down, 'p':p_down, 'direction':'DOWN', 'cliff_d':cliff_d}
        #if p_up smaller
        elif p_up<p_down:
            MWU_all[gg] = {'U1':U1_up, 'U2':U2_up, 'p':p_up, 'direction':'UP', 'cliff_d':cliff_d} 
        #tie? keep both
        else:
            MWU_all[gg] = {'U1':'', 'U2':'', 'p':p_up, 'direction':'NA', 'cliff_d':cliff_d}  
    
outfile = open(MWU_pkl,'wb')
pickle.dump(MWU_all, outfile)
outfile.close()

print('\nMann-Whitney U, up regulated in ecDNA (down regulated in nonecDNA): ', len(results_ecDNA_up))
print('Mann-Whitney U, down regulated in ecDNA (up regulated in nonecDNA): ', len(results_ecDNA_down))
print('Mann-Whitney U all: ', len(mannwhitneyu_list)) #10,848

with open(results_fn,'w') as outfile:
    outfile.write('\t'.join(['#Gene_name|Gene_id','Mann-Whitney_U','Mann-Whitney_U_pvalue','cliff_d'])+'\n')
    for gg in TCGA_all_gene:
        
        #mannwhitneyu (smaller value, higher rank)
        mannwhitu = ''
        mannwhitu_p = ''
        cliff_d = ''
        cliff_d_sign = ''
 
        if gg in results_ecDNA_up:

            mannwhitu = 'upregulated(ecDNA)'
            mannwhitu_p = results_ecDNA_up[gg]['p']
            cliff_d = results_ecDNA_up[gg]['cliff_d']
            
            cliff_d_sign = sign_num(cliff_d)
            if cliff_d_sign!='positive':
                raise ValueError('up', cliff_d)
            
        elif gg in results_ecDNA_down:

            mannwhitu = 'downregulated(ecDNA)'
            mannwhitu_p = results_ecDNA_down[gg]['p']
            cliff_d = results_ecDNA_down[gg]['cliff_d']            

            cliff_d_sign = sign_num(cliff_d)
            if cliff_d_sign!='negative':
                raise ValueError('down', cliff_d)
                
        else:
            mannwhitu = 'NotSignificant|'+MWU_all[gg]['direction']
            mannwhitu_p = MWU_all[gg]['p']
            cliff_d = MWU_all[gg]['cliff_d']
   
        outfile.write('\t'.join([gg, mannwhitu, repr(mannwhitu_p), repr(cliff_d)])+'\n')


############################## per tumor #######################################

MWU_all = {tumor:{} for tumor in TCGA_tumors}
with open(Tumor_results_fn,'w') as tumor_results_out:
    tumor_results_out.write('\t'.join(['#tumor', 'Gene_name|Gene_id', 'Mann-Whitney_U','Mann-Whitney_U_pvalue','cliff_d'])+'\n')
    
    for tumor in TCGA_tumors:
        print('\n')
        print(tumor)
        
        results_ecDNA_up = {}
        results_ecDNA_down = {}
        mannwhitneyu_list = set()
        
        tumor_genes_parsed = set()
        
        for gg in TCGA_all_gene:
        
            ecDNA_gene_values = [TCGA_combined_exp_matrix[gg][sample] for sample in TCGA_tumor_samples[tumor]['ecDNA']]
            noecDNA_gene_values = [TCGA_combined_exp_matrix[gg][sample] for sample in TCGA_tumor_samples[tumor]['Non_ecDNA']]
                        
            nx, ny = len(ecDNA_gene_values), len(noecDNA_gene_values)

            #if less than 8 ecdna or non-ecdna samples, should not run mwu
            if nx<8 or ny<8:
                continue
            
            tumor_genes_parsed.add(gg)
            
            cliff_d, res = cliffs_delta(ecDNA_gene_values, noecDNA_gene_values)
            
            #ecDNA up, noecDNA down
            U1_up, p_up = mannwhitneyu(ecDNA_gene_values, noecDNA_gene_values, alternative='greater')
            U2_up = nx*ny - U1_up
            if p_up<=alpha:
                results_ecDNA_up[gg] = {'U1':U1_up, 'U2':U2_up, 'p':p_up, 'cliff_d':cliff_d} 
                mannwhitneyu_list.add(gg)
                MWU_all[tumor][gg] = {'U1':U1_up, 'U2':U2_up, 'p':p_up, 'direction':'UP', 'cliff_d':cliff_d} 
        
            #ecDNA down, noecDNA up
            U1_down, p_down = mannwhitneyu(ecDNA_gene_values, noecDNA_gene_values, alternative='less')
            U2_down = nx*ny - U1_down    
            if p_down<=alpha:
                results_ecDNA_down[gg] = {'U1':U1_down, 'U2':U2_down, 'p':p_down, 'cliff_d':cliff_d} 
                mannwhitneyu_list.add(gg)
                if gg in MWU_all[tumor]:
                    print(gg, 'UP & DOWN', p_up, p_down)
                    raise ValueError()
                MWU_all[tumor][gg] = {'U1':U1_down, 'U2':U2_down, 'p':p_down, 'direction':'DOWN', 'cliff_d':cliff_d} 
            
            #if both MWU tests no-significant, direction of non-sig is smaller p-value
            if gg not in MWU_all[tumor]:
                #if p_down smaller
                if p_down<p_up:
                    MWU_all[tumor][gg] = {'U1':U1_down, 'U2':U2_down, 'p':p_down, 'direction':'DOWN', 'cliff_d':cliff_d}
                #if p_up smaller
                elif p_up<p_down:
                    MWU_all[tumor][gg] = {'U1':U1_up, 'U2':U2_up, 'p':p_up, 'direction':'UP', 'cliff_d':cliff_d} 
                #tie? keep both
                else:
                    MWU_all[tumor][gg] = {'U1':'', 'U2':'', 'p':p_up, 'direction':'NA', 'cliff_d':cliff_d}  
        
        tumor_genes_parsed = sorted(tumor_genes_parsed)
        
        print('tumor_genes_parsed: ', len(tumor_genes_parsed))
        
        if tumor_genes_parsed:
            print('\nMann-Whitney U, up regulated in ecDNA (down regulated in nonecDNA): ', len(results_ecDNA_up))
            print('Mann-Whitney U, down regulated in ecDNA (up regulated in nonecDNA): ', len(results_ecDNA_down))
            print('Mann-Whitney U all: ', len(mannwhitneyu_list))
            
        
            for gg in tumor_genes_parsed:
                #mannwhitneyu (smaller value, higher rank)
                mannwhitu = ''
                mannwhitu_p = ''
                cliff_d = ''
                cliff_d_sign = ''
         
                if gg in results_ecDNA_up:
                    mannwhitu = 'upregulated(ecDNA)'
                    mannwhitu_p = results_ecDNA_up[gg]['p']
                    cliff_d = results_ecDNA_up[gg]['cliff_d']

                    cliff_d_sign = sign_num(cliff_d)
                    if cliff_d_sign!='positive':
                        raise ValueError('up', cliff_d)                    

                elif gg in results_ecDNA_down:
                    mannwhitu = 'downregulated(ecDNA)'
                    mannwhitu_p = results_ecDNA_down[gg]['p']
                    cliff_d = results_ecDNA_down[gg]['cliff_d']

                    cliff_d_sign = sign_num(cliff_d)
                    if cliff_d_sign!='negative':
                        raise ValueError('down', cliff_d)

                else:
                    mannwhitu = 'NotSignificant|'+MWU_all[tumor][gg]['direction']
                    mannwhitu_p = MWU_all[tumor][gg]['p']
                    cliff_d = MWU_all[tumor][gg]['cliff_d']
           
                tumor_results_out.write('\t'.join([tumor, gg, mannwhitu, repr(mannwhitu_p), repr(cliff_d)])+'\n')


outfile = open(Tumor_MWU_pkl,'wb')
pickle.dump(MWU_all, outfile)
outfile.close()