# -*- coding: utf-8 -*-
import os
import numpy as np
import argparse

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

cliff_d_threshold = 0.147
post_lfc_threshold = np.log2(1.1)
 
#MWU cliff's delta (non-parametric)
Cliff_delta_dir = os.path.join(base_dir, 'data', matrix_type, 'Cliff_Delta')
MWU_all_fn = os.path.join(Cliff_delta_dir, 'All_TCGA_genes_cliffd.txt')
#posterior LFC effect size (negative binomial)
#output from running DeSeq2 (design = ~batch+condition)
if matrix_type=='AC049_TCGA_KNN_imputed_min_ecdna_3_matrix':
    DESeq2_dir = os.path.join(base_dir, 'data', 'RSEM_estimated_counts', 'DESeq_DE_ecdna3')
if matrix_type=='AC049_TCGA_KNN_imputed_matrix':
    DESeq2_dir = os.path.join(base_dir, 'data', 'RSEM_estimated_counts', 'DESeq_DE')
deseq_output_fn = os.path.join(DESeq2_dir, 'TCGA_all', \
                               'condition_ecDNA_vs_Non_ecDNA_results_batch_corr_ALL.tsv')

gene_directionality_dir = os.path.join(base_dir, 'data', matrix_type,\
                                       'gene_directionality')
if not os.path.isdir(gene_directionality_dir):
    os.makedirs(gene_directionality_dir)

def cliff_delta_magnitude(n):
    magnitude = ''
    if abs(n)<0.33:
        magnitude = 1
    elif 0.33<=abs(n)<0.474:
        magnitude = 2
    elif abs(n)>=0.474:
        magnitude = 3        
    return magnitude

def deseq_lfc_magnitude(n):
    magnitude = ''
    if abs(n)<np.log2(1.5):
        magnitude = 1
    elif np.log2(1.5)<=abs(n)<np.log2(2):
        magnitude = 2
    elif abs(n)>=np.log2(2):
        magnitude = 3        
    return magnitude

def sign_num(num):
    if num==0:
        return('0')
    if abs(num)/num==1:
        return('UP')
    elif abs(num)/num==-1:
        return('DOWN')

entrez_id_to_full = {}
non_passing_genes = {}
gene_direction = {}
gene_effect_size = {}

################ cliff's delta effect size ################
cliff_delta_results = {}
MWU = {}
with open(MWU_all_fn,'r') as infile:
    for i,line in enumerate(infile):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            gg = sp[colidx['#Gene_name|Gene_id']]
            
            symbol, entrez = gg.split('|')
            entrez_id_to_full[entrez] = gg
            
            MWU_pvalue = sp[colidx['Mann-Whitney_U_pvalue']]
            MWU_direction = sp[colidx['Mann-Whitney_U']]
            MWU[gg] = {'pvalue':MWU_pvalue,
                       'direction':MWU_direction}
            
            cliff_d = float(sp[colidx['cliff_d']])
            direction = sign_num(cliff_d)
            gene_effect_size[gg] = {'cliff_d':cliff_d}
            
            if abs(cliff_d)<cliff_d_threshold:
                non_passing_genes[gg] = {'cliff_d':cliff_d}
                continue
            
            cliff_delta_results[gg] = direction
            gene_direction[gg] = set([direction])
            

print('\ncliff_delta_results: ', len(cliff_delta_results))
cliff_delta_UP = set([x for x in cliff_delta_results if\
                      cliff_delta_results[x]=='UP'])
cliff_delta_DOWN = set([x for x in cliff_delta_results if\
                        cliff_delta_results[x]=='DOWN'])
    
print('cliff_delta_UP: ', len(cliff_delta_UP))
print('cliff_delta_DOWN: ', len(cliff_delta_DOWN))

################ deseq2 ################
deseq_lfc_results = {}
with open(deseq_output_fn,'r') as infile:
    for i,line in enumerate(infile):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x.replace('"',''):xi+1 for xi,x in enumerate(sp)}
            colidx['Gene|GeneId'] = 0 
        else:
            gg = sp[colidx['Gene|GeneId']].replace('"','')

            symbol, entrez = gg.split('|')
            entrez_id_to_full[entrez] = gg
            if symbol == '?':
                symbol = ''
                gg = symbol+'|'+entrez

            posteriorLFC = float(sp[colidx['posteriorLFC']])   #from lfcShrink
            direction = sign_num(posteriorLFC)
            
            if gg not in gene_effect_size:
                gene_effect_size[gg] = {}
            gene_effect_size[gg]['posteriorLFC'] = posteriorLFC
            
            if abs(posteriorLFC)<post_lfc_threshold:
                if gg not in non_passing_genes:
                    non_passing_genes[gg] = {'posteriorLFC':posteriorLFC}
                else:
                    non_passing_genes[gg]['posteriorLFC'] = posteriorLFC
                continue            
            
            deseq_lfc_results[gg] = direction
            if gg not in gene_direction:
                gene_direction[gg] = set()
            gene_direction[gg].add(direction)
            
print('\ndeseq_lfc_results: ', len(deseq_lfc_results))
DeSeq_UP = set([x for x in deseq_lfc_results if deseq_lfc_results[x]=='UP'])
DeSeq_DOWN = set([x for x in deseq_lfc_results if deseq_lfc_results[x]=='DOWN'])
print('DeSeq_UP: ', len(DeSeq_UP))
print('DeSeq_DOWN: ', len(DeSeq_DOWN))

UP_intersection = DeSeq_UP.intersection(cliff_delta_UP)
DOWN_intersection = DeSeq_DOWN.intersection(cliff_delta_DOWN)
print('\nUP_intersection: ', len(UP_intersection))
print('DOWN_intersection: ', len(DOWN_intersection))

print('\nnon_passing_genes: ', len(non_passing_genes))

conflicting_direction = set()
re_decide_direction = {}
for gene in gene_direction:
    if len(gene_direction[gene])!=1:
        conflicting_direction.add(gene)
 
        cliffd = cliff_delta_magnitude(gene_effect_size[gene]['cliff_d'])
        lfc = deseq_lfc_magnitude(gene_effect_size[gene]['posteriorLFC'])
        
        #if magnitude tied, can't decide, remove
        if cliffd == lfc:
            continue
        else:
            if cliffd>lfc:
                re_decide_direction[gene] = cliff_delta_results[gene]
            elif lfc>cliffd:
                re_decide_direction[gene] = deseq_lfc_results[gene]
                
genes_passing_both_metrics = set(cliff_delta_results.keys()).intersection(set(deseq_lfc_results.keys()))
print('\ngenes_passing_both_metrics: ', len(genes_passing_both_metrics))

print('gene_direction: ', len(gene_direction))        

print('\nconflicting_direction genes: ', len(conflicting_direction)) 
print('\ngenes w/ conflicting direction, redecided based on magnitude: ', len(re_decide_direction))

UP_combined = set()
UP_combined.update(cliff_delta_UP)
UP_combined.update(DeSeq_UP)
UP_combined = UP_combined-conflicting_direction

DOWN_combined = set()
DOWN_combined.update(cliff_delta_DOWN)
DOWN_combined.update(DeSeq_DOWN)
DOWN_combined = DOWN_combined-conflicting_direction

if re_decide_direction:
    for gg in re_decide_direction:
        if re_decide_direction[gg]=='UP':
            UP_combined.add(gg)
        elif re_decide_direction[gg]=='DOWN':
            DOWN_combined.add(gg)

print('\n##########')
print('UP_combined: ', len(UP_combined))   
print('DOWN_combined: ', len(DOWN_combined))

conflicting_direction = conflicting_direction-set(re_decide_direction.keys())
print('conflicting_direction genes: ', len(conflicting_direction)) 

mwu_direction_list = ['NotSignificant|DOWN',
'downregulated(ecDNA)',
'upregulated(ecDNA)',
'NotSignificant|UP',
'NotSignificant|NA']

map_direction = {'downregulated(ecDNA)':'DOWN',
                 'upregulated(ecDNA)':'UP'}
      
################################################################
All_genes_with_direction = set()
All_genes_with_direction.update(UP_combined)
All_genes_with_direction.update(DOWN_combined)
print('\n\nAll_genes_with_direction: ', len(All_genes_with_direction))
print('all genes: ', len(gene_effect_size))

with open(os.path.join(gene_directionality_dir, 'cliffd_lfc_values.txt'),'w') as outfile:
    outfile.write('\t'.join(['Direction', 'Gene', 'Symbol', 'EntrezID',\
                             'Cliff_delta','Posterior_LFC'])+'\n')
    for gg in All_genes_with_direction:
        
        direction = ''
        if gg in UP_combined:
            direction = 'UP'
        if gg in DOWN_combined:
            direction = 'DOWN'
            
        Cliff_delta = gene_effect_size[gg]['cliff_d']
        
        Posterior_LFC = ''
        if 'posteriorLFC' in gene_effect_size[gg]:
            Posterior_LFC = gene_effect_size[gg]['posteriorLFC']
        
        symbol, entrez = gg.split('|') 
        outfile.write('\t'.join([direction, gg, symbol, entrez,\
                                 repr(Cliff_delta), repr(Posterior_LFC)])+'\n')
    
    for gg in conflicting_direction:
        if gg in All_genes_with_direction:
            raise ValueError(gg)
            
        direction = 'conflicting_direction'
        symbol, entrez = gg.split('|') 
        Cliff_delta = gene_effect_size[gg]['cliff_d']
        Posterior_LFC = gene_effect_size[gg]['posteriorLFC']
        
        outfile.write('\t'.join([direction, gg, symbol, entrez,\
                                 repr(Cliff_delta), repr(Posterior_LFC)])+'\n')

    for gg in non_passing_genes:
        if gg in All_genes_with_direction:
            continue
        
        direction = 'negligible_effect_size'
        symbol, entrez = gg.split('|') 
        Cliff_delta = gene_effect_size[gg]['cliff_d']
        Posterior_LFC = gene_effect_size[gg]['posteriorLFC']
        
        outfile.write('\t'.join([direction, gg, symbol, entrez,\
                                 repr(Cliff_delta), repr(Posterior_LFC)])+'\n')