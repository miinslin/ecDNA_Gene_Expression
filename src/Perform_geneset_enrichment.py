# -*- coding: utf-8 -*-
import os
import argparse
import pandas as pd
from geneset_enrichment_function import geneset_enrichment, parse_GO_db_terms
import itertools 
from sklearn.metrics import cohen_kappa_score
from boruta_output_parser import parse_boruta_featurestxt

base_dir = os.path.abspath(os.path.join(__file__, "../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--matrix", dest = "matrix",\
                    default = os.path.join(base_dir, 'data', \
                    'RNAseq_expression_matrix',\
                    'AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt'),\
                    help="Gene expression matrix path.")
args = parser.parse_args()

num_runs = 200
percentage_of_runs = 0.05
min_runs = int(num_runs*percentage_of_runs)

TCGA_matrix_fn = os.path.normpath(args.matrix)
matrix_type = os.path.basename(TCGA_matrix_fn).replace('.txt','')

Gene_set_enrichment_dir = os.path.join(base_dir, 'data', matrix_type, 'Gene_Set_Enrichment')
if not os.path.isdir(Gene_set_enrichment_dir):
    os.makedirs(Gene_set_enrichment_dir)

#GOBP mappings
GOdb_terms = parse_GO_db_terms()
goterm_to_id = GOdb_terms['goterm_to_id']

#Core genes (i.e., gene identified in at least min_runs boruta runs)
Boruta_runs_dir = os.path.join(base_dir, 'data', matrix_type, 'Boruta_Trials')
Bortua_genes_dir = os.path.join(base_dir, 'data', matrix_type, 'Boruta_Genes')
boruta_results = parse_boruta_featurestxt(Boruta_runs_dir, Bortua_genes_dir, num_runs)
Core_genes = sorted([x.split('|')[1] for x in boruta_results if len(boruta_results[x])>=min_runs])
print('Core_genes: ', len(Core_genes)) 

#CorEx genes
au_threshold = 95
coexpression_dir = os.path.join(base_dir, 'data', matrix_type, 'co_expression')
corex_gene_summary = os.path.join(coexpression_dir, 'au'+str(au_threshold)+'_cluster_boruta_genes.txt')
corex_gene_summary_df = pd.read_csv(corex_gene_summary, sep='\t', header=0, dtype='str')
corex_gene_summary_df = corex_gene_summary_df.astype({'#Trials(200)_by_cluster':float, '#Trials(200)':float})
CorEx_genes_df = corex_gene_summary_df.loc[corex_gene_summary_df['#Trials(200)_by_cluster']>=min_runs]
CorEx_genes = CorEx_genes_df['Gene|GeneId'].tolist()
CorEx_genes_cluster_run_counts = CorEx_genes_df['#Trials(200)_by_cluster'].tolist()
CorEx_genes_boruta_run_counts = CorEx_genes_df['#Trials(200)'].tolist()
CorEx_genes = {g.split('|')[1]:{'cluster_trial_counts':CorEx_genes_cluster_run_counts[i],'gene_trial_counts':CorEx_genes_boruta_run_counts[i]} for i,g in enumerate(CorEx_genes)}
print('CorEx_genes: ', len(CorEx_genes))

#gene directionality
gene_directionality_fn = os.path.join(base_dir, 'data', matrix_type,\
                                       'gene_directionality', 'cliffd_lfc_values.txt')
UP_combined = set()
DOWN_combined = set()
entrez_id_to_full = {}
gene_effect_size = {}
with open(gene_directionality_fn,'r') as infile:
    for i,line in enumerate(infile):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            Direction = sp[colidx['Direction']]    
            Gene = sp[colidx['Gene']]
            EntrezID = sp[colidx['EntrezID']]
            if EntrezID not in CorEx_genes:
                continue
            entrez_id_to_full[EntrezID] = Gene   
            if Direction == 'UP':
                UP_combined.add(EntrezID)
            if Direction == 'DOWN':
                DOWN_combined.add(EntrezID)
            gene_effect_size[EntrezID] = {'cliff_d':sp[colidx['Cliff_delta']],
                                          'posteriorLFC':sp[colidx['Posterior_LFC']]}
if not gene_effect_size:
    raise ValueError('Error: No genes included...')
print('UP_combined: ', len(UP_combined))
print('DOWN_combined: ', len(DOWN_combined))

All_genes = set()
All_genes.update(UP_combined)
All_genes.update(DOWN_combined)
All_genes = sorted([(entrez_id_to_full[x], x) for x in All_genes])
All_genes = [x[1] for x in All_genes]

up_regulated_gene_entrez_fn = os.path.join(Gene_set_enrichment_dir, 'UP_genes_entrez.txt')
with open(up_regulated_gene_entrez_fn,'w') as outfile:
    for entrez in UP_combined:
        outfile.write(entrez+'\n')

down_regulated_gene_entrez_fn = os.path.join(Gene_set_enrichment_dir, 'DOWN_genes_entrez.txt')
with open(down_regulated_gene_entrez_fn,'w') as outfile:
    for entrez in DOWN_combined:
        outfile.write(entrez+'\n')

up_down_gene_fn = os.path.join(Gene_set_enrichment_dir, 'UP_DOWN_genes.txt')
with open(up_down_gene_fn,'w') as outfile:
    outfile.write('\t'.join(['Gene', 'Symbol', 'EntrezID', 'Direction', 'Trial_count_by_gene(200)', 'Trial_count_by_cluster(200)','Cliff_delta','Posterior_LFC'])+'\n')
    for entrez in All_genes:
        direction = ''
        if entrez in UP_combined:
            direction = 'UP'
        if entrez in DOWN_combined:
            direction = 'DOWN'
            
        Cliff_delta = 'NA'
        if 'cliff_d' in gene_effect_size[entrez]:
            Cliff_delta = gene_effect_size[entrez]['cliff_d']
        Posterior_LFC = 'NA'
        if 'posteriorLFC' in gene_effect_size[entrez]:
            Posterior_LFC = gene_effect_size[entrez]['posteriorLFC']
        gg = entrez_id_to_full[entrez]
        symbol,entrez = gg.split('|')
        
        run_count_by_gene = ''
        run_count_by_cluster = ''
        if entrez in CorEx_genes:
            run_count_by_gene = CorEx_genes[entrez]['gene_trial_counts']
            run_count_by_cluster = CorEx_genes[entrez]['cluster_trial_counts']
        
        outfile.write('\t'.join([gg, symbol, entrez, direction, \
                                 str(run_count_by_gene), str(run_count_by_cluster),\
                                 repr(Cliff_delta), repr(Posterior_LFC)])+'\n')
        
        
#N is the total number of gene universe (all measurable genes in rnaseq matrix)
TCGA_matrix_df = pd.read_csv(TCGA_matrix_fn, sep='\t', header=0, dtype='str')
measurable_genes = TCGA_matrix_df['Gene_name|Gene_id'].tolist()
#print('\nmeasurable_genes: ', len(measurable_genes)) #16309
N = len(measurable_genes)
 
#geneset enrichment for up-regulated genes
UP_geneset_enrichment = geneset_enrichment(N, UP_combined, entrez_id_to_full, Core_genes, Gene_set_enrichment_dir, direction = 'UP')
#geneset enrichment for down-regulated genes
DOWN_geneset_enrichment = geneset_enrichment(N, DOWN_combined, entrez_id_to_full, Core_genes, Gene_set_enrichment_dir, direction = 'DOWN')


####################################################################
initial_group_threshold = 1
seeding_group_condition_2 = 0.5
 
geneset_enrichment_results = {'UP':UP_geneset_enrichment['geneset_enrichment_results'],
                              'DOWN':DOWN_geneset_enrichment['geneset_enrichment_results']}
geneset_gene_lists = {'UP':UP_geneset_enrichment['geneset_gene_lists'],\
                      'DOWN':DOWN_geneset_enrichment['geneset_gene_lists']}

def group_genesets(direction, kappa_score_level_threshold, multiple_linkage_threshold):
    print('##'+direction)

    sub_gene_set_enrichment_dir = os.path.join(Gene_set_enrichment_dir, 'groups',\
                                               'kappa_score_'+str(kappa_score_level_threshold)+'_multiple_linkage_'+str(multiple_linkage_threshold))
    if not os.path.isdir(sub_gene_set_enrichment_dir):
        os.makedirs(sub_gene_set_enrichment_dir)
      
    grouped_geneset_fn = os.path.join(sub_gene_set_enrichment_dir, direction+'_geneset_groups_category.txt')   
    grouped_geneset_list_fn = os.path.join(sub_gene_set_enrichment_dir, direction+'_geneset_groups_list.txt')  
       
    DEG_genesets = sorted(geneset_enrichment_results[direction].keys()) #99
    DEG_overlap_genes_all = sorted(geneset_gene_lists[direction])
    DEG_matrix = {}
    for geneset in DEG_genesets:
        geneset_genes = geneset_enrichment_results[direction][geneset]['geneset_genes']
        geneset_v = [1 if entrez in geneset_genes else 0 for entrez in DEG_overlap_genes_all]
        DEG_matrix[geneset[5:]] = geneset_v

    DEG_genesets = [x[5:] for x in DEG_genesets]
    
    #write matrix to file
    matrix_fn = os.path.join(sub_gene_set_enrichment_dir, direction+'_geneset_genes_matrix.txt')
    with open(matrix_fn,'w') as outfile:
        outfile.write('geneset'+'\t'+'\t'.join([entrez_id_to_full[x] for x in DEG_overlap_genes_all])+'\n')
        for geneset in DEG_matrix:
            outfile.write(geneset+'\t'+'\t'.join([str(x) for x in DEG_matrix[geneset]])+'\n')
    
    #genes in enriched genesets
    print('DEG_genesets: ', len(DEG_genesets))
    # print('DEG_overlap_genes_all: ', len(DEG_overlap_genes_all))

    # STEP 1: compute term-term kappa scores
    kappa_scores = {}
    for i in range(len(DEG_genesets)):
        g1_name = DEG_genesets[i]
        kappa_scores[g1_name] = {}
        for j in range(len(DEG_genesets)):
            g2_name = DEG_genesets[j]
            g1 = DEG_matrix[g1_name]
            g2 = DEG_matrix[g2_name]
            kappa_scores[g1_name][g2_name] = cohen_kappa_score(g1, g2)

    # STEP 2: create initial seeding groups
    initial_seeds = []
    outlier_genesets = []

    for g1 in kappa_scores:#99
        close_genesets = [g2 for g2 in kappa_scores[g1] if kappa_scores[g1][g2]>=kappa_score_level_threshold]
        
        close_geneset_DEG_overlap_genes = set()
        for g in close_genesets:
            close_geneset_DEG_overlap_genes.update(geneset_enrichment_results[direction]['GOBP_'+g]['overlap_genes'])
            
        if not set(close_genesets)-set([g1]):
            outlier_genesets.append(g1)
            continue
        else:
            condition_2_pairs = list(itertools.combinations(close_genesets, 2))
            condition_2_values = []
            c=0
            for x,y in condition_2_pairs:
                if x==y:
                    continue
                value = kappa_scores[x][y]
                c+=1
                if value>=kappa_score_level_threshold:
                    condition_2_values.append((x,y,value))
 
            condition_2_percentage = float(len(condition_2_values))/float(c)
            if condition_2_percentage>=seeding_group_condition_2:
                initial_seeds.append((len(close_geneset_DEG_overlap_genes), set(close_genesets)))
    
    # print('initial_seeds: ', len(initial_seeds))
 
    #Step3:  Iteratively merging above qualified seeds:
    grouped_seeds = sorted(initial_seeds, reverse=True)
    grouped_seeds = [x[1] for x in grouped_seeds]
    
    while grouped_seeds:
        grouped_i = set()
        grouped_members = set() 
        ungrouped_i = range(len(grouped_seeds))
        ungrouped_i_pairs = list(itertools.combinations(ungrouped_i, 2))
        for x,y in ungrouped_i_pairs:
            seed_1 = grouped_seeds[x]
            seed_2 = grouped_seeds[y]
            #>50% of smaller group
            smaller_seed_group = min([len(seed_1), len(seed_2)])
            seed_intersection = seed_1.intersection(seed_2)
            overlap_p = float(len(seed_intersection)/smaller_seed_group)
            if overlap_p>=multiple_linkage_threshold:
                grouped_i.add(x)
                grouped_i.add(y)
                grouped_members.update(seed_1)
                grouped_members.update(seed_2)
                break
        #update
        if grouped_members:
            grouped_seeds = [grouped_seeds[i] for i in range(len(grouped_seeds)) if i not in grouped_i]
            grouped_seeds.append(grouped_members)
        else:
            break

    grouped_seeds = sorted(grouped_seeds)
    
    # #check genesets overlapped in multiple groups
    # print('check genesets overlapped in multiple groups:')
    # multi_group_count = {}
    # for xi,group in enumerate(grouped_seeds):
    #     for s in group:
    #         if s not in multi_group_count:
    #             multi_group_count[s] = set()
    #         multi_group_count[s].add(xi)
    # for s in multi_group_count:
    #     if len(multi_group_count[s])>1:
    #         print(s, multi_group_count[s])
       
    #representative based on pvalue
    with open(grouped_geneset_fn,'w') as outfile, open(grouped_geneset_list_fn, 'w') as outfile_list:
        outfile.write('\t'.join(['index', 'rep_geneset', 'rep_padj', 'group_n', 'group_members', 'group_members_goid'])+'\n')
        for i,group in enumerate(grouped_seeds):
            group_n = len(group)
            
            group_scores = sorted([(geneset_enrichment_results[direction]['GOBP_'+member]['padj'], member) for member in group])
            group_rep = group_scores[0][1]
            group_rep_padj = geneset_enrichment_results[direction]['GOBP_'+group_rep]['padj']
            group = [x[1] for x in group_scores]
            goids = [goterm_to_id[s] for s in group]
            outfile.write('\t'.join([str(i), group_rep, repr(group_rep_padj), str(group_n), ';'.join(group), ';'.join(goids)])+'\n')
            outfile_list.write('#rep:'+group_rep+'\n'+'\n'.join(['\t'+x[1]+'\t'+repr(x[0])+'\t'+goterm_to_id[x[1]]+'\t'+\
                                                                     ';'.join(sorted([entrez_id_to_full[gg] for gg in geneset_enrichment_results[direction]['GOBP_'+x[1]]['overlap_genes']]))\
                                                                         for x in group_scores])+'\n'+'#END\n')
        #write outliers
        for j,geneset in enumerate(outlier_genesets):
            rep_padj = geneset_enrichment_results[direction]['GOBP_'+geneset]['padj']
            goid = goterm_to_id[geneset]
            gene_list = ';'.join(sorted([entrez_id_to_full[gg] for gg in geneset_enrichment_results[direction]['GOBP_'+geneset]['overlap_genes']]))
            outfile.write('\t'.join(['',geneset, repr(rep_padj), str(1), '', goid])+'\n')
            outfile_list.write('\t'.join([geneset, repr(rep_padj), goid, gene_list])+'\n')

    return None

up_results = group_genesets('UP', 0.6, 0.5)
down_results = group_genesets('DOWN', 0.5, 0.25)