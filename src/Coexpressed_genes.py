# -*- coding: utf-8 -*-
import os
import argparse
from boruta_output_parser import parse_boruta_featurestxt

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--boruta_dir", dest = "boruta",\
                    default = os.path.join(base_dir, 'data', \
                    'AC049_TCGA_KNN_imputed_min_ecdna_3_matrix',\
                    'Boruta_Runs'),\
                    help="Path to Boruta trial output files.")
parser.add_argument("-n", "--num_boruta_trials", dest = "num_runs",\
                    default = 200,\
                    help="Number of Boruta Trials")
args = parser.parse_args()

min_runs = 10
num_runs = int(args.num_runs)
print('Number of Boruta Trials: ', num_runs)
au_threshold = 95
gene_column = '#Gene_name|Gene_id'

Boruta_runs_dir = os.path.normpath(args.boruta)
matrix_type = os.path.basename(os.path.abspath(os.path.join(Boruta_runs_dir ,"..")))

coexpression_dir = os.path.join(base_dir, 'data', matrix_type, 'co_expression')
pvclust_au = os.path.join(coexpression_dir, 'au_'+str(au_threshold)+'.txt')

pvclust_au_boruta_summary = os.path.join(coexpression_dir, 'au'+str(au_threshold)+'_cluster_boruta_genes.txt')
pvclust_au_boruta_indicator = os.path.join(coexpression_dir, 'au_'+str(au_threshold)+'_clusters.txt')

Bortua_genes_dir = os.path.join(base_dir, 'data', matrix_type, 'Boruta_Genes')
if not os.path.isdir(Bortua_genes_dir):
    os.makedirs(Bortua_genes_dir)

#pcvlust clusters where au>=95
pvclust_edges = {}
with open(pvclust_au,'r') as infile:
    for i,line in enumerate(infile):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            gene = sp[colidx['gene']]
            edge = sp[colidx['edge']]
            edge_au = sp[colidx['edge_au']]
            edge_si = sp[colidx['edge_si']]
            if edge not in pvclust_edges:
                pvclust_edges[edge] = {'au':edge_au,'si':edge_si,'members':set()}
            pvclust_edges[edge]['members'].add(gene)
print('pvclust_edges: ', len(pvclust_edges)) #1374 edges

#Boruta results
boruta_results = parse_boruta_featurestxt(Boruta_runs_dir, Bortua_genes_dir, num_runs)

Boruta_genes = set()
Core_genes = set()
for gg in boruta_results:
    if len(boruta_results[gg])>=min_runs:
        Core_genes.add(gg)
    Boruta_genes.add(gg)
print('Core_genes: ', len(Core_genes)) #408
print('Boruta_genes: ', len(Boruta_genes))

#compute runs by cluster
edge_boruta_runs = {}
boruta_genes_in_pvclust_clusters = set()
for edge in pvclust_edges:
    cluster_members = pvclust_edges[edge]['members']
    cluster_boruta_genes = cluster_members.intersection(Boruta_genes)
    boruta_genes_in_pvclust_clusters.update(cluster_boruta_genes)
    for gene in cluster_boruta_genes:
        if edge not in edge_boruta_runs:
            edge_boruta_runs[edge] = set()
        edge_boruta_runs[edge].update(boruta_results[gene])

BorutaGeneNotInCluster = Boruta_genes-boruta_genes_in_pvclust_clusters
for gene in BorutaGeneNotInCluster:
    edge_boruta_runs[gene] = boruta_results[gene]

edge_runcount = {}
gene_to_edge = {}
for edge in edge_boruta_runs:
    num_runs_by_cluster = len(edge_boruta_runs[edge])
    if num_runs_by_cluster < min_runs:
        continue
    edge_runcount[edge] = num_runs_by_cluster
    if edge in BorutaGeneNotInCluster:
        pvclust_edges[edge] = {'au': 'NA', 'si': 'NA', 'members': set([edge])}
    for m in pvclust_edges[edge]['members']:
        gene_to_edge[m] = edge

edges_passing_threshold = set(edge_runcount.keys())
print('edge_passing_threshold: ', len(edges_passing_threshold))
# edge_passing_threshold:  354

remove_edges = set(pvclust_edges.keys()) - edges_passing_threshold
for edge in remove_edges:
    del pvclust_edges[edge]

if set(pvclust_edges.keys())!=edges_passing_threshold:
    raise ValueError('edge inconsistency...')

CorEx_genes = set(gene_to_edge.keys())
print('CorEx_genes: ', len(CorEx_genes))

edge_to_clusterid = sorted([(len(pvclust_edges[edge]['members']),
                          len(pvclust_edges[edge]['members'].intersection(Core_genes)), edge) for edge in pvclust_edges], reverse=True)
edge_to_clusterid = {entry[2]:str(i+1) for i,entry in enumerate(edge_to_clusterid)}

with open(pvclust_au_boruta_indicator,'w') as outfile:
    outfile.write('\t'.join(['Cluster#', 'au', 'si', 'trial_count_by_cluster',
                             '#Boruta_genes', '#Core_genes', '#CorEx_genes',
                             'Boruta_genes', 'Core_genes', 'CorEx_genes'])+'\n')
    for edge in edge_to_clusterid:
        members = pvclust_edges[edge]['members'] #includes boruta genes
        boruta_in_cluster = members.intersection(Boruta_genes)
        core_in_cluster = members.intersection(Core_genes)
        runcount = edge_runcount[edge]
        au_value = pvclust_edges[edge]['au']
        si_value = pvclust_edges[edge]['si']
        outfile.write('\t'.join([edge_to_clusterid[edge], au_value, si_value, str(runcount),\
                                 str(len(boruta_in_cluster)), str(len(core_in_cluster)), str(len(members)),
                                 ';'.join(sorted(boruta_in_cluster)),  ';'.join(sorted(core_in_cluster)),
                                 ';'.join(sorted(members))])+'\n')

with open(pvclust_au_boruta_summary, 'w') as outfile:
    outfile.write('#' + '\t'.join(['Cluster#', 'Gene|GeneId', '#Trials(200)', '#Trials(200)_by_cluster']) + '\n')
    for gene in gene_to_edge:
        edge = gene_to_edge[gene]
        num_runs_by_gene = '0'
        if gene in Boruta_genes:
            num_runs_by_gene = len(boruta_results[gene])
        num_runs_by_cluster = edge_runcount[edge]
        outfile.write('\t'.join([edge_to_clusterid[edge], gene, str(num_runs_by_gene), str(num_runs_by_cluster)]) + '\n')