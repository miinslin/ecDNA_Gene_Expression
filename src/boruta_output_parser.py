import pandas as pd
import os

def parse_runs_summary(boruta_summary_fn):
    df = pd.read_csv(boruta_summary_fn, sep='\t', header=0, dtype='str')
    boruta_gene_runids = {}
    for i,entry in df.iterrows():
        gene = entry['#Gene|GeneId']
        runs = entry['Trials'].split(';')
        boruta_gene_runids[gene] = runs
    return boruta_gene_runids

def parse_boruta_featurestxt(Boruta_runs_dir, Bortua_genes_dir, num_trials):
    gene_column = '#Gene_name|Gene_id'

    boruta_summary_fn = os.path.join(Bortua_genes_dir, 'Boruta_' + str(num_trials) + '_trials_summary.txt')

    ##########################################
    boruta_runs_genes = {}
    for b_run in range(0, num_trials, 1):
        data_run = 'data_' + str(b_run)

        BorutaPy_features_fn = os.path.join(Boruta_runs_dir, data_run, data_run + '_BorutaPy_features.txt')

        with open(BorutaPy_features_fn, 'r') as infile:
            for i, line in enumerate(infile):
                sp = line.rstrip('\n').split('\t')
                if i == 0:
                    colidx = {x: xi for xi, x in enumerate(sp)}
                else:
                    gg = sp[colidx[gene_column]]
                    two_step_corr = sp[colidx['two_step_corr']]
                    if two_step_corr != 'True':
                        continue
                    if gg not in boruta_runs_genes:
                        boruta_runs_genes[gg] = set()
                    boruta_runs_genes[gg].add(data_run)

    with open(boruta_summary_fn, 'w') as outfile:
        outfile.write('\t'.join(['#Gene|GeneId', '#Trials(' + str(num_trials) + ')', 'Trials']) + '\n')
        for gene in sorted(boruta_runs_genes):
            data_runs = boruta_runs_genes[gene]
            outfile.write('\t'.join([gene, str(len(data_runs)), ';'.join(sorted(data_runs))]) + '\n')

    return parse_runs_summary(boruta_summary_fn)