# -*- coding: utf-8 -*-
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
import numpy as np
import os

base_dir = os.path.abspath(os.path.join(__file__, "../.."))

def parse_GO_db_terms():
    GO_db_terms = os.path.join(base_dir, 'data', 'reference', 'gene_ontology', 'go-basic.id_to_name.txt')
    
    """## **Go term mappings**"""
    # GO_id	name	namespace	alt_ids
    # GO:0000001	mitochondrion inheritance	biological_process	
    
    goterm_to_id = {}
    goid_to_name = {}
    with open(GO_db_terms,'r') as infile:
        for i,line in enumerate(infile):
            sp = line.rstrip('\n').split('\t')
            if i==0:
                colidx = {x:xi for xi,x in enumerate(sp)}
            else: 
                goterm_id = sp[colidx['GO_id']]
                original_name = sp[colidx['name']]
                alt_ids = sp[colidx['alt_ids']].split(';')
                all_ids = [goterm_id]
                all_ids.extend(alt_ids)
                namespace = sp[colidx['namespace']]
                if namespace == 'biological_process':
                    namespace = 'GOBP'
                else:
                    continue
                
                goterm_name = ''.join(['_' if x==' ' else x.upper() for x in sp[colidx['name']].replace(',','').replace("'",' ').replace('-',' ').replace('/',' ')]).replace('__','_')

                if 'OBSOLETE_' in goterm_name:
                    goterm_name = goterm_name.replace('OBSOLETE_', '')
                if 'obsolete ' in original_name:
                    original_name = original_name.replace('obsolete ','')
                
                goterm_to_id[goterm_name] = goterm_id #use main id
                for g_id in all_ids:
                    goid_to_name[g_id] = original_name

    
    #print('goid_to_name: ', len(goid_to_name)) #32,690
    
    goterm_to_id['CENP_A_CONTAINING_CHROMATIN_ORGANIZATION'] = 'GO:0061641'
    goterm_to_id['CHROMATIN_ASSEMBLY_OR_DISASSEMBLY'] = 'GO:0006333'
    goterm_to_id['DNA_DEPENDENT_DNA_REPLICATION'] = 'GO:0006261'
    goterm_to_id['DNA_DEPENDENT_DNA_REPLICATION_MAINTENANCE_OF_FIDELITY'] = 'GO:0045005'
    goterm_to_id['DNA_PACKAGING'] = 'GO:0006323'
    goterm_to_id['DNA_REPLICATION_INDEPENDENT_CHROMATIN_ORGANIZATION'] = 'GO:0034724'
    goterm_to_id['NEGATIVE_REGULATION_OF_DNA_DEPENDENT_DNA_REPLICATION'] = 'GO:2000104'
    goterm_to_id['REGULATION_OF_DNA_DEPENDENT_DNA_REPLICATION'] = 'GO:0090329'
    goterm_to_id['NUCLEAR_TRANSCRIBED_MRNA_POLY_A_TAIL_SHORTENING'] = 'GO:0000289'
    goterm_to_id['MATURATION_OF_LSU_RRNA_FROM_TRICISTRONIC_RRNA_TRANSCRIPT_SSU_RRNA_5_8S_RRNA_LSU_RRNA'] = 'GO:0000463'
    goterm_to_id['RNA_DEPENDENT_DNA_BIOSYNTHETIC_PROCESS'] = 'GO:0006278'
    goterm_to_id['POSITIVE_REGULATION_OF_CELLULAR_PROTEIN_CATABOLIC_PROCESS'] = 'GO:1903364'
    goterm_to_id['POSITIVE_REGULATION_OF_PROTEOLYSIS_INVOLVED_IN_CELLULAR_PROTEIN_CATABOLIC_PROCESS'] = 'GO:1903052'
    goterm_to_id['POSITIVE_REGULATION_OF_DNA_DEPENDENT_DNA_REPLICATION'] = 'GO:2000105'
    goterm_to_id['PROTEIN_MATURATION_BY_4FE_4S_CLUSTER_TRANSFER'] = 'GO:0106035'
    goterm_to_id['MATURATION_OF_5_8S_RRNA'] = 'GO:0000460'
    goterm_to_id['POSITIVE_REGULATION_OF_CHEMOKINE_C_X_C_MOTIF_LIGAND_2_PRODUCTION'] = 'GO:2000343'
    
    return {'goid_to_name':goid_to_name,'goterm_to_id':goterm_to_id}

def geneset_enrichment(N, query, entrez_id_to_full, Core_genes, Gene_set_enrichment_dir, direction = ''):
    
    msigdb_fn = os.path.join(base_dir, 'data', 'reference', 'GSEA_msigdb', 'c5.go.bp.v7.5.1.entrez.gmt')
    
    GOdb_terms = parse_GO_db_terms()
    goid_to_name = GOdb_terms['goid_to_name']
    goterm_to_id = GOdb_terms['goterm_to_id']

    term_entrez_ids = {}
    with open(msigdb_fn, 'r') as infile:
        for line in infile:
            sp = line.strip().split('\t')
            NAME = sp[0]
            ids = sp[2:]
            term_entrez_ids[NAME] = ids    
    
    genes_n = len(query)
    
    #print('N=', N)
    #print('n=', genes_n)
    
    #print('gene enrichment. hypergeometric. one-sided Fishers exact test...')
    gset_pvalues = []
    for geneset in term_entrez_ids:
        geneset_genes = set(term_entrez_ids[geneset])
        #m is the number of genes in the set from MSigDB
        m = len(geneset_genes)
        #k is the number of genes in the intersection of the query set with a set from MSigDB
        k_ = len(geneset_genes.intersection(query))
        if k_!=0:            
            table = np.array([[k_, m-k_],
                              [genes_n-k_, N+k_-genes_n-m]])            
            oddsr, p = fisher_exact(table, alternative='greater')            
            gset_pvalues.append((p, geneset, geneset_genes.intersection(query)))

    if not gset_pvalues:
        return 'none'
 
    #print('pvalues adjusted for multiple hypothesis testing to limit FDR...')
    fdr_alpha = 0.05
    #print('fdr_alpha=',fdr_alpha)
    padj_p = 0.05
    #print('padj_p=',padj_p)

    fdr_correction_ = fdrcorrection([x[0] for x in gset_pvalues], alpha=fdr_alpha)
    corrected_p_values = fdr_correction_[1]
    
    geneset_enrichment_results = {}
    non_enriched_genesets = {}
    enriched_geneset_order = []
    non_enriched_geneset_order = []
    geneset_gene_lists = set()
    for i,p in enumerate(gset_pvalues):
        geneset = gset_pvalues[i][1]
        raw_pvalue = gset_pvalues[i][0]
        overlap_genes = gset_pvalues[i][2]
        Geneset_genes = term_entrez_ids[geneset]
        padj = corrected_p_values[i]

        if padj<padj_p:
            geneset_enrichment_results[geneset] = {'geneset_genes':Geneset_genes, 'overlap_genes':overlap_genes, 'padj':padj, 'raw_pvalue':raw_pvalue}
            enriched_geneset_order.append((padj,geneset))
            geneset_gene_lists.update(overlap_genes)        
        else:
            non_enriched_genesets[geneset] = {'geneset_genes':Geneset_genes, 'overlap_genes':overlap_genes, 'padj':padj, 'raw_pvalue':raw_pvalue}
            non_enriched_geneset_order.append((padj,geneset))

    enriched_geneset_order = [x[1] for x in sorted(enriched_geneset_order)]
    genesets_enriched_fn = os.path.join(Gene_set_enrichment_dir, direction+'_genes_enriched_GOBP_terms.txt')    
    with open(genesets_enriched_fn,'w') as outfile:
        outfile.write('\t'.join(['id','term','p_value','padj', 'gene_set_size', 'DEG_overlap_k', 'DEG_overlap_k_CORE', 'DEG_overlap_genes'])+'\n')
        for geneset in enriched_geneset_order:
            entry = geneset_enrichment_results[geneset]
            DEG_overlap_genes_CORE = []
            if not Core_genes:
                DEG_overlap_genes_CORE = []
            else:
                DEG_overlap_genes_CORE = [entrez_id_to_full[entrez] for entrez in entry['overlap_genes'] if entrez in Core_genes]
            DEG_overlap_genes = sorted([entrez_id_to_full[entrez] for entrez in entry['overlap_genes']])
            gobp_id = goterm_to_id[geneset.replace('GOBP_','')]
            lower_case_term = goid_to_name[gobp_id]
            entry = '\t'.join([gobp_id, lower_case_term,\
                                repr(entry['raw_pvalue']), repr(entry['padj']),\
                                str(len(entry['geneset_genes'])),\
                                str(len(DEG_overlap_genes)),\
                                str(len(DEG_overlap_genes_CORE)),\
                                ';'.join(DEG_overlap_genes)])
            outfile.write(entry+'\n')  

    non_enriched_geneset_order = [x[1] for x in sorted(non_enriched_geneset_order)]
    genesets_NOT_enriched_fn = os.path.join(Gene_set_enrichment_dir, direction+'_genes_NOT_enriched_GOBP_terms.txt')
    with open(genesets_NOT_enriched_fn,'w') as outfile:
        outfile.write('\t'.join(['id','term','p_value','padj', 'gene_set_size', 'DEG_overlap_k', 'DEG_overlap_k_CORE', 'DEG_overlap_genes'])+'\n')
        for geneset in non_enriched_geneset_order:
            entry = non_enriched_genesets[geneset]
            DEG_overlap_genes_CORE = []
            if not Core_genes:
                DEG_overlap_genes_CORE = []
            else:
                DEG_overlap_genes_CORE = [entrez_id_to_full[entrez] for entrez in entry['overlap_genes'] if entrez in Core_genes]
            DEG_overlap_genes = sorted([entrez_id_to_full[entrez] for entrez in entry['overlap_genes']])
            gobp_id = ''
            if geneset.replace('GOBP_','') in goterm_to_id:
                gobp_id = goterm_to_id[geneset.replace('GOBP_','')]
            lower_case_term = goid_to_name[gobp_id]
            entry = '\t'.join([gobp_id, lower_case_term,\
                                repr(entry['raw_pvalue']), repr(entry['padj']),\
                                str(len(entry['geneset_genes'])),\
                                str(len(DEG_overlap_genes)),\
                                str(len(DEG_overlap_genes_CORE)),\
                                ';'.join(DEG_overlap_genes)])
            outfile.write(entry+'\n')

    return {'geneset_enrichment_results':geneset_enrichment_results,
            'geneset_gene_lists':geneset_gene_lists,
            'non_enriched_genesets':non_enriched_genesets}
