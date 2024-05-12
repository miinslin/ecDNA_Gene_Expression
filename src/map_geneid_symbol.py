# -*- coding: utf-8 -*-
import os

base_dir = os.path.abspath(os.path.join(__file__, "../.."))

hgnc_complete_set_fn = os.path.join(base_dir, 'data', 'reference', \
                                    'gene_id_mappings', 'hgnc_complete_set_7.24.2018.txt')


def map_gene_ids(protein_coding_bool=False):
    human_entrez = set()
    locus_group_types = set()
    hgnc_entrez = {}
    hgnc_symbol = {}
    hgnc_symbol_ALL = {}
    hgnc_entrez_ALL = {}
    prev_symbols_ALL = {}
    alias_symbols_ALL = {}
    with open(hgnc_complete_set_fn, 'r', encoding='utf-8') as infile:
        for i, line in enumerate(infile):
            sp = line.rstrip('\n').split('\t')
            if i == 0:
                colidx = {x: xi for xi, x in enumerate(sp)}
                # symbol entrez_id  ensembl_gene_id
            else:
                symbol = sp[colidx['symbol']]
                entrez_id = sp[colidx['entrez_id']]
                human_entrez.add(entrez_id)
                ensembl_gene_id = sp[colidx['ensembl_gene_id']]
                locus_group = sp[colidx['locus_group']]
                locus_group_types.add(locus_group)
                name = sp[colidx['name']]
                prev_symbols = sp[colidx['prev_symbol']].strip().split('|')
                alias_symbols = sp[colidx['alias_symbol']].strip().replace('"', '').split('|')

                hgnc_symbol_ALL[symbol.upper()] = {'entrez_id': entrez_id, 'symbol': symbol, 'name': name,
                                                   'ensembl_gene_id': ensembl_gene_id, 'locus_group': locus_group}
                hgnc_entrez_ALL[entrez_id] = {'symbol': symbol, 'name': name, 'ensembl_gene_id': ensembl_gene_id,
                                              'locus_group': locus_group}

                for x in prev_symbols:
                    if x.strip():
                        prev_symbols_ALL[x.upper()] = {'entrez_id': entrez_id, 'symbol': symbol, 'name': name,
                                                       'ensembl_gene_id': ensembl_gene_id, 'locus_group': locus_group}

                for x in alias_symbols:
                    if x.strip():
                        alias_symbols_ALL[x.upper()] = {'entrez_id': entrez_id, 'symbol': symbol, 'name': name,
                                                        'ensembl_gene_id': ensembl_gene_id, 'locus_group': locus_group}

                if protein_coding_bool == True:
                    if locus_group != 'protein-coding gene':
                        continue
                hgnc_entrez[entrez_id] = {'symbol': symbol, 'name': name, 'ensembl_gene_id': ensembl_gene_id}
                hgnc_symbol[symbol] = {'entrez_id': entrez_id, 'name': name, 'ensembl_gene_id': ensembl_gene_id}

    return {'hgnc_symbol_ALL': hgnc_symbol_ALL, 'hgnc_entrez_ALL': hgnc_entrez_ALL,
            'prev_symbols_ALL': prev_symbols_ALL, 'alias_symbols_ALL': alias_symbols_ALL,
            'human_entrez_ids': human_entrez,
            'hgnc_entrez': hgnc_entrez, 'hgnc_symbol': hgnc_symbol}