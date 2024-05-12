import os
import numpy as np

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))
reference_dir = os.path.join(base_dir, 'data', 'reference')
sample_info_fn = os.path.join(reference_dir, 'TCGA', 'AC049_cbioportal_sample_info.txt')

def retrieve_sample_info(min_ecdna=0, remove_metastatic=True):
    cbioportal_sample_info = {}
    AC049_counts = {}
    class_counts = {'ecDNA': set(), 'Non_ecDNA': set()}
    with open(sample_info_fn, 'r') as infile:
        for i, line in enumerate(infile):
            sp = line.rstrip('\n').split('\t')
            if i == 0:
                colidx = {x: xi for xi, x in enumerate(sp)}
            else:
                ecDNA_status = sp[colidx['#ecDNA_status']]  # Non_ecDNA, ecDNA
                Tumor = sp[colidx['Tumor']]
                Sample_barcode = sp[colidx['Sample_barcode']]
                SampleType = sp[colidx['SampleType']]  # Primary Solid Tumor, Metastatic
                if remove_metastatic == True:
                    if SampleType == 'Metastatic':
                        continue
                full_barcode = sp[colidx['full_barcode']]
                Center = sp[colidx['Center']]
                platform = sp[colidx['platform']]
                Adjustment = sp[colidx['Adjustment']]
                class_counts[ecDNA_status].add(Sample_barcode)
                cbioportal_sample_info[Sample_barcode] = {'ecDNA_status': ecDNA_status,
                                                          'Tumor': Tumor,
                                                          'full_barcode': full_barcode,
                                                          'SampleType': SampleType,
                                                          'Center': Center,
                                                          'platform': platform,
                                                          'Adjustment': Adjustment}
                if Tumor not in AC049_counts:
                    AC049_counts[Tumor] = {}
                if ecDNA_status not in AC049_counts[Tumor]:
                    AC049_counts[Tumor][ecDNA_status] = set()
                AC049_counts[Tumor][ecDNA_status].add(Sample_barcode)

    TCGA_tumors = AC049_counts.keys()
    samples_to_remove = set()
    if min_ecdna > 0:
        # remove samples from tumors with less than min_ecdna ecDNA(+) samples
        TCGA_tumors = set()
        for Tumor in AC049_counts:
            if 'ecDNA' not in AC049_counts[Tumor]:
                samples_to_remove.update(AC049_counts[Tumor]['Non_ecDNA'])
                class_counts['Non_ecDNA'] = class_counts['Non_ecDNA'] - AC049_counts[Tumor]['Non_ecDNA']
                continue
            ecDNA_samples = len(AC049_counts[Tumor]['ecDNA'])
            if ecDNA_samples < min_ecdna:
                samples_to_remove.update(AC049_counts[Tumor]['ecDNA'])
                class_counts['ecDNA'] = class_counts['ecDNA'] - AC049_counts[Tumor]['ecDNA']
                samples_to_remove.update(AC049_counts[Tumor]['Non_ecDNA'])
                class_counts['Non_ecDNA'] = class_counts['Non_ecDNA'] - AC049_counts[Tumor]['Non_ecDNA']
                continue
            TCGA_tumors.add(Tumor)

    for Sample_barcode in samples_to_remove:
        del cbioportal_sample_info[Sample_barcode]

    TCGA_tumors = sorted(TCGA_tumors)
    print('\nTCGA_tumors: ', len(TCGA_tumors))
    print('Total # samples: ', len(cbioportal_sample_info))
    print('ecDNA(+): ', len(class_counts['ecDNA']))
    print('ecDNA(-): ', len(class_counts['Non_ecDNA']))

    return {'sample_info': cbioportal_sample_info,
            'class_counts': class_counts,
            'TCGA_tumors': TCGA_tumors}


def read_expression_matrix(fn, key=''):
    use_log = True
    gene_column = 'Gene_name|Gene_id'

    TCGA_combined_exp_matrix = {}
    sample_classification = {'Non_ecDNA': [], 'ecDNA': []}
    TCGA_sample_to_tumor = {}
    TCGA_sample_to_class = {}
    TCGA_all_tumor_samples = []
    with open(fn, 'r') as infile:
        for i, line in enumerate(infile):
            sp = line.rstrip('\n').split('\t')
            if i == 0:
                colidx = {x: xi for xi, x in enumerate(sp)}
                # Gene_name|Gene_id	ecDNA|BLCA|TCGA-BL-A0C8-01
                samples = []
                column_ids = sp[1:]
                for entry in column_ids:
                    ecDNA_status, tumor, sample_barcode = entry.split('|')
                    samples.append(entry)
                    sample_classification[ecDNA_status].append(sample_barcode)
                    TCGA_sample_to_tumor[sample_barcode] = tumor
                    TCGA_sample_to_class[sample_barcode] = ecDNA_status
                    TCGA_all_tumor_samples.append(sample_barcode)
            else:
                gg = sp[colidx[gene_column]]
                if key == 'entrez':
                    gg = gg.split('|')[1]

                TCGA_combined_exp_matrix[gg] = {}
                for entry in samples:
                    ecDNA_status, tumor, sample_barcode = entry.split('|')
                    value = sp[colidx[entry]].replace('"', '').replace("'", '')
                    if use_log == True:
                        value = np.log2(float(value) + 1)
                    else:
                        value = float(value)
                    TCGA_combined_exp_matrix[gg][sample_barcode] = value

    return {'TCGA_combined_exp_matrix': TCGA_combined_exp_matrix,
            'sample_classification': sample_classification,
            'TCGA_sample_to_tumor': TCGA_sample_to_tumor,
            'TCGA_sample_to_class': TCGA_sample_to_class,
            'TCGA_all_tumor_samples': TCGA_all_tumor_samples}
