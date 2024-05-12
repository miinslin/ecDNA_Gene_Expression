# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd
import pickle
from boruta_py import BorutaPy
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import binom
from statsmodels.stats.multitest import fdrcorrection
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
training_testing_dir = os.path.join(base_dir, 'data', matrix_type, 'Training_Testing_datasets')

num_runs = 200

datasets_n = sorted([int(x.split('_')[1]) for x in os.listdir(training_testing_dir)])
start_run = datasets_n[0]

TCGA_matrix_samples_fn = os.path.join(os.path.dirname(TCGA_matrix_fn),\
                                      os.path.basename(TCGA_matrix_fn).replace('_matrix.txt','_matrix_samples.txt'))
if not os.path.isfile(TCGA_matrix_samples_fn):
    err_msg = TCGA_matrix_samples_fn+' does not exist!'
    raise ValueError(err_msg)

Boruta_runs_dir = os.path.join(base_dir, 'data', matrix_type, 'Boruta_Trials')
if not os.path.isdir(Boruta_runs_dir):
    os.makedirs(Boruta_runs_dir)

leave_out_percentage = 0.20
gene_column = 'Gene_name|Gene_id'
use_log = True
#RandomForestClassifier params
parallel = 2
class_weight='balanced_subsample'
max_depth=7
#boruta params
alpha = 0.05
ppf_cutoff = 0.005
stagnant_count_max = 5
tentative_count_min = 50    
boruta_iter = 400

#read gene expression matrix
parse_matrix = read_expression_matrix(TCGA_matrix_fn)
TCGA_combined_exp_matrix = parse_matrix['TCGA_combined_exp_matrix']
sample_classification = parse_matrix['sample_classification']
TCGA_sample_to_tumor = parse_matrix['TCGA_sample_to_tumor']
TCGA_sample_to_class = parse_matrix['TCGA_sample_to_class']
TCGA_all_tumor_samples = parse_matrix['TCGA_all_tumor_samples']

TCGA_all_gene = sorted(TCGA_combined_exp_matrix.keys())
print('#genes: ',len(TCGA_all_gene))

TCGA_all_tumor_samples = sorted(TCGA_all_tumor_samples)
print('#samples: ', len(TCGA_all_tumor_samples))

#run boruta
for b_run in range(start_run, num_runs, 1):
    data_run = 'data_'+str(b_run)
    
    # feed training dataset to boruta (default rf params in boruta script)
    training_samples_fn = os.path.join(training_testing_dir, \
                                       data_run, \
                                       'Training_p'+str(100-int(100*leave_out_percentage))+'_TCGA_samples.txt')
    
    output_dir = os.path.join(Boruta_runs_dir, data_run)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        
    BorutaPy_features_fn = os.path.join(output_dir, data_run+'_BorutaPy_features.txt')
    
    #logfile
    logfile_fn = os.path.join(output_dir, data_run+'.log')
    if os.path.isfile(logfile_fn):
        print('logfile_fn already exists...skipping', data_run)
        continue
   
    logfile = os.open(logfile_fn, os.O_WRONLY|os.O_CREAT)
    line = '\n'.join([data_run,
                      os.path.basename(__file__),
                      TCGA_matrix_fn,
                      TCGA_matrix_samples_fn,
                      training_samples_fn,
                      BorutaPy_features_fn,
                      'class_weight: '+class_weight,
                      'max_depth: '+str(max_depth),
                      'boruta_iter: '+str(boruta_iter),
                      'stagnant_count_max: '+str(stagnant_count_max),
                      'tentative_count_min: '+str(tentative_count_min),
                      'ppf_cutoff: '+str(ppf_cutoff),
                      'alpha: '+str(alpha)+'\n'])
    os.write(logfile, str.encode(line+'\n'))
    
    ############## split training testing datasets ##############
    ecDNA_samples = sample_classification['ecDNA']
    NoSCNA_samples = sample_classification['Non_ecDNA']
    os.write(logfile, str.encode('ecDNA: '+str(len(ecDNA_samples))+' Non-ecDNA: '+str(len(NoSCNA_samples))+'\n'))
    
    os.write(logfile, str.encode('##Reading training, testing data sets...'+'\n'))  

    training_ecDNA_samples = set()
    training_NoSCNA_samples = set()
    with open(training_samples_fn,'r') as infile:
        for i,line in enumerate(infile):
            sp = line.rstrip().split('\t')
            if i ==0:
                header = {x:xi for xi,x in enumerate(sp)}
            else:
                sampleID = sp[header['sampleID']]
                classification = sp[header['classification']] 
                if classification == 'Non_ecDNA':
                    training_NoSCNA_samples.add(sampleID)
                if classification == 'ecDNA':
                    training_ecDNA_samples.add(sampleID)

    training_samples = set()
    training_samples.update(training_ecDNA_samples)
    training_samples.update(training_NoSCNA_samples)

    os.write(logfile, str.encode('remaining '+str(100-100.0*leave_out_percentage)+'% data: '+str(len(training_samples))+'\n'))
    os.write(logfile, str.encode('ecDNA:'+str(len(training_ecDNA_samples))+' Non-ecDNA:'+str(len(training_NoSCNA_samples))+'\n'))
    
    training_ecDNA_samples = list(training_ecDNA_samples)
    training_NoSCNA_samples = list(training_NoSCNA_samples)
    training_samples = list(training_samples)
    
    ############## Boruta feature selection (on Training data) ##############
    os.write(logfile, str.encode('## Boruta on training samples...'+'\n'))
    
    X_orig_fn = os.path.join(output_dir, data_run+'_X_train_dict.pkl')
    X_train_fn = os.path.join(output_dir, data_run+'_X_train.pkl')
    X_train_fn_columns_fn = os.path.join(output_dir, data_run+'_X_train_column_names.pkl')
    
    if not os.path.isfile(X_train_fn):
        X = {g:[] for g in TCGA_all_gene}
        for gene in TCGA_all_gene:
            for sample in training_samples:
                value = TCGA_combined_exp_matrix[gene][sample]
                X[gene].append(value)
    
        outfile = open(X_orig_fn,'wb')
        pickle.dump(X, outfile)
        outfile.close()
 
        X = pd.DataFrame(X)
        X.to_pickle(X_train_fn)
        X_columns = X.columns.to_list()

        outfile = open(X_train_fn_columns_fn,'wb')
        pickle.dump(X_columns, outfile)
        outfile.close()
    
    X = pd.read_pickle(X_train_fn)    
    gene_name_list = pd.read_pickle(X_train_fn_columns_fn)
    
    y = np.array([1 if TCGA_sample_to_class[x]=='ecDNA' else 0 for x in training_samples])
    
    boruta_fn = os.path.join(output_dir, data_run+'_boruta.log') 
    
    if not os.path.isfile(boruta_fn):
        print('###Running boruta...')
        forest = RandomForestClassifier(class_weight=class_weight,\
                                        max_depth=max_depth,\
                                        random_state=1, n_jobs=parallel)

        boruta = BorutaPy(estimator = forest,\
                          n_estimators='auto',\
                          max_iter = boruta_iter, verbose=2,\
                          logfile_fn = boruta_fn, pre_stop=True, \
                          stagnant_count_max = stagnant_count_max,\
                          tentative_count_min = tentative_count_min,\
                          random_state=1)
            
        boruta.fit(np.array(X),y)

    print('###Parsing boruta log file...')
    boruta_features = []
    green_area = []
    blue_area = []
    iteration_history = {}
    active_feature_indexes = ''
    Confirmed = ''
    Tentative = ''
    Rejected = ''
    bonferroni_corrected_threshold = ''
    uncorr_pvalues_accept = []
    uncorr_pvalues_reject = []
    bh_adjusted_pvalues_accept = []
    bh_adjusted_pvalues_reject = []
    # dec_reg[active_features[to_accept]] = 1
    # dec_reg[active_features[to_reject]] = -1
    #1 uncorr pvalue < bonferroni correction 
    #2 adjusted pvalue
    #3 if both true, update status
    iteration = ''
    max_interations = ''
    last_iteration_values = []
    hit_reg = []
    dec_reg = []
    Last_iteration = ''
    
    with open(boruta_fn,'r') as infile:
        for line in infile:
            if '#max_iterations:' in line:
                max_interations = int(line.strip().split(':')[1].strip())
                
            if '#iteration:' in line:
                iteration = int(line.strip().split('#iteration: ')[1])
                iteration_history[iteration] = {}
            # all >=0 , confirmed + tentative
            if 'active_feature_indexes (dec_reg)' in line:
                active_feature_indexes = [int(x.strip()) for x in line.strip().split(':')[1].split(';')]
                iteration_history[iteration]['active_feature_indexes'] = active_feature_indexes 
            
            if 'bonferroni_corrected_threshold:' in line:
                bonferroni_corrected_threshold = float(line.strip().split(':')[1].strip())
                iteration_history[iteration]['bonferroni_corrected_threshold'] = bonferroni_corrected_threshold 
            if 'uncorr_pvalues to_accept' in line:
                uncorr_pvalues_accept = [float(x.strip()) for x in line.strip().split(':')[1].split(';')]
                iteration_history[iteration]['uncorr_pvalues_accept'] = uncorr_pvalues_accept 
            if 'uncorr_pvalues to_reject' in line:
                uncorr_pvalues_reject = [float(x.strip()) for x in line.strip().split(':')[1].split(';')]
                iteration_history[iteration]['uncorr_pvalues_reject'] = uncorr_pvalues_reject 
            
            if 'bh_adjusted_pvalues to_accept' in line:
                bh_adjusted_pvalues_accept = [float(x.strip()) for x in line.strip().split(':')[1].split(';')]
                iteration_history[iteration]['bh_adjusted_pvalues_accept'] = bh_adjusted_pvalues_accept 
            if 'bh_adjusted_pvalues to_reject' in line:
                bh_adjusted_pvalues_reject = [float(x.strip()) for x in line.strip().split(':')[1].split(';')]
                iteration_history[iteration]['bh_adjusted_pvalues_reject'] = bh_adjusted_pvalues_reject 
    
            if 'Iteration:' in line:
                Last_iteration = int(line.strip().split(':')[1].strip().split("/")[0].strip())
            
            if 'Confirmed:' in line:
                Confirmed = int(line.strip().split(':')[1].strip())
                iteration_history[iteration]['Confirmed'] = Confirmed 
            if 'Tentative:' in line:
                Tentative = int(line.strip().split(':')[1].strip())
                iteration_history[iteration]['Tentative'] = Tentative 
            if 'Rejected:' in line:
                Rejected = int(line.strip().split(':')[1].strip())
                iteration_history[iteration]['Rejected'] = Rejected
            
            if '#hit_reg (counts how many times' in line:
                hit_reg = [int(x.strip()) for x in line.strip().split(':')[1].split(';')]
            
            if '#dec_reg (holds the decision about each feature' in line:
                dec_reg = [int(x.strip()) for x in line.strip().split(':')[1].split(';')]
            
            if '#green_area' in line:
                green_area = [xi for xi,x in enumerate(line.strip().split(':')[1].split(';')) if x.strip()=='True']
    
            if '#blue_area' in line:
                blue_area = [xi for xi,x in enumerate(line.strip().split(':')[1].split(';')) if x.strip()=='True']
    
    print('features in the green area:', len(green_area))
    print('features in the blue area:', len(blue_area))
    os.write(logfile, str.encode('features in the green area:'+str(len(green_area))+'\n'))
    os.write(logfile, str.encode('features in the blue area:'+str(len(blue_area))+'\n'))

    print('last iteration: ',Last_iteration, ' out of ', max_interations)
    os.write(logfile, str.encode('last iteration: '+str(Last_iteration)+'\n'))
    
    bluegreen_area = set(green_area[0:])
    bluegreen_area.update(blue_area)
    
    n,p = Last_iteration, 0.5
    
    reject_cutoff = binom.ppf(ppf_cutoff, n, p)
    boruta_features = {}
    for xi,x in enumerate(hit_reg):
        if x>reject_cutoff:
            g = gene_name_list[xi]
            if xi in green_area:
                boruta_features[g] = {'area':'Confirmed','hits':x}
            elif xi in blue_area:
                boruta_features[g] = {'area':'Tentative','hits':x}
    print('boruta_features: ',len(boruta_features))
        
    bonferroni_corrected_pvalue = alpha/float(n)
    
    uncorrected_p_values = {}
    for g in boruta_features:
        hits = boruta_features[g]['hits']
        uncorrected_p_values[g] = binom.sf(hits - 1, n, p)        
    uncorrected_p_values = uncorrected_p_values.items()
    
    fdr_correction_ = fdrcorrection([x[1] for x in uncorrected_p_values], alpha=alpha)
    
    to_accept = np.array(fdr_correction_[0])
    to_accept2 = np.array([x[1] for x in uncorrected_p_values]) <= bonferroni_corrected_pvalue
    to_accept *= to_accept2
    to_accept_list = to_accept.tolist()

    corrected_p_values = fdr_correction_[1]
    computed_boruta_pvalues = {}
    for xi,x in enumerate(uncorrected_p_values):
        g = x[0]
        uncorr_p = x[1]
        two_step_corr = str(to_accept_list[xi])
        computed_boruta_pvalues[g] = {'uncorr_p':uncorr_p,'two_step_corr': two_step_corr}
    
    Boruta_genelist = set()
    with open(BorutaPy_features_fn,'w') as outfile:
        outfile.write('#'+'\t'.join([gene_column, 'decision', 'hits', 'uncorr_pvalue','two_step_corr'])+'\n')
        for gg in boruta_features:
            p1 = repr(computed_boruta_pvalues[gg]['uncorr_p'])
            two_step_corr = computed_boruta_pvalues[gg]['two_step_corr']
            decision = boruta_features[gg]['area']
            hit_n = repr(boruta_features[gg]['hits'])
            outfile.write('\t'.join([gg, decision, hit_n, p1, two_step_corr])+'\n')
            Boruta_genelist.add(gg)
    
    Boruta_genelist = sorted(Boruta_genelist)
    print(str(len(Boruta_genelist))+' Boruta_genelist')
    os.write(logfile, str.encode(str(len(Boruta_genelist))+' Boruta_genelist'+'\n'))
    
    os.close(logfile)

    #remove pkl files
    os.remove(X_orig_fn)
    os.remove(X_train_fn)
    os.remove(X_train_fn_columns_fn)   