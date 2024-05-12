# -*- coding: utf-8 -*-
import os
import argparse
import random

base_dir = os.path.abspath(os.path.join(__file__ ,"../.."))

parser = argparse.ArgumentParser()
parser.add_argument("-D", "--matrix", dest = "matrix",\
                    default = os.path.join(base_dir, 'data', \
                    'RNAseq_expression_matrix',\
                    'AC049_TCGA_KNN_imputed_min_ecdna_3_matrix.txt'),\
                    help="Gene expression matrix path.")
args = parser.parse_args()

TCGA_matrix_fn = os.path.normpath(args.matrix)
num_runs = 200
leave_out_percentage = 0.20

matrix_type = os.path.basename(TCGA_matrix_fn).replace('.txt','')
training_testing_dir = os.path.join(base_dir, 'data', matrix_type, 'Training_Testing_datasets')
if not os.path.isdir(training_testing_dir):
    print('Creating ', training_testing_dir)
    os.makedirs(training_testing_dir)

TCGA_matrix_samples_fn = os.path.join(os.path.dirname(TCGA_matrix_fn),\
                                      os.path.basename(TCGA_matrix_fn).replace('_matrix.txt','_matrix_samples.txt'))
if not os.path.isfile(TCGA_matrix_samples_fn):
    err_msg = TCGA_matrix_samples_fn+' does not exist!'
    raise ValueError(err_msg)
 
sample_classification = {'Non_ecDNA':[],'ecDNA':[]}
TCGA_sample_to_tumor = {}
TCGA_sample_to_class = {}
TCGA_all_tumor_samples = []
with open(TCGA_matrix_samples_fn,'r') as infile:
    for i,line in enumerate(infile):
        sp = line.rstrip('\n').split('\t')
        if i==0:
            colidx = {x:xi for xi,x in enumerate(sp)}
        else:
            ecDNA_status = sp[colidx['#ecDNA_status']]
            tumor = sp[colidx['Tumor']]
            sample_barcode = sp[colidx['Sample_barcode']]
            sample_classification[ecDNA_status].append(sample_barcode)
            TCGA_sample_to_tumor[sample_barcode] = tumor
            TCGA_sample_to_class[sample_barcode] = ecDNA_status
            TCGA_all_tumor_samples.append(sample_barcode)

TCGA_all_tumor_samples = sorted(TCGA_all_tumor_samples)
print('All samples: ', len(TCGA_all_tumor_samples))
print('ecDNA(+): ', len(sample_classification['ecDNA']))
print('ecDNA(-): ', len(sample_classification['Non_ecDNA']))

for b_run in range(num_runs):
    
    data_run = 'data_'+str(b_run)
        
    output_dir = os.path.join(training_testing_dir, data_run)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    training_samples_fn = os.path.join(output_dir,\
                         'Training_p'+str(100-int(100*leave_out_percentage))+'_TCGA_samples.txt')
    testing_samples_fn = os.path.join(output_dir,\
                         'Testing_p'+str(int(100*leave_out_percentage))+'_TCGA_samples.txt')
       
    logfile_fn = os.path.join(output_dir, data_run+'_leave_out_'+str(leave_out_percentage)+'p.log')
    
    logfile = os.open(logfile_fn, os.O_WRONLY|os.O_CREAT)
    line = '\n'.join([data_run,
                      os.path.basename(__file__),
                      'TCGA_all_tumor_samples: '+str(len(TCGA_all_tumor_samples)),
                      testing_samples_fn,
                      training_samples_fn,
                      'leave out '+str(100.0*leave_out_percentage)+'% data for testing\n'])
    os.write(logfile, str.encode(line+'\n'))
    
    ############## split training testing datasets ##############
    ecDNA_samples = sample_classification['ecDNA']
    NoSCNA_samples = sample_classification['Non_ecDNA']
    os.write(logfile, str.encode('ecDNA: '+str(len(ecDNA_samples))+\
                                 ' Non-ecDNA: '+str(len(NoSCNA_samples))+'\n'))
    
    #randomly select testing samples
    testing_ecDNA_samples = random.sample(ecDNA_samples,int(leave_out_percentage*len(ecDNA_samples)))
    testing_NoSCNA_samples = random.sample(NoSCNA_samples,int(leave_out_percentage*len(NoSCNA_samples)))
    testing_set = set()
    testing_set.update(testing_ecDNA_samples)
    testing_set.update(testing_NoSCNA_samples)
    
    os.write(logfile, str.encode(str(100.0*leave_out_percentage)+\
                                 '% data left out for testing: '+str(len(testing_set))+'\n'))
    os.write(logfile,str.encode('ecDNA: '+str(len(testing_ecDNA_samples))+\
                                ' Non-ecDNA: '+str(len(testing_NoSCNA_samples))+'\n'))
    
    with open(testing_samples_fn,'w') as outfile:
        outfile.write('\t'.join(['sampleID','tumor','classification'])+'\n')
        for x in testing_set:
            outfile.write('\t'.join([x, TCGA_sample_to_tumor[x], TCGA_sample_to_class[x]])+'\n')
    
    #remaining samples as training
    training_ecDNA_samples = set(ecDNA_samples)-testing_set
    training_NoSCNA_samples = set(NoSCNA_samples)-testing_set
    training_samples = set()
    training_samples.update(training_ecDNA_samples)
    training_samples.update(training_NoSCNA_samples)
    os.write(logfile, str.encode('remaining '+str(100-100.0*leave_out_percentage)+\
                                 '% data: '+str(len(training_samples))+'\n'))
    os.write(logfile, str.encode('ecDNA:'+str(len(training_ecDNA_samples))+\
                                 ' Non-ecDNA:'+str(len(training_NoSCNA_samples))+'\n'))
    
    with open(training_samples_fn,'w') as outfile:
        outfile.write('\t'.join(['sampleID','tumor','classification'])+'\n')
        for x in training_samples:
            outfile.write('\t'.join([x, TCGA_sample_to_tumor[x], TCGA_sample_to_class[x]])+'\n')
    
    os.close(logfile)