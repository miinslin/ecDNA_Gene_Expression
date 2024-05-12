# -*- coding: utf-8 -*-
import os
import random
import numpy as np
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
#stratified k fold split data
from sklearn.model_selection import StratifiedKFold
#grid search cv for hyper parameter optimization
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
#scores
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score

leave_out_percentage = 0.2
parallel = 3
randomsearch_iter = 50

gene_column = '#Gene_name|Gene_id'
scoring_metric = 'f1'
k_folds = 5

def validate_genes(training_testing_dir, model_performance_dir, TCGA_data_sub, genelist, num_runs, gene_type, num_genes, TCGA_combined_exp_matrix, TCGA_all_tumor_samples):

    for b_run in range(0, num_runs, 1):
    
        data_run = 'data_'+str(b_run)
        print(data_run)    
        
        output_dir = os.path.join(model_performance_dir, data_run)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
            
        model_important_features_fn = os.path.join(output_dir, data_run+'_model_important_features.txt')
        if os.path.isfile(model_important_features_fn):
            print('model_important_features_fn already exists...skipping', data_run)
            continue
        
        if gene_type == 'random':
            TCGA_all_gene = sorted(TCGA_combined_exp_matrix.keys())
            #choose random genes, write to file
            genelist = random.sample(TCGA_all_gene, num_genes)
            genelist = sorted(genelist)
            
            random_gene_fn = os.path.join(output_dir, data_run+'_random_genes.txt')
            with open(random_gene_fn,'w') as outfile:
                outfile.write('\n'.join(genelist))
            
            TCGA_data_sub = {}
            for sample in TCGA_all_tumor_samples:
                values = [TCGA_combined_exp_matrix[g][sample] for g in genelist]
                TCGA_data_sub[sample] = values
        
        training_samples_fn = os.path.join(training_testing_dir, data_run,\
                                           'Training_p'+str(100-int(100*leave_out_percentage))+'_TCGA_samples.txt')
        testing_samples_fn = os.path.join(training_testing_dir, data_run,\
                                          'Testing_p'+str(int(100*leave_out_percentage))+'_TCGA_samples.txt')
        #logfile
        logfile_fn = os.path.join(output_dir, data_run+'.log')
        
        logfile = os.open(logfile_fn, os.O_WRONLY|os.O_CREAT)
        line = '\n'.join([data_run,
                          os.path.basename(__file__), 
                          testing_samples_fn,
                          training_samples_fn,
                          'genelist: '+str(len(genelist)),
                          model_important_features_fn,
                          'randomsearch_iter: '+ str(randomsearch_iter),
                          'runcount: '+str(num_runs),                
                          'k_folds: '+str(k_folds)+'\n'])
        os.write(logfile, str.encode(line+'\n'))
    
        ############### Testing, Training Samples ##############
        os.write(logfile, str.encode('##Reading training, testing data sets...'+'\n'))      
       
        TCGA_sample_to_class = {}
       
        testing_set = []
        with open(testing_samples_fn,'r') as infile:
            for i,line in enumerate(infile):
                sp = line.rstrip().split('\t')
                if i ==0:
                    header = {x:xi for xi,x in enumerate(sp)}
                else:
                    sampleID = sp[header['sampleID']]
                    classification = sp[header['classification']] 
                    testing_set.append(sampleID)
                    TCGA_sample_to_class[sampleID] = classification
        os.write(logfile, str.encode(str(100.0*leave_out_percentage)+'% data left out for testing: '+str(len(testing_set))+'\n'))
    
        training_samples = []
        with open(training_samples_fn,'r') as infile:
            for i,line in enumerate(infile):
                sp = line.rstrip().split('\t')
                if i ==0:
                    header = {x:xi for xi,x in enumerate(sp)}
                else:
                    sampleID = sp[header['sampleID']]
                    classification = sp[header['classification']] 
                    training_samples.append(sampleID)
                    TCGA_sample_to_class[sampleID] = classification
        os.write(logfile, str.encode('remaining '+str(100-100.0*leave_out_percentage)+'% data: '+str(len(training_samples))+'\n'))
    
        ############### Training, Boruta selected genes ##############
        
        X_train = [TCGA_data_sub[s] for s in training_samples]
        X_train = np.array([np.array(xi) for xi in X_train])
        #ecDNA 1, noSCNA 0
        y_train = np.array([1 if TCGA_sample_to_class[x]=='ecDNA' else 0 for x in training_samples])
    
        X_test = [TCGA_data_sub[s] for s in testing_set]
        X_test = np.array([np.array(xi) for xi in X_test])
        y_test = np.array([1 if TCGA_sample_to_class[x]=='ecDNA' else 0 for x in testing_set])
        
        
        randomized_scv_stats = {}  
        grid_scv_stats = {}
        
        #RandomizedSearchCV 
        os.write(logfile, str.encode('##RandomizedSearchCV...'+'\n'))   
        print('\n##RandomizedSearchCV')
        # Number of trees in random forest
        n_estimators = [int(x) for x in np.linspace(100, 2000, num = 10)]
        # Maximum number of levels in tree
        max_depth = [int(x) for x in np.linspace(10, 100, num = 10)]
        max_depth.append(None)
        # # Minimum number of samples required to split a node
        min_samples_split = [2, 5, 10]
        # # Minimum number of samples required at each leaf node
        min_samples_leaf = [1, 2, 4]
        # # Method of selecting samples for training each tree
        bootstrap = [True]
        # Create the random grid
        random_params = {'n_estimators': n_estimators,
                       'max_depth': max_depth,
                       'min_samples_split': min_samples_split,
                       'min_samples_leaf': min_samples_leaf,
                       'bootstrap': bootstrap}    
        
        rf = RandomForestClassifier(class_weight='balanced_subsample', n_jobs=parallel, random_state=1)
        rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_params,\
                                       scoring = scoring_metric,\
                                       n_iter = randomsearch_iter, cv = k_folds, \
                                       verbose=2, random_state=1, n_jobs = parallel)
            
        RandomizedSearchCV_result = rf_random.fit(X_train, y_train)
        RandomizedSearchCV_best_estimator = RandomizedSearchCV_result.best_estimator_
        RandomizedSearchCV_best_params = RandomizedSearchCV_result.best_params_  
        RandomizedSearchCV_model_feature_importances = RandomizedSearchCV_result.best_estimator_.feature_importances_
            
        random_yhat = RandomizedSearchCV_best_estimator.predict(X_test)
        
        acc = accuracy_score(y_test, random_yhat)
        roc_auc = roc_auc_score(y_test, random_yhat)
        precision = precision_score(y_test, random_yhat)
        recall = recall_score(y_test, random_yhat)
        f1 = f1_score(y_test, random_yhat)
        
        randomized_scv_stats['f1'] = float(f1)
        randomized_scv_stats['recall'] = float(recall)
        randomized_scv_stats['precision'] = float(precision)
        randomized_scv_stats['acc'] = float(acc)
        randomized_scv_stats['roc_auc'] = float(roc_auc)
            
        tn, fp, fn, tp = confusion_matrix(y_test, random_yhat).ravel()
        
        print('>acc=%.3f, rocauc=%.3f, precision=%.3f, recall=%.3f, f1=%.3f' % (acc, roc_auc, precision, recall, f1))
        print('>TN=%.3f, FP=%.3f, FN=%.3f, TP=%.3f' % (tn, fp, fn, tp)+'\n')
    
        os.write(logfile,str.encode('best_estimator_: '+str(RandomizedSearchCV_best_estimator)+'\n'))
        os.write(logfile,str.encode('best_params_: '+str(RandomizedSearchCV_best_params)+'\n'))
        os.write(logfile,str.encode('RF train accuracy: '+repr(round(RandomizedSearchCV_result.score(X_train, y_train),3))+'\n'))
        os.write(logfile,str.encode('RF test accuracy: '+repr(round(RandomizedSearchCV_result.score(X_test, y_test),3))+'\n\n'))
        os.write(logfile,str.encode('>acc=%.3f, rocauc=%.3f, precision=%.3f, recall=%.3f, f1=%.3f' % (acc, roc_auc, precision, recall, f1)+'\n'))
        os.write(logfile,str.encode('>TN=%.3f, FP=%.3f, FN=%.3f, TP=%.3f' % (tn, fp, fn, tp)+'\n'))
    
        #GridSearchCV 
        os.write(logfile, str.encode('##GridSearchCV...'+'\n'))   
        print('\n##GridSearchCV')
        
        #n estimators
        rand_n_est = RandomizedSearchCV_best_params['n_estimators']
        rand_n_est = [rand_n_est-150, rand_n_est, rand_n_est+150]
        rand_n_est = [x for x in rand_n_est if x>=100]
        
        #max depth
        rand_max_dep = RandomizedSearchCV_best_params['max_depth']
        if rand_max_dep!=None:
            rand_max_dep = [rand_max_dep-10, rand_max_dep, rand_max_dep+10]
            rand_max_dep = [x for x in rand_max_dep if x>=10]
            if len(rand_max_dep)!=3:
                add_V = rand_max_dep[-1]+10
                rand_max_dep.append(add_V)
        else:
            rand_max_dep = [rand_max_dep]
            
        gridsearch_params = {'n_estimators': rand_n_est,
                   'max_depth': rand_max_dep,
                   'min_samples_split': min_samples_split,
                   'min_samples_leaf': min_samples_leaf,
                   'bootstrap': bootstrap}
        
        rf = RandomForestClassifier(class_weight='balanced_subsample', n_jobs=parallel, random_state=1)
        cv_ = StratifiedKFold(n_splits=k_folds, shuffle=True, random_state=1)
        grid_rf = GridSearchCV(estimator = rf, param_grid = gridsearch_params,\
                              scoring = scoring_metric,\
                              cv=cv_,\
                              n_jobs=parallel, verbose=1)
            
        GridSearchCV_result = grid_rf.fit(X_train, y_train)
        GridSearchCV_best_estimator = GridSearchCV_result.best_estimator_
        GridSearchCV_best_params = GridSearchCV_result.best_params_
        GridSearchCV_model_feature_importances = GridSearchCV_result.best_estimator_.feature_importances_
        
        grid_yhat = GridSearchCV_best_estimator.predict(X_test)
        
        acc = accuracy_score(y_test, grid_yhat)
        roc_auc = roc_auc_score(y_test, grid_yhat)
        precision = precision_score(y_test, grid_yhat)
        recall = recall_score(y_test, grid_yhat)
        f1 = f1_score(y_test, grid_yhat)
    
        grid_scv_stats['f1'] = float(f1)
        grid_scv_stats['recall'] = float(recall)
        grid_scv_stats['precision'] = float(precision)
        grid_scv_stats['acc'] = float(acc)
        grid_scv_stats['roc_auc'] = float(roc_auc)
    
        tn, fp, fn, tp = confusion_matrix(y_test, grid_yhat).ravel()
        
        print('>acc=%.3f, rocauc=%.3f, precision=%.3f, recall=%.3f, f1=%.3f' % (acc, roc_auc, precision, recall, f1))
        print('>TN=%.3f, FP=%.3f, FN=%.3f, TP=%.3f' % (tn, fp, fn, tp)+'\n')
    
        os.write(logfile,str.encode('best_estimator_: '+str(GridSearchCV_best_estimator)+'\n'))
        os.write(logfile,str.encode('best_params_: '+str(GridSearchCV_best_params)+'\n'))
        os.write(logfile,str.encode('RF train accuracy: '+repr(round(GridSearchCV_result.score(X_train, y_train),3))+'\n'))
        os.write(logfile,str.encode('RF test accuracy: '+repr(round(GridSearchCV_result.score(X_test, y_test),3))+'\n\n'))
        os.write(logfile,str.encode('>acc=%.3f, rocauc=%.3f, precision=%.3f, recall=%.3f, f1=%.3f' % (acc, roc_auc, precision, recall, f1)+'\n'))
        os.write(logfile,str.encode('>TN=%.3f, FP=%.3f, FN=%.3f, TP=%.3f' % (tn, fp, fn, tp)+'\n'))

        feats = {} # a dict to hold feature_name: feature_importance
        if randomized_scv_stats['f1']>grid_scv_stats['f1']:
            for feature, importance in zip(genelist, RandomizedSearchCV_model_feature_importances):
                feats[feature] = importance #add the name/value pair 
        else:
            for feature, importance in zip(genelist, GridSearchCV_model_feature_importances):
                feats[feature] = importance #add the name/value pair               
            
        ranked_features = sorted(feats.items(),key = lambda x: x[1], reverse=True)   
        with open(model_important_features_fn,'w') as outfile:
            outfile.write('#'+'\t'.join(['Gene|GeneId', 'feature_importance'])+'\n')
            for entry in ranked_features:
                gg = entry[0]
                importance_value = entry[1]
                outfile.write('\t'.join([gg, repr(importance_value)])+'\n')    
    
        os.close(logfile)