#!/usr/bin/env python
# coding: utf-8

# In[1]:


#args =  ["classifier.py","/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "CRC_k6&CRC_k7", "kmer", "BatchCorrected",  "bin_crc_normal",1,0,20,1,1] 

#python MINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc "CRC_k6&CRC_k7" kmer BatchCorrected bin_crc_normal 1 0 20 1 1
# python MINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc "AGP_k6&CRC_k7" kmer BatchCorrected bin_crc_normal 1 0 20 1 1

# In[2]:


## data = np.random.rand(100,80)
# labels = np.random.binomial(size=100, n=1, p= 0.5)

# parameter_dict = {'n_estimators':[10,100],'criterion': ['entropy'],'min_samples_leaf': [10],'max_features':[0.3],'min_samples_split': [2, 5, 10],}
# clf = GridSearchCV(rf, parameter_dict,scoring="roc_auc")
# clf.fit(data, labels)
# clf.best_params_
# y_pred = clf.predict_proba(data)
# y_true = labels
# already_trained_auc = roc_auc_score(y_true = y_true, y_score = y_pred[:,1])
# already_trained_auc
# reducing n_estimators reduced AUC on train, increasing min_samples_leaf reduced AUC on train

import sys
import pandas as pd
import utils
import numpy as np
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.ensemble import RandomForestClassifier
from sklearn import model_selection 
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, cross_val_score
import statsmodels.formula.api as sm
from sklearn.metrics import roc_auc_score
from collections import Counter
from timeit import default_timer as timer
import pickle



args = sys.argv


# In[3]:


greater_folder = args[1] # what folder do you save your different datasets in
study_names = args[2].split("&")  # what is the name of the dataset (tells the program which folder to check)
data_type = args[3] # type of data. kmer vs OTU

prefix_name = args[4] # what is the prefix of the file name
column_of_interest = args[5] # what is the phenotype you are predicting (use the same name in the column of the metadata you want to predict), this programs reads from metadata.txt

norm_input = bool(int(args[6]))
map_with_accession = bool(int(args[7]))
num_pcs = 20
num_pcs = int(args[8])
special_name = args[9]

if len(args) > 10:
    label_pos_or_neg = int(args[10]) # do you want to treat CRC as positive class or negative class? 
    target_label = args[11] # phenotype representing positive class or negative class? eg. CRC eg. H
    print(target_label)
else:
    label_pos_or_neg = 1
    target_label = 1
use_domain_pheno = False # for when running raw to compare to domain pheno
data_folders = [greater_folder + "/data/" + study_name + "/" for study_name in study_names]   



# In[4]:


#########################################################################
###### COMMENTARY: load data from your k-mer matrix, load metadata ######
#########################################################################

feature_table_dict = utils.load_feature_table(data_folders,data_type = data_type)
metadata = pd.read_csv(data_folders[0] + "metadata.txt",delimiter="\t")

if norm_input:
    for d in range(len(study_names)):
        #feature_table_dict[d] = normalize(np.array(feature_table_dict[d].transpose()), axis = 1, norm = 'l1')
        temp = pd.DataFrame(normalize(feature_table_dict[d].transpose(), axis = 1, norm = 'l1').transpose())
        temp.index = feature_table_dict[d].index
        temp.columns = feature_table_dict[d].columns
        feature_table_dict[d] = temp


# In[ ]:





# In[ ]:


if "AGP" in study_names[0]:
    for d in range(len(study_names)):
        tissue_samples = metadata.index[metadata['body_habitat.x'] == "UBERON:feces"]

        feature_table_dict[d] = feature_table_dict[d][tissue_samples]
        metadata = metadata.loc[tissue_samples]


print(Counter(metadata[column_of_interest]))

print("pos label")
print(target_label)
if len(args) > 10:
    if label_pos_or_neg == 1:
        print("positive")
        metadata[column_of_interest] = utils.binarize_labels_mod(metadata[column_of_interest],none_labels = ["not applicable",float("Nan"),'not provided'],pos_labels =[target_label])
    elif label_pos_or_neg == 0:
        metadata[column_of_interest] = utils.binarize_labels_mod(metadata[column_of_interest],none_labels = ["not applicable",float("Nan"),'not provided'],neg_labels =[target_label])

print(Counter(metadata[column_of_interest]))


# remove na samples
print("first dimension")
print(metadata.shape)
#print(np.isnan(metadata[column_of_interest]))
print(metadata[column_of_interest][0:5])
non_nan_samples = metadata.index[np.invert(np.isnan(metadata[column_of_interest]))]
print("non na samples")
print(non_nan_samples[0:4])
metadata = metadata.loc[non_nan_samples]
print(metadata.shape)

for d in range(len(study_names)):
    feature_table_dict[d] = feature_table_dict[d][non_nan_samples]
    


#print(non_nan_samples)
#[x for x in [0,1,2]]
#########################################################################
###### COMMENTARY:  efining labels and binarize if not already     ######
#########################################################################



labels = metadata[column_of_interest]




# In[ ]:


def pca_regression(y,X):
    model = sm.OLS(y,X)
    results = model.fit()
    predictedValues = results.predict()
    residuals = y - predictedValues
    return(residuals)

def RF_grid_search(data,labels,param_dict):
    rf = RandomForestClassifier()
    clf = GridSearchCV(rf, param_dict,scoring="roc_auc")
    clf.fit(data, labels)

    print("Best parameters set found on development set:")
    best_params = clf.best_params_
    return clf,best_params

def RF_cv(data,labels,param_dict):
    clf = RandomForestClassifier(max_depth=5, random_state=0,n_estimators = param_dict['n_estimators'],            criterion = param_dict['criterion'],min_samples_leaf = param_dict['min_samples_leaf'],                           max_features = param_dict['max_features'])
    results = cross_val_score(clf,X=data,y=labels,scoring="roc_auc")
    return(results)

   


# In[ ]:


# outline
# for data_table:
    # test train split
    # for train
    # for PC number:
        # regress out PCs
        # for grid cell:
            # get accuracy:
        # get max accuracy
    # get max accuracy across PCs
    # for test
    # get the test accuracy with best parameters
    # save best grid and best PCs


# In[ ]:


# get PC scores
pc_table_dict = dict()
feature_table_np = dict()
labels_np = dict()
pca = PCA(n_components=num_pcs,svd_solver='randomized')
for d in range(len(study_names)):
    temp = feature_table_dict[d].transpose()
    pca.fit(temp)
    pc_table_dict[d] = pca.transform(temp)
    
    feature_table_np[d] = np.array(feature_table_dict[0])
    labels_np[d] = np.array(labels)
    
  


# In[ ]:


n_splits = 5
n_repeats = 1
rskf = model_selection.RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=123)

if "maxfeat10" in special_name:
    parameter_dict = {'n_estimators':[20],'criterion': ['entropy'],\
    'min_samples_leaf': [10],'max_features':[0.10],'min_samples_split': [5],'max_depth':[5]}
if "AGP" in special_name:
    parameter_dict = {'n_estimators':[20,100,1000],'criterion': ['entropy','gini'],\
    'min_samples_leaf': [10,15],'max_features':[0.10,0.30],'min_samples_split': [5],'max_depth':[5,10]}

# parameter_dict = {'n_estimators':[10,50],'criterion': ['entropy'],\
# 'min_samples_leaf': [2,5],'max_features':[0.3],'min_samples_split': [2, 5],}    


# In[ ]:





# In[ ]:


# Regress out PCs
all_datasets_dict =  dict()
# For each dataset (kmer size)

start = timer()


for d in range(len(study_names)): # range(1):#
    dataset_start= timer()
    
    results_dict = dict()
    
    X = feature_table_np[d].transpose()
    y = labels_np[d]
    pc_scores = pc_table_dict[d] # get pc scores
    na_mask = pd.isna(y)
    
    X = X[~na_mask,:]
    y = y[~na_mask]
    pc_scores = pc_scores[~na_mask,:]
    # for each test train split in 5 fold cross validation
    
    train_it = 0
    for train_index, test_index in rskf.split(X, y):
        
        test_train_start = timer()
        #print(train_index)
        #print(test_index)
        
        X_train, X_test = X[train_index,], X[test_index,]
        y_train, y_test = y[train_index], y[test_index]
        pc_scores_train, pc_scores_test =  pc_scores[train_index], pc_scores[test_index]
        
        if train_it == 0: 
            results_dict["number samples"] = []
        results_dict["number samples"].append(X_train.shape[0])
        # for each PC we regress out 
        for p in range(num_pcs): #range(3):#
            if train_it == 0: 
                results_dict["PC" + str(p)] = dict()
                results_dict["PC" + str(p)]['train_best_params'] = dict()
                results_dict["PC" + str(p)]['train_auc_trained'] = []
                results_dict["PC" + str(p)]['mean_train_cv_auc'] = []
                results_dict["PC" + str(p)]['mean_test_cv_auc'] = []
                results_dict["PC" + str(p)]['test_auc_trained'] = []
               
            if p == 0:
                X_train_corrected = X_train
                X_test_corrected = X_test
            else:
                X_train_corrected = pca_regression(X_train,pc_scores_train[:,0:p])
                X_test_corrected = pca_regression(X_test,pc_scores_test[:,0:p])
            
            # perform grid search on train
            best_train_model, best_params = RF_grid_search(X_train_corrected, y_train,parameter_dict)
            print("finished grid search: " + "train it " + str(train_it) + ", PC" + str(p))
            y_train_pred_prob = best_train_model.predict_proba(X_train_corrected)
            results_dict["PC" + str(p)]['train_best_params'][train_it] = best_params
            already_trained_auc = roc_auc_score(y_true = y_train, y_score = y_train_pred_prob[:,1])
            results_dict["PC" + str(p)]['train_auc_trained'].append(already_trained_auc)
            print("already train Rf " + str(already_trained_auc))
            
            newly_trained_auc = RF_cv(X_train_corrected,y_train,best_params)
            
            results_dict["PC" + str(p)]['mean_train_cv_auc'].append(np.mean(newly_trained_auc))
            print("mean train Rf " + str(np.mean(newly_trained_auc)))
            
            
            
            # perform cv result on test with best param
            test_RF = RF_cv(X_test_corrected,y_test,best_params)
            results_dict["PC" + str(p)]['mean_test_cv_auc'].append(np.mean(test_RF))
            print("mean test Rf " + str(np.mean(test_RF)))
            
            y_test_pred_prob = best_train_model.predict_proba(X_test_corrected)
            already_trained_test_auc = roc_auc_score(y_true = y_test, y_score = y_test_pred_prob[:,1])
            
            print("already train model test RF" + str(already_trained_test_auc))
            results_dict["PC" + str(p)]['test_auc_trained'].append(already_trained_test_auc)

            pickle.dump(all_datasets_dict , open( data_folders[0] +"_" + special_name + "_MINERVA_tt_grid.pkl", "wb" ) )
             
            
            
        train_it += 1
        test_train_end = timer()
        print("Finished one test train split")
        print(test_train_end - test_train_start)

    all_datasets_dict["dataset" + str(d)] = results_dict
    
    dataset_end = timer()
    print("Finished one dataset")
    print(dataset_end - dataset_start)

# ...
end = timer()
print(end - start)
        
pickle.dump(all_datasets_dict , open( data_folders[0] +"_" + special_name + "_MINERVA_tt_grid.pkl", "wb" ) )
            
        
    


