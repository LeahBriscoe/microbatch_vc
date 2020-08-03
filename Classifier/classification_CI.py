
# example run of this script


 # Goal: Classify samples as NORMAL(or adenoma) vs CRC


# python ./classifier.py /u/home/b/briscoel/gapped_kmer YuCRC patternABC kmer kmer_matrix bin_crc_adenomaORnormal patternABC 1 1 5 100 "entropy" 5 0.3 1 CRC
# the matrix you provide should have k-mer in row, sample in column

# your directory should be set up like this:
# greater_folder (e.g. /u/home/b/briscoel/gapped_kmer)
# > "data" (e.g. /u/home/b/briscoel/gapped_kmer/data)
# >> study_name  (e.g. /u/home/b/briscoel/gapped_kmer/data/YuCRC)
# >>> matrix_folder (e.g. /u/home/b/briscoel/gapped_kmer/data/YuCRC/patternABC)
# >>>> "kmer_matrix_patternABC.txt"  (e.g. /u/home/b/briscoel/gapped_kmer/data/YuCRC/patternABC/kmer_matrix_patternABC.txt)
# >>>> "metadata.txt" (e.g. /u/home/b/briscoel/gapped_kmer/data/YuCRC/patternABC/metadata.txt)

# greater_folder = /u/home/b/briscoel/gapped_kmer
# study_name = YuCRC
# matrix_folder = patternABC
# data_type = kmer
# prefix_name = kmer_matrix
# column_of_interest = bin_crc_adenomaORnormal
# methods = patternABC

# norm_input = 1
# map_with_accession = 1
# n_repeats = 5

# random forest parameters
# number_estimators_rf = 100

# if len(args) > 11:
#     label_pos_or_neg = 1
#     target_label = CRC



from sklearn import model_selection 
import pandas as pd
import utils
import numpy as np
from scipy import interpolate
import os
import seaborn as sns
import random
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.model_selection import train_test_split,KFold
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn import metrics
from sklearn.metrics import mean_squared_error 
from sklearn.linear_model import LinearRegression
from sklearn.metrics import accuracy_score
import math

import sys
from collections import Counter
occurrences = lambda s, lst: (i for i,e in enumerate(lst) if e == s)

args = sys.argv
print(args)
print(len(args))
greater_folder = args[1] # what folder do you save your different datasets in
study_name = args[2] # what is the name of the dataset (tells the program which folder to check)
matrix_folder = args[3] #which folder is the matrix saved in
data_type = args[4] # type of data. kmer vs OTU

prefix_name = args[5] # what is the prefix of the file name
column_of_interest = args[6] # what is the phenotype you are predicting (use the same name in the column of the metadata you want to predict), this programs reads from metadata.txt
methods = args[7].split("&") # you can specify prediction with multiple matrices by separating with &. Just ignore the '&'

norm_input = bool(int(args[8]))
map_with_accession = bool(int(args[9]))
n_repeats = int(args[10]) # number of different folds for cross validation. 5 or 10 is good. 


# RANDOM FOREST PARAMETERS
number_estimators_rf = int(args[11]) # number of estimators in random forest
criterion_input = args[12] 
min_samples_leaf_input = int(args[13])
max_features_input = float(args[14])
# END RANDOM FOREST PARAMETERS

if len(args) > 15:
    label_pos_or_neg = int(args[15]) # do you want to treat CRC as positive class or negative class? 
    target_label = args[16] # phenotype representing positive class or negative class? eg. CRC eg. H
    print(target_label)
else:
    label_pos_or_neg = 1
    target_label = 1
use_domain_pheno = False # for when running raw to compare to domain pheno
data_folder = greater_folder + "/data/" + study_name + "/"   
plot_folder = greater_folder + "/plots/" + study_name + "/" #+ 

#########################################################################
###### COMMENTARY: load data from your k-mer matrix, load metadata ######
#########################################################################

methods_dict = utils.load_data(data_folder,prefix_name,methods,batch_column = matrix_folder)
metadata = pd.read_csv(data_folder + "metadata.txt",delimiter="\t")



#########################################################################
###### COMMENTARY: get the column of interest: DiseaseStatus and convert it to binary labels 0 or 1 ######
#########################################################################

if map_with_accession:
    metadata.index = metadata['accession']

if column_of_interest == "antibiotic" and "AGP" in study_name:
    bin_antibiotic = utils.binarize_labels_mod(metadata["antibiotic_history"],pos_labels =['Year','Month','6 months','Week'],none_labels = ["Not provided",float("Nan"),'not provided'])
    #Counter(metadata["antibiotic"])
    metadata[column_of_interest] = bin_antibiotic
    pos_label = 1 #"Healthy"#'1-2' #'0-0.5'#'Omnivore' # '0-1.5'
else:
    if len(args) > 15:
        if label_pos_or_neg == 1:
            print("positive")
            bin_column_of_interest = utils.binarize_labels_mod(metadata[column_of_interest],none_labels = ["not applicable",float("Nan"),'not provided'],pos_labels =[target_label])
        elif label_pos_or_neg == 0:
            bin_column_of_interest = utils.binarize_labels_mod(metadata[column_of_interest],none_labels = ["not applicable",float("Nan"),'not provided'],neg_labels =[target_label])
        else: 
            bin_column_of_interest = utils.multiclassarize_labels_mod(metadata[column_of_interest],none_labels = ["not applicable",float("Nan"),'not provided'])
        
    else:

        bin_column_of_interest = utils.binarize_labels_mod(metadata[column_of_interest],none_labels = ["not applicable",float("Nan"),'not provided'])
    metadata[column_of_interest] = bin_column_of_interest
    pos_label = 1 #"Healthy"#'1-2' #'0-0.5'#'Omnivore' # '0-1.5


#########################################################################
###### COMMENTARY: Make the two classes even sized ######
#########################################################################

print(Counter(metadata[column_of_interest] ))
names = ["Random Forest","Naive Bayes"]
all_methods_metrics = dict()
all_methods_means = pd.DataFrame(index = methods, columns= names)
all_methods_auc_stats = dict()

bootstrap_prop = 0.80
for method in methods:
    random.seed(30)
    print(method)
    all_methods_auc_stats[method] = dict()
    #bootstrap_sample_size = int(bootstrap_prop * methods_dict[method].shape[1])
    #sampled_columns = random.sample(list(methods_dict[method].columns), bootstrap_sample_size)
    # class respective subsampling
    sampled_columns = []
    #for i in range(2):

    if use_domain_pheno:
        np_interest =  np.array(metadata['domain_pheno'] )
        unique_categories = np.unique(np_interest)
        
    else:
        np_interest =  np.array(metadata[column_of_interest])
        unique_categories = np.unique(np_interest)
        #print(unique_categories)
        #if len(args) <= 12:
        unique_categories = [cat for cat in unique_categories if not math.isnan(cat)]

    
    
    print(unique_categories)
    
    print("column of interest")
    #print(metadata[column_of_interest])   


    for i in unique_categories:
        #print(i)
        if use_domain_pheno:
            eligible_columns_all = metadata['domain_pheno'][metadata['domain_pheno'] == i].index.values
        else:
            eligible_columns_all = metadata[column_of_interest][metadata[column_of_interest] == i].index.values
        eligible_columns = [col for col in eligible_columns_all if col in methods_dict[method].columns]
        #print(eligible_columns)
        if "sub" in column_of_interest:
            X_host = eligible_columns
            y_host = list(metadata.loc[eligible_columns]["host_id"])
            new_eligibility = []
            for s in list(set(y_host)):
                new_eligibility.extend(random.sample(list(occurrences(s,y_host)),1))
            eligible_columns = list(np.array(X_host)[new_eligibility])
            bootstrap_prop= 1

        if use_domain_pheno:
            category_counter = dict(Counter(metadata['domain_pheno']  ))
            categorical_counts = [category_counter[key] for key in category_counter.keys()]
        else:
            category_counter = dict(Counter(metadata[column_of_interest]  ))
            categorical_counts = [category_counter[key] for key in category_counter.keys()]
            #print(category_counter)
            #if len(args) <= 12:
            categorical_counts = [category_counter[key] for key in category_counter.keys() if not math.isnan(key) ]
        #print(category_counter)


        #print(categorical_counts)
        bootstrap_sample_size  = np.min(categorical_counts)
        bootstrap_sample_size  = np.min([bootstrap_sample_size,len(eligible_columns)])
        #print("bootstrap_sample_size ")
        #print(bootstrap_sample_size )
        #print("length elidible_columns")
        #print(len(eligible_columns))
        #print("bootstrap size" + str(bootstrap_sample_size))

        # if label_pos_or_neg == 3:
        # else:
        #     bootstrap_sample_size = int(bootstrap_prop * len(eligible_columns))
        sampled_columns += random.sample(list(eligible_columns), bootstrap_sample_size)
        #print(len(sampled_columns))


    ###############################################################################
    ######  COMMENTARY: 
    ######  Make each sample sum to 1, so that regardless of number of reads per sample, they can be compared
    ######  Remove samples that have NA disease status
    ###############################################################################


    if norm_input:
        X  = normalize(np.array(methods_dict[method][sampled_columns].transpose()), axis = 1, norm = 'l1')
    else:
        X = np.array(methods_dict[method][sampled_columns].transpose())
    
    y = np.array(metadata.loc[sampled_columns][column_of_interest]) 

    y_counter = dict(Counter(y))
    print(y_counter)

    na_mask = pd.isna(y)
    y = y[~na_mask]
    X = X[~na_mask,:]
    print("Shapes")
    print(y.shape)
    print(X.shape)

    ###############################################################################
    ######  COMMENTARY: 
    ######  Calculate cross validation folds (train and test samples) with RepeatStratified KFold
    ######  Run the classifier
    ###############################################################################


    n_splits = 5
    rskf = model_selection.RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=123)
    h = .02  # step size in the mesh
    metric_classifier = pd.DataFrame(index = list(range(0,10)), columns= names)
#     classifiers = [
#         KNeighborsClassifier(3),
#         SVC(kernel="linear", C=0.025),
#         SVC(gamma=2, C=1),
#         GaussianProcessClassifier(1.0 * RBF(1.0)),
#         DecisionTreeClassifier(max_depth=5),
#         RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
#         MLPClassifier(alpha=1, max_iter=1000),
#         AdaBoostClassifier(),
#         GaussianNB(),
#         QuadraticDiscriminantAnalysis()]   
#     "Naive Bayes","AdaBoost","RBF SVM","Linear SVM"
    classifiers = [
        RandomForestClassifier(max_depth=5, random_state=0,n_estimators = number_estimators_rf,\
            criterion = criterion_input,min_samples_leaf = min_samples_leaf_input,max_features = max_features_input),
        GaussianNB()]#,AdaBoostClassifier(),SVC(kernel="linear", C=0.025),SVC(gamma=2, C=1),KNeighborsClassifier(3)]
#         MLPClassifier(alpha=1, max_iter=1000),
#         AdaBoostClassifier()]
    # iterate over classifiers
    # RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    classifier_it = 0
    for name, clf in zip(names, classifiers):
        if classifier_it:
            print("hyperparam")
            print(min_samples_leaf_input)
        print(names[classifier_it])
        all_methods_auc_stats[method][names[classifier_it]]  = dict()
        #print("TRAIN:", train_index, "TEST:", test_index)
        #auc_list = [] #pd.DataFrame(index = list(range(n_splits * n_repeats)))
        #tpr_list = [] #pd.DataFrame(index = list(range(n_splits * n_repeats)))
        #fpr_list = [] #pd.DataFrame(index = list(range(n_splits * n_repeats)))
        cv_it = 0
        auc_all = []
        accuracy_all = []
        y_tr_all = []
        y_pr_all = []
        y_tr = []
        y_pr = []
        #importances_matrix = pd.DataFrame(index=methods_dict[method].index
        importances_matrix = pd.DataFrame(index=methods_dict[method].index,columns = range(n_splits * n_repeats))


        training_iteration = 0
        for train_index, test_index in rskf.split(X, y):
            #print(cv_it)
            X_train, X_test = X[train_index,], X[test_index,]
            y_train, y_test = y[train_index], y[test_index]
            clf.fit(X_train, y_train)

            if classifier_it == 0: # if random forest
                print("rf here")
               
                importances = clf.feature_importances_
                print(type(importances))
                #print(list(importances)[0:10])
                #importances = pd.Series(importances)
                #importances.index = methods_dict[method].index
                #print("it" + str(training_iteration))
                importances_matrix[training_iteration] = list(importances)
                #pd.concat([importances_matrix,importances],axis=1)

            training_iteration += 1
            # importances = importances.sort_values(axis=0, ascending=False)
            # nonzero_count = sum(importances > 0)
            # print(importances[0:(nonzero_count + 1)])
            # print("rank of batch")
            # print(list(importances.index).index('batch_labels_factor2'))
            # print("rank of batch 3")
            # print(list(importances.index).index('batch_labels_factor3'))

            # for f in range(X.shape[1]):
            #     print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
            #np.savetxt(data_folder +  data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_rf_importance.csv", importances, delimiter=',')
                    

            score = clf.score(X_test, y_test)
            #print(score)
            metric_classifier.loc[cv_it,names[classifier_it]] = score
            y_scores = clf.predict(X_test)
            if label_pos_or_neg == 3 and target_label == "1or0":
                y_scores = [0 if k[-1] == '0' else 1 for k in y_scores]
                y_test = [0 if k[-1] == '0' else 1 for k in y_test]

                category_counter = dict(Counter(y_scores  ))
                print("y score counter")
                print(category_counter )

                print(y_scores[0:20])

                category_counter = dict(Counter(y_test  ))
                print("y test counter")
                print(category_counter )
                print(y_test[0:20])

            
            
            accuracy_all.append(accuracy_score(y_test, y_scores))

            if label_pos_or_neg != 3 or target_label == "1or0":
                print("get auc")
            
                y_scores_bin = clf.predict_proba(X_test).transpose()[1]
                y_test_bin = utils.binarize_labels(y_test,pos_label =pos_label)
                y_tr.append(y_test_bin)
                y_pr.append(y_scores_bin)
                y_tr_all.extend(y_test_bin)
                y_pr_all.extend(y_scores_bin)
                auc = metrics.roc_auc_score(y_test_bin, y_scores_bin)
                #fpr, tpr, thresholds = metrics.roc_curve(y_test_bin, y_scores_bin)
                #fpr_list.append(fpr)
                #tpr_list.append(tpr)
                auc_all.append(auc)
            #auc_matrix.loc[cv_it,:] = auc
            #fpr_matrix.loc[cv_it,:] = fpr
            #tpr_matrix.loc[cv_it,:] = tpr       
            cv_it += 1

        all_methods_auc_stats[method][names[classifier_it]]['accuracy'] = accuracy_all
        if label_pos_or_neg != 3 or target_label == "1or0":

            fpr_all, tpr_all, thresholds_all = metrics.roc_curve(y_tr_all, y_pr_all)
            tpr_i = []
            for s in range(20):
                fpr, tpr, thresholds = metrics.roc_curve(y_tr[s], y_pr[s])
                if not np.isnan(tpr[0]):
                    tpr_i.append(interpolate.interp1d(fpr, tpr, 'nearest')(fpr_all))       
            all_methods_auc_stats[method][names[classifier_it]]['fpr_all'] = fpr_all#auc_matrix
            all_methods_auc_stats[method][names[classifier_it]]['tpr_all'] = tpr_all#auc_matrix
            all_methods_auc_stats[method][names[classifier_it]]['tpr_i'] = tpr_i #auc_matrix
            all_methods_auc_stats[method][names[classifier_it]]['auc_all'] = auc_all #auc_matrix
            #all_methods_auc_stats[method][names[classifier_it]]['feature_importance'] = importances_matrix
            #all_methods_auc_stats[method][names[classifier_it]]['auc'] = auc_list #auc_matrix
            #all_methods_auc_stats[method][names[classifier_it]]['fpr'] = fpr_list #fpr_matrix
            #all_methods_auc_stats[method][names[classifier_it]]['tpr'] = tpr_list #tpr_matrix
        if classifier_it == 0: 
            importances_matrix.to_csv(data_folder + data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_RF_importance.csv") 

        classifier_it += 1
        
    all_methods_metrics[method] = metric_classifier
    all_methods_means.loc[method,:] = np.array(pd.DataFrame.mean(metric_classifier,axis =0))   
    pickle.dump(all_methods_metrics , open( data_folder +  data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_classification_metrics.pkl", "wb" ) )
    pickle.dump(all_methods_auc_stats , open( data_folder + data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_classification_auc.pkl", "wb" ) )
    pickle.dump(all_methods_means , open( data_folder + data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_classification_means.pkl", "wb" ) )
    


