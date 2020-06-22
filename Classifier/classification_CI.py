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
from sklearn.preprocessing import StandardScaler
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

import sys
from collections import Counter
occurrences = lambda s, lst: (i for i,e in enumerate(lst) if e == s)
args = sys.argv
print(args)
print(len(args))
greater_folder = args[1]
study_name = args[2]
prefix_name = args[3]
column_of_interest = args[4]
methods = args[5].split("&")
n_repeats = int(args[6])
data_type = args[7]
batch_def_folder = args[8]
if len(args) > 9:
    label_pos_or_neg = int(args[9])
    target_label = args[10]
    print(target_label)
else:
    label_pos_or_neg = 1
    target_label = 1
data_folder = greater_folder + "/data/" + study_name + "/"   
plot_folder = greater_folder + "/plots/" + study_name + "/" #+ 
methods_dict = utils.load_data(data_folder,prefix_name,methods,batch_column = batch_def_folder)
metadata = pd.read_csv(data_folder + "metadata.txt",delimiter="\t")
if column_of_interest == "antibiotic" and "AGP" in study_name:
    bin_antibiotic = utils.binarize_labels_mod(metadata["antibiotic_history"],pos_labels =['Year','Month','6 months','Week'],none_labels = ["Not provided",float("Nan"),'not provided'])
    #Counter(metadata["antibiotic"])
    metadata[column_of_interest] = bin_antibiotic
    pos_label = 1 #"Healthy"#'1-2' #'0-0.5'#'Omnivore' # '0-1.5'
# if column_of_interest == "antibiotic" and "Hispanic" in study_name:
#     bin_antibiotic = utils.binarize_labels_mod(metadata["antibiotic"],pos_labels =['1'],none_labels = ["Not provided",float("Nan"),'not provided'])
#     #Counter(metadata["antibiotic"])
#     metadata[column_of_interest] = bin_antibiotic
#     pos_label = 1 #"Healthy"#'1-2' #'0-0.5'#'Omnivore' # '0-1.5'
elif column_of_interest == "age_of_reloc":
    bin_column_of_interest = utils.binarize_labels_mod(metadata["agegroup_c6_v2.x"],pos_labels =['1','2','3','4'],none_labels = ["not applicable",float("Nan"),'not provided'])
    metadata[column_of_interest] = bin_column_of_interest
    column_of_interest = column_of_interest
    pos_label = 1 #"Healthy"#'1-2' #'0-0.5'#'Omnivore' # '0-1.5'
elif "trimester" in column_of_interest:
    #"1st trimester"
    if "trimester1" in column_of_interest:
        trimester_selected = "1st trimester"
    elif "trimester2" in column_of_interest: 
        trimester_selected = "2nd trimester"
    elif "trimester3" in column_of_interest:
        trimester_selected = "3rd trimester"
    eligible_rows = metadata['trimester'][metadata['trimester'] == trimester_selected].index.values
    metadata = metadata.loc[eligible_rows]
    bin_column_of_interest = utils.binarize_labels_mod(metadata["preg_outcome"],none_labels = ["not applicable",float("Nan"),'not provided'],pos_labels =[target_label])
    metadata[column_of_interest] = bin_column_of_interest
    column_of_interest = column_of_interest
    pos_label = 1   
    for method in methods:
        eligible_columns_temp = [col for col in eligible_rows if col in methods_dict[method].columns]
        methods_dict[method] = methods_dict[method][eligible_columns_temp]  
elif column_of_interest == "preg_outcome_sub":
    bin_column_of_interest = utils.binarize_labels_mod(metadata["preg_outcome"],none_labels = ["not applicable",float("Nan"),'not provided'],pos_labels =[target_label])
    metadata[column_of_interest] = bin_column_of_interest
    pos_label = 1 #"Healthy"#'1-2' #'0-0.5'#'Omnivore' # '0-1.5
else:
    if len(args) > 9:
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
print(Counter(metadata[column_of_interest] ))
names = ["Random Forest","Naive Bayes"]
all_methods_metrics = dict()
all_methods_means = pd.DataFrame(index = methods, columns= names)
all_methods_auc_stats = dict()
n_splits = 3
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
    np_interest =  np.array(metadata[column_of_interest])
    unique_categories = np.unique(np_interest)
    for i in unique_categories:
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
        if label_pos_or_neg == 3:
            category_counter = dict(Counter(metadata[column_of_interest]  ))
            categorical_counts = [category_counter[key] for key in category_counter.keys()]
            bootstrap_sample_size  = np.min(categorical_counts)
            print("bootstrap size" + str(bootstrap_sample_size))
        else:
            bootstrap_sample_size = int(bootstrap_prop * len(eligible_columns))
        sampled_columns += random.sample(list(eligible_columns), bootstrap_sample_size)
        print(len(sampled_columns))
    X = np.array(methods_dict[method][sampled_columns].transpose())
    y = np.array(metadata.loc[sampled_columns][column_of_interest]) 
    na_mask = pd.isna(y)
    y = y[~na_mask]
    X = X[~na_mask,:]
    print("Shapes")
    print(y.shape)
    print(X.shape)
    
    rskf = model_selection.RepeatedStratifiedKFold(n_splits=5, n_repeats=n_repeats, random_state=123)
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
        RandomForestClassifier(max_depth=5, random_state=0),
        GaussianNB()]#,AdaBoostClassifier(),SVC(kernel="linear", C=0.025),SVC(gamma=2, C=1),KNeighborsClassifier(3)]
#         MLPClassifier(alpha=1, max_iter=1000),
#         AdaBoostClassifier()]
    # iterate over classifiers
    # RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    classifier_it = 0
    for name, clf in zip(names, classifiers):
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
        for train_index, test_index in rskf.split(X, y):
            #print(cv_it)
            X_train, X_test = X[train_index,], X[test_index,]
            y_train, y_test = y[train_index], y[test_index]
            clf.fit(X_train, y_train)
            score = clf.score(X_test, y_test)
            #print(score)
            metric_classifier.loc[cv_it,names[classifier_it]] = score
            y_scores = clf.predict(X_test)
            if label_pos_or_neg == 3 and target_label == "1or0":
                y_scores = [0 if k[-1] == '0' else 1 for k in y_scores]
                y_test = [0 if k[-1] == '0' else 1 for k in y_test]
            
            
            accuracy_all.append(accuracy_score(y_test, y_scores))

            if label_pos_or_neg != 3 or target_label == "1or0":
            
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
            #all_methods_auc_stats[method][names[classifier_it]]['auc'] = auc_list #auc_matrix
            #all_methods_auc_stats[method][names[classifier_it]]['fpr'] = fpr_list #fpr_matrix
            #all_methods_auc_stats[method][names[classifier_it]]['tpr'] = tpr_list #tpr_matrix
        classifier_it += 1
    all_methods_metrics[method] = metric_classifier
    all_methods_means.loc[method,:] = np.array(pd.DataFrame.mean(metric_classifier,axis =0))   
    pickle.dump(all_methods_metrics , open( data_folder +  data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_classification_metrics.pkl", "wb" ) )
    pickle.dump(all_methods_auc_stats , open( data_folder + data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_classification_auc.pkl", "wb" ) )
    pickle.dump(all_methods_means , open( data_folder + data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_classification_means.pkl", "wb" ) )




