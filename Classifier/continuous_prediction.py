from sklearn import model_selection 
import pandas as pd
import utils_bmi
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
from sklearn import linear_model
import sys

args = sys.argv
print(args)
print(len(args))
greater_folder = args[1]
study_name = args[2]
prefix_name = args[3]
column_of_interest = args[4]
methods = args[5].split("&")
n_cvs = int(args[6])
data_type = args[7]
batch_def_folder = args[8]
lin_model = args[9].strip()


data_folder = greater_folder + "/data/" + study_name + "/"
plot_folder = greater_folder + "/plots/" + study_name + "/" #+ 
methods_dict = utils_bmi.load_data(data_folder + batch_def_folder + "/",prefix_name,methods,data_type)

metadata = pd.read_csv(data_folder + "metadata.txt",delimiter="\t")
metadata["continuous_var"] = [float('Nan') if i == 'not applicable' or i == 'not provided' or 
                              i == "Not provided" or i== "Unspecified" or i == "Not applicable" 
                              else float(i) for i in list(metadata[column_of_interest])]

all_methods_stats = dict()



print("loaded and cleaned metadata")

bootstrap_prop = 0.80
# FIX THE SEED

for cv in range(n_cvs):
    
    for method in methods:
        
        random.seed(cv)
        if cv == 0:
            all_methods_stats[method] = dict()
            all_methods_stats[method]['pearson'] = []
            all_methods_stats[method]['mse'] = []
            
        
        print(method)
        X = []
        y = []


        
        bootstrap_sample_size = int(bootstrap_prop * methods_dict[method].shape[1])
        sampled_columns = random.sample(list(methods_dict[method].columns), bootstrap_sample_size)
        X = np.array(methods_dict[method][sampled_columns].transpose())
        y = np.array(metadata.loc[sampled_columns]["continuous_var"])
        y = [float("nan") if y_i == "Not provided" else float(y_i) for y_i in y]
        na_mask = pd.isna(y)
        y = np.array(y)
        y = y[~na_mask]
        X = X[~na_mask,:]

        #rskf = model_selection.RepeatedStratifiedKFold(n_splits= n_splits, n_repeats=n_repeats, random_state=123)
        rskf = model_selection.KFold(n_splits=5, random_state=123)
        
        cv_pearson = []
        cv_mse = []
        
        for train_index, test_index in rskf.split(X, y):
            #print(cv_it)
            X_train, X_test = X[train_index,], X[test_index,]
            y_train, y_test = y[train_index], y[test_index]

            if lin_model == "L1":
                reg = linear_model.Lasso(alpha=0.5).fit(X_train, y_train)
            else:

                reg = linear_model.LinearRegression().fit(X_train, y_train)
            pred = reg.predict(X_test)
            all_methods_stats[method]['pearson'].append(np.corrcoef(x=list(y_test),y=list(pred))[0,1])
            all_methods_stats[method]['mse'].append(mean_squared_error(list(y_test),list(pred) ))


        #all_methods_stats[method]['pearson'] = cv_pearson
        #all_methods_stats[method]['mse'] = cv_mse
        
            pickle.dump( all_methods_stats[method], open( data_folder + data_type + "_" + prefix_name + "_" + column_of_interest + "_" + method + "_lin_model_" + lin_model + "_pearson_and_mse.pkl", "wb" ) )
        





