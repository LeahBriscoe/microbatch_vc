#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


# In[2]:


data_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/'
plot_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/'



subfile = ""
select_columns_bool = False
shortened_shortened = False
trans_vec = False




phen_pretty = "colorectal cancer status" #"colorectal cancer status"
phen_type = "class"


folder_names = ["YuCRC","YuCRC"]
file_names = ['patternABC.txt","patterneven.txt"] 
select_columns_bool = True
select_labels = ["Even-Spaced","Odd-Spaced"] # how do you want the labels to look in the final plot 
metric_word = "auc_all" #what metric do you want to look at
chosen_classifier = ["Naive Bayes","Naive Bayes"] #["Naive Bayes" for i in range(len(file_names))] == ["Naive Bayes","Naive Bayes","Naive Bayes","Naive Bayes","Naive Bayes","Naive Bayes"]
title = 'AUC for prediction of ' + phen_pretty #'AUC for prediction of antibiotic history'#'Pearson correlation of predicted body mass index (BMI)'#'AUC for prediction of antibiotic history'#' #AUC for prediction of antibiotic history' #
pair_test_all = True
limit_spec =(0.5,1) # y-axis range
shortened=True
not_rotate=False
data_type = "kmer"

classifier = ["Naive Bayes","Random Forest"]#Random Forest"#"#'Regression'#"Random Forest"

#title = 'AUC for prediction of ' + phen_pretty #'AUC for prediction of antibiotic history'#'Pearson correlation of predicted body mass index (BMI)'#'AUC for prediction of antibiotic history'#' #AUC for prediction of antibiotic history' #
special_name = phen[0] #"Antibiotic history class prediction"#"BMI prediction" #"Antibiotic history class prediction"#"BMI prediction"# "Antibody prediction"# 


df_metric = pd.DataFrame()


# each method/file you are comparing will be a column in df_metric

for f in range(len(file_names)):
    
    
    filename_temp = data_folder + folder_names[f] + "/" + file_names[f] + "_classification_auc.pkl"
    print(filename_temp)
    if os.path.isfile(filename_temp):
        
        data_temp = pickle.load( open( filename_temp ,"rb"))
    
        #methods_lists.append([k for k in data1.keys()]) 
        
        if metric_word == "auc_all" or metric_word == "accuracy":
            index = folder_names[f] + file_names[f]
            #index = folder_names[f]
            #index = nice_names[f] + "_first_" + str(numPc[f])
            
            df_metric[index] = pd.Series(data_temp[file_names[f]][classifier[f]][metric_word])
            #nice_names[f] + "_first_" + str(numPc[f])
        else:
            #index = folder_names[f] + file_names[f]
            index = file_names[f]
            df_metric[index] = pd.Series(data_temp[metric_word])
    else: # if file not found, it will put 0 for that column
        print("not file")
        print(file_names[f])
        df_metric[nice_names[f] + file_names[f]] = 0
        df_metric[nice_names[f] + file_names[f]] = 0
        
    


select_columns = df_metric.columns #[ i + 'clr_pca_regress_out_no_scale_first10' for i in nice_names ] 

current_palette = sns.color_palette()
plot_color = current_palette.as_hex()

df_metric.mean(axis=0)


import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation

sns.set(style="whitegrid")


x = "Method"
y = "AUC"

ax = sns.boxplot(data=df_metric,palette=plot_color)
sns.set_context("paper", font_scale=1.5)   

add_stat_annotation(ax, data=df_metric,
    box_pairs=[(df_metric.columns[0], df_metric.columns[1]),
              (df_metric.columns[0], df_metric.columns[2])],
    test='t-test_ind', text_format='star', loc='outside', verbose=2)
if not_rotate:
    ax.set_xticklabels(labels = select_labels)
else:
    ax.set_xticklabels(rotation=90,labels = select_labels)
ax.set(ylim=limit_spec)
ax.set_title(title)
plt.savefig(plot_folder + chosen_classifier + data_type + "_" + trans[0] + '_boxplots_' + classifier_ofc + '_' + metric_word + "_" + special_name + '_2.pdf',bbox_inches='tight')



from scipy.stats import ttest_ind


ttest_ind(df_metric[df_metric.columns[0]], df_metric[df_metric.columns[1]])

