import pandas as pd
import numpy as np
import scipy.stats

def load_data(data_folder,prefix_name, methods,batch_column):
    
    method_dict = dict()
    for m in methods:
        print(m)
        batch_corrected_matrix = pd.read_csv(str(data_folder + batch_column + "/" + prefix_name + "_" + m +".txt"),delimiter="\t")
        method_dict[m] = batch_corrected_matrix
    return method_dict

def load_pcscores(data_folder,method,data_type):

    pcscores = pd.read_csv(str(data_folder +"/" + "pcascore_" + data_type + "_" + method +".csv"))
    return pcscores

def binarize_labels(labels,pos_label):
    return [1 if lab==pos_label else 0 for lab in labels]


import numpy as np
import scipy.stats
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

def binarize_labels_mod(labels,none_labels, pos_labels = [],neg_labels = []):
    new_labels = []
    print("length pos labels")
    print(len(pos_labels))
    for lab in labels:
        if lab in none_labels:
            new_labels.append(None)
        elif len(pos_labels) > 0:
            if str(lab) in pos_labels:
                #print("truth")
                new_labels.append(1)
            else:
                new_labels.append(0)
        elif len(neg_labels) > 0:
            if str(lab) in neg_labels:
                new_labels.append(0)
            else:
                new_labels.append(1)
        else:
            new_labels.append(lab)

    return new_labels


def multiclassarize_labels_mod(labels,none_labels):
    new_labels = []
    print("length pos labels")
    print(len(pos_labels))
    for lab in labels:
        if lab in none_labels:
            new_labels.append(None)
        else:
            new_labels.append(lab)

    return new_labels



