import pandas as pd
import numpy as np
import scipy.stats

def load_data(data_folder,methods,data_type):
    
    method_dict = dict()
    for m in methods:
        batch_corrected_matrix = pd.read_csv(str(data_folder +"/" + "kmer_BatchCorrected_" + m +".txt"),delimiter="\t")
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