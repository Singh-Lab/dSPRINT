#Imports
import pandas as pd
import numpy as np
import pickle
from os import getcwd

#Classifier imports
from sklearn.neighbors import KNeighborsClassifier 
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC


# Neural Net imports
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.optim.lr_scheduler import LambdaLR, ReduceLROnPlateau
from sklearn.model_selection import RepeatedStratifiedKFold

#General learning framework imports
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import auc, roc_auc_score, precision_recall_curve
from sklearn.utils import shuffle

from CV_funcs import calc_CV_idx_iterative, calc_CV_whole_domain_iterative
from NN_classes import Net, Net_tune, curr_device
#from generate_models_dict_global_auprc import Net

#CV splits dictionary
pfam_version = "31"
datafile_date = "08.06.18"
prec_th_str = "dna0.5_rna0.5_ion0.75"
folds_num = 5
curr_dir = getcwd()
exac_dir = curr_dir[:curr_dir.find("ExAC")]

models_req_scaling = ["SVM", "KNN", "Logistic", "NN"]
#====================================================================================================================#

def compute_per_domain_auc(y_test, pred_probs, pred_idx, classifier):
    """
    Compute the average per_domain auc and auprc for the test set
    """
    
    y_test_copy = y_test.copy(deep=True)
    y_test_copy["pred_probs"] = pred_probs
    
    domain_auc_list = []
    domain_auprc_list = []
    domain_auprc_ratio_list = []
    domain_name_list = []
    domain_pos_num_list = []
    domain_neg_num_list = []
    
    idx = y_test.index
    y_test_copy["domain_name"] = [x[:x.rfind("_")] for x in idx]
    domains_list = y_test_copy["domain_name"].unique().tolist()
        
    for domain_name in domains_list:
        
        #Get only the domain positions
        domain_df = y_test_copy[y_test_copy["domain_name"] == domain_name]

        #Find the binding and non-binding positions of this domain 
        bind_list = domain_df[domain_df["label"] == 1].index
        bind_idx = [int(x[len(domain_name)+1:]) for x in bind_list]
        bind_num = len(bind_idx)
        non_bind_list = domain_df[domain_df["label"] == 0].index
        non_bind_idx = [int(x[len(domain_name)+1:]) for x in non_bind_list]
        non_bind_num = len(non_bind_idx)
        if (bind_num == 0 or non_bind_num == 0):
            #No positions of one of the classes "binding/non-binding" - skipping"
            continue
            
        #Add number of positives and number of negatives
        domain_pos_num_list.append(bind_num)
        domain_neg_num_list.append(non_bind_num)
        #Compute domain AUC
        domain_auc = roc_auc_score(domain_df["label"], domain_df["pred_probs"])
        domain_auc_list.append(domain_auc)
        #Compute domain AUPRC
        precision, recall, thresholds = precision_recall_curve(domain_df["label"], domain_df["pred_probs"])
        domain_auprc = auc(recall, precision)
        domain_auprc_list.append(domain_auprc)
        #Add positives fraction to list
        pos_frac_ratio = bind_num/float(domain_df.shape[0])
        #Add ratio of AUPRC and positives fraction to list
        domain_auprc_ratio_list.append(domain_auprc/float(pos_frac_ratio))
        #Add domain name for AUC/AUPRC/Ratio tables
        domain_name_list.append(domain_name)
        
    #Compute the means for the lists 
    domain_auc_mean = np.mean(domain_auc_list)
    domain_auprc_mean = np.mean(domain_auprc_list)
    domain_auprc_ratio_mean = np.mean(domain_auprc_ratio_list)
    
    return (domain_auc_mean, domain_auprc_mean, domain_auprc_ratio_mean)
