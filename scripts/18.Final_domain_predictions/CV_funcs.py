import pandas as pd
import numpy as np
import pickle

def add_domain_name_from_table_idx(table):
    
    table_domain_name = table.copy(deep=True)
    table_domain_name["domain_name"] = [x[:x.rfind("_")] for x in table_domain_name.index]
    
    return table_domain_name
#====================================================================================================================#

def calc_CV_idx_domains(X, splits_dict):
    """
    Calculate the train and test indices in each of the folds.
    Returns a dictionary with the test and train indices for each fold.
    """
    
    cv_idx = []
    #read the CV splits dictionary
 
        
    X_domain_name = add_domain_name_from_table_idx(X)
    all_idx = X_domain_name.index.tolist()
    
    for group_num in splits_dict.keys():
        
        #Adding the positions of all the group domains
        group_domains = splits_dict[group_num]["domains"]
        group_test_idx = []
        for domain in group_domains:
            domain_idx = X_domain_name[X_domain_name['domain_name'].str.match(domain)].index.tolist()
            group_test_idx.extend(domain_idx)
            
        #Adding the complementary indices as the training
        group_train_idx = [x for x in all_idx if x not in group_test_idx]
        
        cv_idx.append({"train": group_train_idx, "test": group_test_idx})
    
    return cv_idx
#====================================================================================================================#

def calc_CV_idx_iterative(X, splits_dict):
    """
    Calculate the train and test indices in each of the folds.
    Returns a dictionary with the test and train indices for each fold.
    """
    
    cv_idx = []
    #read the CV splits dictionary
 
    all_idx = X.index.tolist()
    
    for group_num in splits_dict.keys():
        
        #Defining the test idx as the intesrsection of group indices and X indices
        group_test_idx = list(set(all_idx) & set(splits_dict[group_num]["positions"]))
            
        #Adding the complementary indices as the training
        group_train_idx = [x for x in all_idx if x not in group_test_idx]
        
        cv_idx.append({"train": group_train_idx, "test": group_test_idx})
    
    return cv_idx
#====================================================================================================================#

def calc_CV_whole_domain_iterative(X, splits_dict):
    """
    Calculate the train and test indices in each of the folds for whole-domain indices.
    Returns a dictionary with the test and train indices for each fold.
    """
    
    cv_idx = []
    #read the CV splits dictionary
 
    all_idx = X.index.tolist()
    
    for group_num in splits_dict.keys():
        
        #Defining the test idx as the intesrsection of group indices and X indices
        group_test_idx = list(set(all_idx) & set(splits_dict[group_num]["domains"]))
            
        #Adding the complementary indices as the training
        group_train_idx = [x for x in all_idx if x not in group_test_idx]
        
        cv_idx.append({"train": group_train_idx, "test": group_test_idx})
    
    return cv_idx