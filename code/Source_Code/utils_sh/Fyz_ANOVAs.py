import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as sci
import warnings

import sys

sys.path.append("/scratch/l.lexi/WAPIAW2024/Source_Code")
# Fyzeen local
sys.path.append("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/code/Source_Code")

from visualization import map_yeo




#are all column names the same? We don't have information on what the individual idp types are atm
#need a column that is IDP type- i.e. glasser, schaefer, dkt. Ty is adding this, everyone will need to rerun their data
dataset_list = ['UKB', 'HCP_YA', 'HCP_Aging', 'HCP_Development', 'ABCD', 'ANXPE']

#omnibus needs to have a merged dataset. 3 vars- network (7 yeo + subcort
#one for datasets
#outcome var is coefficient

#new approach is a function that can handle an arbitrary number of factors, presented in list
#exaple_function(data,factor_list[dataset, network, etc])

def one_way_anovas(data,outcome, factors):
    '''

    :param data: the dataset, as a pandas dataframe
    :param outcome: the predicted variable; usually coefficient in this study! this is the column name
    :param factors: a LIST of column names that you want to be the independent variables.
    :return: the anova output in csv form
    '''
    #read factor list
    #loop and do one way anovas
    for factor in factors:
        factor_coef = []
        factor_levels = data[factor].unique()
        for level in factor_levels:

            add = data.loc[data[factor] == level, outcome]
            factor_coef.append(add)

        a = sci.f_oneway(*factor_coef)
        a.to_csv(f"{factor}_anova_output.csv")
        #do run way anova and save csv

    #repeat





##### FYZ AND SAMUEL CHANGED BELOW ######

def concat_aseg_dkt(dkt_fpath, aseg_fpath):
    aseg_df = pd.read_csv(aseg_fpath, header=0, index_col=0)
    dkt_df = pd.read_csv(dkt_fpath, header=0, index_col=0)
    dkt_df

    return pd.concat([dkt_df, aseg_df])

def concat_datasets(dataset_list=["ABCD", "ANXPE", "hcp_aging", "hcp_dev", "hcp_young", "UKB"], path_to_reg_outputs="/scratch/l.lexi/WAPIAW2024/data_for_figures"):
    '''
    ** This code is only written to work on structural IDPs atm **
    Inputs 
    ----------
    dataset_list: list
        List of datasets to concatenate into one large dataframe:
        [ABCD, ANXPE, hcp_aging, hcp_dev, hcp_young, UKB]
    '''
    for i, dataset in enumerate(dataset_list):
        dkt_fpath = os.path.join(path_to_reg_outputs, f"{dataset}_dkt_allconf.csv")
        aseg_fpath = os.path.join(path_to_reg_outputs, f"{dataset}_aseg_allconf.csv")
        concatenated = concat_aseg_dkt(dkt_fpath, aseg_fpath)

        if i==0:
            out = concatenated
        else:
            out = pd.concat([out, concatenated])
    
    #dictionary = pd.read_csv("/scratch/l.lexi/WAPIAW2024/Source_Code/maps/structural_dictionary.csv") 
    #func_dictionary = pd.read_csv("/scratch/l.lexi/WAPIAW2024/Source_Code/maps/functional_dictionary.csv")
    # Fyz local
    dictionary =  pd.read_csv("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/code/Source_Code/maps/structural_dictionary.csv") 
    func_dictionary = pd.read_csv("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/code/Source_Code/maps/functional_dictionary.csv")

    # we reset the `out` index because we need the indices for each row to be unique for map_yeo() to work
    out = map_yeo(out.reset_index(), dictionary, func_dictionary)

    return out
        

def threeWayANOVA(concatenated_df, network_col="Yeo_name"):
    return








#could make two way anova work on giant dataset- it would use the dataset column to filter out extraneous rows!!
#def threeWayANOVA(data,outcome,independent_variables):
    #purpose- runs a three way anova with networks, dataset, and IDP type as factors
    #Inputs
    #data- csv containing all data from all studies, concatenated vertically
    #outcome- the name of the column containing the networks
    #independent_variables-
#    return


df = concat_datasets(path_to_reg_outputs="/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/data_for_figures")
print(df)




#impact of depression phenotypes
#each dataset also needs anova run on it alone
#1st factor is phenotypes- get from unique values? or hard code per dataset
#second factor is IDP-type, which Ty is adding to the data now

#impact of study sample

#impact of confounds
#3 way. factors are confounds, dataset, and network
#what does collapse over phenotypes mean

