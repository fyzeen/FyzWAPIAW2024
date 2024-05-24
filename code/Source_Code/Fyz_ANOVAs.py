import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as sci
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import MultiComparison
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
    
    aseg_df["IDP_type"] = "aseg"
    dkt_df["IDP_type"] = "dkt"

    return pd.concat([dkt_df, aseg_df])

def concat_datasets(dataset_list=["ABCD", "ANXPE", "hcp_aging", "hcp_dev", "hcp_young", "UKB"], path_to_reg_outputs="/scratch/l.lexi/WAPIAW2024/data_for_figures", confounds="allconf"):
    '''
    ** This code is only written to work on structural IDPs atm **
    Inputs 
    ----------
    dataset_list: list
        List of datasets to concatenate into one large dataframe:
        [ABCD, ANXPE, hcp_aging, hcp_dev, hcp_young, UKB]
    path_to_reg_outputs: str or path
        Path to all regresion outputs with a standardized naming convention
    confounds: str
        if "both", will concatenate tables for both confound and no-confound refression
        if "allconf", will concatenate tables for ONLY regressions using confounds
        if "noconf", will concatenate tables for ONLY regressions using NO confounds
    '''
    for i, dataset in enumerate(dataset_list):
        if confounds == "both":
            dkt_fpath_allconf = os.path.join(path_to_reg_outputs, f"{dataset}_dkt_allconf.csv")
            aseg_fpath_allconf = os.path.join(path_to_reg_outputs, f"{dataset}_aseg_allconf.csv")
            concatenated_allconf = concat_aseg_dkt(dkt_fpath_allconf, aseg_fpath_allconf)
            concatenated_allconf["confounds"] = "allconf"

            dkt_fpath_noconf = os.path.join(path_to_reg_outputs, f"{dataset}_dkt_noconf.csv")
            aseg_fpath_noconf = os.path.join(path_to_reg_outputs, f"{dataset}_aseg_noconf.csv")
            concatenated_noconf = concat_aseg_dkt(dkt_fpath_noconf, aseg_fpath_noconf)
            concatenated_noconf["confounds"] = "noconf"

            concatenated = pd.concat([concatenated_allconf, concatenated_noconf]) 

        else:
            dkt_fpath = os.path.join(path_to_reg_outputs, f"{dataset}_dkt_{confounds}.csv")
            aseg_fpath = os.path.join(path_to_reg_outputs, f"{dataset}_aseg_{confounds}.csv")
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
        

def threeWayANOVA(df, type=2):
    '''
    This function will perform a 3-way anova (dependent variable: 'coef_', independent variables: ['Yeo_names', 'Dataset', and 'IDP_type'])

    Input 
    ----------
    df: pandas.DataFrame
        Concatenated df with all the betas for all the regressors as outputted by concat_datasets()
    
    type: int 
        Integer (either 1, 2, or 3) to define the "type" of ANOVA -- more information here: https://www.r-bloggers.com/2011/03/anova-%e2%80%93-type-iiiiii-ss-explained/ and https://www.statsmodels.org/dev/generated/statsmodels.stats.anova.anova_lm.html#statsmodels.stats.anova.anova_lm

    Output
    ----------
    anova_out_df: pandas.DataFrame
        A pandas dataframe with summary fo the ANOVA
    ----------

    '''
    df = df.loc[(df["regressor"] == "IDP (feature)") | (df["regressor"] == "IDP (dkt)") | (df["regressor"] == "IDP (aseg)")]
    #perform three-way ANOVA
    model = ols("""coef_ ~ C(Dataset) + C(Yeo_name) +
                C(Dataset):C(IDP_type) + C(Dataset):C(Yeo_name) + C(IDP_type):C(Yeo_name) +
                C(Dataset):C(IDP_type):C(Yeo_name)""", data=df).fit()
    
    return sm.stats.anova_lm(model, typ=type)








#could make two way anova work on giant dataset- it would use the dataset column to filter out extraneous rows!!
#def threeWayANOVA(data,outcome,independent_variables):
    #purpose- runs a three way anova with networks, dataset, and IDP type as factors
    #Inputs
    #data- csv containing all data from all studies, concatenated vertically
    #outcome- the name of the column containing the networks
    #independent_variables-
#    return


df = concat_datasets(path_to_reg_outputs="/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/data_for_figures")

df = df.loc[(df["regressor"] == "IDP (feature)") | (df["regressor"] == "IDP (dkt)") | (df["regressor"] == "IDP (aseg)")]

#perform three-way ANOVA
model = ols("""coef_ ~ C(Dataset) + C(Yeo_name) + C(Dataset):C(Yeo_name)""", data=df).fit()

print("\n\n Two-Way ANOVA with Yeo_name and Dataset (no absolute value) \n\n")
print((sm.stats.anova_lm(model, typ=2)))

# post hoc for Yeo_name
mc = MultiComparison(df['coef_'], df['Yeo_name'])
posthoc_results = mc.tukeyhsd()
print("\n\n Post Hoc for Yeo_name no absolute value \n\n")
print(posthoc_results)

# post_hoc for Dataset
mc = MultiComparison(df['coef_'], df['Dataset'])
posthoc_results = mc.tukeyhsd()
print("\n\n Post Hoc for Dataset (no absolute value) \n\n")
print(posthoc_results)

########### absolute value analyses ############
abs_df = df.copy()
abs_df["coef_"] = abs_df["coef_"].abs()

#perform three-way ANOVA
model = ols("""coef_ ~ C(Dataset) + C(Yeo_name) + C(Dataset):C(Yeo_name)""", data=abs_df).fit()

print("\n\n Two-Way ANOVA with Yeo_name and Dataset - absolute value \n\n")
print((sm.stats.anova_lm(model, typ=2)))

# post hoc for Yeo_name 
mc = MultiComparison(abs_df['coef_'], abs_df['Yeo_name'])
posthoc_results = mc.tukeyhsd()
print("\n\n Post Hoc for Yeo_name - absolute value \n\n")
print(posthoc_results)

# post_hoc for Dataset
mc = MultiComparison(abs_df['coef_'], abs_df['Dataset'])
posthoc_results = mc.tukeyhsd()
print("\n\n Post Hoc for Dataset - absolute value \n\n")
print(posthoc_results)


####### WITH CONFOUNDS VS NO CONFOUNDS ######
confounds_df = concat_datasets(path_to_reg_outputs="/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/data_for_figures", confounds="both")
confounds_df = confounds_df.loc[(confounds_df["regressor"] == "IDP (feature)") | (df["regressor"] == "IDP (dkt)") | (df["regressor"] == "IDP (aseg)")]

#perform two-way ANOVA
model = ols("""coef_ ~ C(Dataset) + C(confounds) + C(Dataset):C(confounds)""", data=confounds_df).fit()

print("\n\n Two-Way ANOVA with Dataset and confounds (no absolute value) \n\n")
print((sm.stats.anova_lm(model, typ=2)))

# perform two-way anova
model = ols("""coef_ ~ C(Yeo_name) + C(confounds) + C(Yeo_name):C(confounds)""", data=confounds_df).fit()

print("\n\n Two-Way ANOVA with Yeo_name and confounds (no absolute value) \n\n")
print((sm.stats.anova_lm(model, typ=2)))

# perform three-way ANOVA
model = ols("""coef_ ~ C(Dataset) + C(Yeo_name) + C(confounds) +
            C(Dataset):C(Yeo_name) + C(Dataset):C(confounds) + C(Yeo_name):C(confounds)""", data=df).fit()

print("\n\n Three-Way ANOVA with Yeo_name Dataset and confounds (no absolute value) \n\n")
print((sm.stats.anova_lm(model, typ=2)))

# post_hoc for Dataset
mc = MultiComparison(confounds_df['coef_'], confounds_df['confounds'])
posthoc_results = mc.tukeyhsd()
print("\n\n Post Hoc for confounds (no absolute val) \n\n")
print(posthoc_results)



#tukey_results = pairwise_tukeyhsd()

#impact of depression phenotypes
#each dataset also needs anova run on it alone
#1st factor is phenotypes- get from unique values? or hard code per dataset
#second factor is IDP-type, which Ty is adding to the data now

#impact of study sample

#impact of confounds
#3 way. factors are confounds, dataset, and network
#what does collapse over phenotypes mean




# MultiComparison is can be imported using `from statsmodels.stats.multicomp import MultiComparison`
def post_hoc(df, factor_var, response_var="coef_"):
    mc = MultiComparison(df[response_var], df[factor_var])
    posthoc_results = mc.tukeyhsd()
    return posthoc_results