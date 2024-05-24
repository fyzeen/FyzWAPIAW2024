import os
import numpy as np
import pandas as pd


## SAVE OUTPUT
##################################################################################################################################
def _get_regression_dicts(regmodel, df_cols, metadata, param_type, debug=False): 
    # initialize dict
    reg_dict = dict.fromkeys(df_cols)

    if debug:
        ### debugging code ###
        print("The regression-output formatting dictionary is initialized with keys:", reg_dict)
        ### debugging code ###

    params = getattr(regmodel, param_type)

    assert params.shape[0] == len(metadata["phenotypes"]), "number of fit phenotypes do not match number of regression coefficient rows"
    assert params.shape[1] == (len(metadata["confounds"])+1), "number of confounds are inconsistent with number of regression coefficient columns"
    
    reg_dictlist = []
    # assign constants; these will stay static but will be consistently re-appended
    reg_dict["IDP_name"] = regmodel.feature_names_in_[0]
    reg_dict["Dataset"] = metadata["dataset"]
    reg_dict["confounds"] = metadata["confounds"]
    reg_dict[param_type] = param_type
    reg_dict["score (R2)"] = regmodel.R2_score

    for x_idx in range(params.shape[1]):
        # make copy to avoid overwriting
        reg_dict = dict(reg_dict)
        # assign feature-dependent values
        feat_name = regmodel.feature_names_in_[x_idx]
        if regmodel.feature_names_in_[x_idx] in metadata["confounds"]:
            reg_dict["regressor"] = regmodel.feature_names_in_[x_idx]
        else:
            reg_dict["regressor"] = "IDP (feature)"
            
        for y_idx in range(params.shape[0]):
            # make copy to avoid overwriting
            reg_dict = dict(reg_dict)
            # assign phenotype-dependent values
            phenotype = metadata["phenotypes"][y_idx]
            #### we need to figure out if we wish to regress each phenotype separately!
            if debug:
                ### debugging code ###
                print(f"appending fit data for feature \'{feat_name}\' to phenotype \'{phenotype}\',")
                print(f"this corresponds to index pair [y,x]=[{y_idx},{x_idx}]")
                ### debugging code ###

            reg_dict["phenotype"] = phenotype

            # assign feature-phenotype-dependent values
            reg_dict[param_type] = params[y_idx, x_idx]
            reg_dict["p_raw"] = regmodel.p_raw[y_idx, x_idx]
            
            # append new dictionary to list of dictionaries
            reg_dictlist.append(reg_dict)
            if debug:
                ### debugging code ###
                print(f"New dictionary to append: {reg_dict}")
                ### debugging code ###

    if debug:
        ### debugging code ###
        print("List of dictionaries collecting data for single regression model:", reg_dictlist)
        ### debugging code ###
    return reg_dictlist

def _correct_pvals(fit_params, param_type, fit_intercept=False):
    return corr_pvals

def _get_metadata(X_fpath, Y_fpath, cnfvars_fpath=None):
    dataset_name = os.path.dirname(X_fpath.split("Data/")[1])
    phenotype_list = pd.read_csv(Y_fpath, index_col = 0, nrows=0).columns.values.tolist()
    conf_list = pd.read_csv(cnfvars_fpath, index_col = 0, nrows=0).columns.values.tolist()
    metadata = {
            "dataset": dataset_name,
            "phenotypes": phenotype_list,
            "confounds": conf_list
            }
    return metadata
##################################################################################################################################
