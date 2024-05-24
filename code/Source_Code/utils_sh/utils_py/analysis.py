import pygam
import numpy as np
import statsmodels.api as sm
from scipy import stats
from sklearn import linear_model
from statsmodels.stats.diagnostic import het_breuschpagan


##################################################################################################################################
# convert fit type specification string into usable function -- move to outside module for bookkeeping?
def _get_regmodel(argument):
    _switcher = {
        "lstsq": linear_model.LinearRegression(fit_intercept=False),
        "LiGAM": pygam.LinearGAM(),
        "ridge": linear_model.Ridge(alpha=0.5),
        "lasso": linear_model.Lasso(alpha=0.1),
        "elnet": linear_model.ElasticNet(random_state=0),
        "bayes": linear_model.BayesianRidge()
     }
    _regmodel = _switcher.get(argument, lambda argument: print(f"Unknown regressor: {argument}. If you are using a custom function, please add it to the _switcher dictionary."))
    return _regmodel

# Breusch-Pagan Lagrange Multiplier test for heteroscedasticity:
# Tests the hypothesis that the residual variance does not depend on the variables in X
def _do_hetBP_test(regmodel, X, Y, debug=False):
    resid = Y - regmodel.predict(X)
    # Breusch-Pagan test implementation requires appendation of a constant column to (exogenous) data
    X_app = np.append(X, np.ones((X.shape[0],1)), 1)
    bp_test = het_breuschpagan(resid, X_app)
    bp_results = {
            'Lagrange Multiplier statistic': bp_test[0],
            'p-value': bp_test[1],
            'f-value': bp_test[2],
            'f_p-value': bp_test[3]
                  }

    setattr(regmodel, 'Breusch_Pagan_Test', bp_results)



# compute raw pvals -- assumes normal residuals and an intercept of 0 (i.e., data should be centered)
# if data is NOT centered and an intercept is fit, then X --> np.append(np.ones((len(X),1)), X, axis=1)
def _get_raw_pvals(regmodel, X, Y, fixed_effect=False, debug=False):
    Y_pred = regmodel.predict(X)

    # relies on assumption of no intercept; otherwise, we would divide by n+1 to calculate MSE
    deg_fdm = max(1,len(X) - len(X[0]))
    mse = np.power(Y - Y_pred,2).sum() / deg_fdm

    param_variance = mse * (np.linalg.pinv(np.dot(X.T,X)).diagonal())
    param_Tstats = regmodel.coef_ / np.sqrt(param_variance)
    pvals = np.asarray([2*(1 - stats.t.cdf(np.abs(T), deg_fdm)) for T in param_Tstats])

    if debug:
        ### debugging code ###
        p_nonzero = np.asarray([p for p in pvals.flatten() if p > 0])
        print("Raw pvals summary:")
        print(f"max = {max(pvals.flatten())}")
        print(f"number_nonzero = {len(p_nonzero)}")
        print(f"nonzero_min = {min(p_nonzero)}")
        print('')
        ### debugging code ###

    setattr(regmodel, 'p_raw', pvals)
    setattr(regmodel, 'model_degrees_of_freedom', deg_fdm)


# choose p-value correction method
def _correct_pvals(output_df, method="fdr"):
    _switcher = {
        "fdr": _fdr_pvals,
     }
    pval_corrector = _switcher.get(method, lambda argument: print(f"Unknown regressor: {argument}. If you are using a custom function, please add it to the _switcher dictionary."))
    output_df = pval_corrector(output_df)
    return output_df


# run family-wise error correction on raw p-values on formatted output from group regression
def _fdr_pvals(full_regression_df, debug=False, verbose=False):
    regressors = full_regression_df["regressor"].unique()
    phenotypes = full_regression_df["phenotype"].unique()
    for regtype in regressors:
        reg_idx = full_regression_df["regressor"] == regtype
        for phenotype in phenotypes:
            pheno_idx = full_regression_df["phenotype"] == phenotype
            reg_pheno_idx = reg_idx & pheno_idx
            reg_pheno_df = full_regression_df[reg_pheno_idx] 
            if debug:
                ### debugging code ###
                print(f"Pre-update values for first row sub-datafrate with regressor \'{regtype}\':")
                print(reg_pheno_df.iloc[1,:])
                ### debugging code ###

            # calculate false discovery rate control corrections to p-values
            p_fdr = stats.false_discovery_control(reg_pheno_df["p_raw"].values)           
            reg_pheno_df["p_fdr"] = p_fdr
            if debug:
                ### debugging code ###
                print(f"Post-update values for first row sub-datafrate with regressor \'{regtype}\':")
                print(reg_pheno_df.iloc[1,:])
                ### debugging code ###
            full_regression_df[reg_pheno_idx] = reg_pheno_df


    if debug:
        ### debugging code ###
        print("Post-FDR correction p-values across all regressions:")
        print(full_regression_df['p_fdr'])
        print(f"List of unique values: {list(full_regression_df['p_fdr'].unique())}")
        ### debugging code ###
    return full_regression_df


# formats regression output into parse-ready dataframe output
def _get_regression_dicts(regmodel, df_cols, metadata, param_type, debug=False): 
    # initialize dict
    reg_dict = dict.fromkeys(df_cols)

    params = getattr(regmodel, param_type)

    if debug:
        ### debugging code ###
        print("The regression-output formatting dictionary is initialized with keys:", reg_dict)
        print(f"For this (arbitrary) fit, parameters have shape {params.shape}")
        ### debugging code ###

    # assert params.shape[0] == len(metadata["phenotypes"]), "number of fit phenotypes do not match number of regression coefficient rows"
    if metadata["confounds"] is not None:
        assert params.shape[1] == (len(metadata["confounds"])+1), "number of confounds are inconsistent with number of regression coefficient columns"
        
    datatype = metadata["data type"]
    reg_dictlist = []
    # assign constants; these will stay static but will be consistently re-appended
    reg_dict["IDP_name"] = regmodel.feature_names_in_[0]
    reg_dict["Dataset"] = metadata["dataset"]
    reg_dict["confounds"] = metadata["confounds"]
    reg_dict[param_type] = param_type
    reg_dict["score (R2)"] = regmodel.R2_score
    #try:
    #    print(regmodel.R2_score)
    #except: 
        #pass
    reg_dict["data type"] = datatype
    if metadata["confounds"] is None:
        reg_dict["regressor"] = f"IDP ({datatype})"

    for x_idx in range(params.shape[1]):
        # make copy to avoid overwriting
        reg_dict = dict(reg_dict)
        # assign feature-dependent values
        feat_name = regmodel.feature_names_in_[x_idx]
        if metadata["confounds"] is not None:
            if regmodel.feature_names_in_[x_idx] in metadata["confounds"]:
                reg_dict["regressor"] = regmodel.feature_names_in_[x_idx]
            else:
                reg_dict["regressor"] = f"IDP ({datatype})"

        for y_idx in range(params.shape[0]):
            # make copy to avoid overwriting
            reg_dict = dict(reg_dict)
            # assign phenotype-dependent values
            phenotype = regmodel.target_name_[y_idx]
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
##################################################################################################################################
