import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler


##################################################################################################################################
# subselects rows and columns the dataframe to retain only those specified, if any are specified.
# rows are specified by subject identifiers, presumed to do be strings
def _subset_data(df, sublist_fpath=None, colvars_fpath=None, debug=False):
    df.index = df.index.map(str)
    df.columns = df.columns.map(str)

    if sublist_fpath:
        with open(sublist_fpath,'r') as fin:
            sublist = fin.read().split()
            
            if debug:
                ### debugging code ###
                print(f"subject list: {sublist}")
                print(f"index list length: {len(df.index)}")
                print(f"intersection: {list(set(sublist) & set(df.index))}")
                print(f"subjects not in index list: {list(set(sublist) - (set(sublist) & set(df.index)))}")
                ### debugging code ###

        # keep only rows corresponding to list of chosen subject (if such a list is given)
        df=df.loc[sublist]

    if colvars_fpath:
        with open(colvars_fpath,'r') as fin:
            colvars = fin.read().split()

        # keep only columns corresponding to list of chosen variables (if such a list is given)
        df=df[colvars]

    return df

# X_fpath is presumed to be the immediate parent dictory of 
def _load_func_idp(X_fpath, sublist_fpath):
    if sublist_fpath:
        with open(sublist_fpath,'r') as fin:
            sublist = fin.read().split()
    else:
        sublist = None

    filelist = os.listdir(X_fpath)
    if sublist:
        sub_filelist = list(set(fname for fname in filelist for eid in sublist if fname.find(eid) != -1))
    else:
        sub_filelist = [fname for fname in filelist if 'sub-' in fname]
    sub_fpaths = [os.path.join(X_fpath, fname) for fname in sub_filelist]

    # pull imaging data
    X_data = np.array([_pull_subj_data(datapath) for datapath in sub_fpaths], dtype=float)
    colnames = [f"fNM_{i}" for i in range(X_data.shape[1])]

    # Format into simple dataframe w/ an EID column and IDP column(s?)
    X_df = pd.DataFrame(data=X_data, index=sublist, columns=colnames)
    X_df.index = X_df.index.map(str)
    X_df.columns = X_df.columns.map(str)

    return X_df


# normalizes (i.e., z-scores) given dataframe
def _normalize_df(df):
    norm_df = (df - df.mean())/df.std()
#    if df.shape[1] < 1:
#        df = df.to_frame()
#    scaler = StandardScaler()
#    df.iloc[:,0:-1] = scaler.fit_transform(df.iloc[:,0:-1].to_numpy())
    return norm_df


# loads data for a single subject given a subject ID number and general datapath string
def _pull_subj_data(datapath):
    data = np.genfromtxt(datapath)
    data = _validate_readin(data, datapath)
    return data

def _validate_readin(data, datapath, debug=False):
    if np.isnan(data).any():
        if datapath.endswith(".csv"):
            data = np.genfromtxt(datapath, delimiter=",")
            if np.isnan(data).any():
                raise Exception("At least one nan value in loaded data -- probably not a load-in error.")
        elif datapath.endswith(".txt"):
            raise Exception("At least one nan value in loaded data -- probably not a load-in error.")
        else:
            raise Exception(f"At least one nan value in loaded data -- did you mean to have a \".{datapath.split('.')[-1]}\" file extension?")


    if len(data.shape) > 1:
        assert len(data.shape) == 2, "classifier expects subject-wise input data to be either a matrix or vector."
        if data.shape == data.T.shape:
            if np.allclose(data, data.T,rtol=1e-4, atol=1e-6):  # changed tolerence
                data = _triu_vals(data)
                data = _handle_corrs(data)
            else:
                data = data.flatten()
        else:
            data = data.flatten()
    
    if debug:
        ### debugging code ###
        print(f"subect {os.path.basename(datapath)} has data of shape:", data.shape)
        ### debugging code ###

    return np.array(data.flatten(), dtype=float)

# returns (flattened) upper right triangle of a square matrix; discards matrix diagonal
def _triu_vals(A):
    n = A.shape[0]
    vals = A[np.triu_indices(n,1)]
    return vals.flatten()

def _handle_corrs(data):
    if (np.abs(data.flatten()) <= 1).all():
        new_data = np.arctanh(data)
    else:
        new_data = data
    return new_data

# sends information about learning problem parameters to output stream
def _show_input_params(X, Y, cnfvars=None):
    print('')
    print('*****Sample & Feature set sizes for variable set*****')
    print('Independent variable size:', X.shape)
    print('Dependent variable size:', Y.shape)
    if len(cnfvars):
        print('Covariates:', cnfvars.shape)
    print('')

# load namespace of parameter types based on inputs
def _get_metadata(X_fpath, Y_fpath, cnfvars_fpath=None, idplist_fpath=None, ptylist_fpath=None, cnflist_fpath=None):
    dataset_name = os.path.dirname(X_fpath.split("Data/")[1])

    if idplist_fpath:
        with open(idplist_fpath,'r') as fin:
            idp_list = fin.read().split()
    else:
        idp_list = pd.read_csv(X_fpath, index_col = 0, nrows=0).columns.values.tolist()

    if ptylist_fpath:
        with open(ptylist_fpath,'r') as fin:
            phenotype_list = fin.read().split()
    else:
        phenotype_list = pd.read_csv(Y_fpath, index_col = 0, nrows=0).columns.values.tolist()

    if cnfvars_fpath:
        if cnflist_fpath:
            with open(cnflist_fpath,'r') as fin:
                conf_list = fin.read().split()
        else:
            conf_list = pd.read_csv(cnfvars_fpath, index_col = 0, nrows=0).columns.values.tolist()


    metadata = {
            "dataset": dataset_name,
            "IDP list": idp_list,
            "phenotypes": phenotype_list,
            "confounds": conf_list
            }
    return metadata


## NEED TO FIX THIS TO WORK WITH str-VALUED INDEX LISTS!!!!
def _enforce_same_subjs(X, Y, conf=None):
    xy_index = set(X.index).intersection(Y.index)
    if np.any(conf):
        index = list(xy_index.intersection(set(conf.index)))
        conf = conf.loc[index]
    else:
        index = list(xy_index)

    X = X.loc[index]
    Y = Y.loc[index]
    return X, Y, conf
########################conf = conf.loc[list(index)]##########################################################################################################
