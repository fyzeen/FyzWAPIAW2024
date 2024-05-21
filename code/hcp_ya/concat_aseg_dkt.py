import pandas as pd

def concat_aseg_dkt(aseg_path, dkt_path):
    aseg_df = pd.read_csv(aseg_path, header=0, index_col=0)
    dkt_df = pd.read_csv(dkt_path, header=0, index_col=0)

    return pd.concat([dkt_df, aseg_df])


aseg_fpath = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/HCP_YA_asegYeo_allConfounds/HCP_YA_asegYeo_allConfounds.csv"
dkt_fpath = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/HCP_YA_dktYeo_allConfounds/HCP_YA_dktYeo_allConfounds.csv"

concatenated = concat_aseg_dkt(aseg_fpath, dkt_fpath)

concatenated.to_csv("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/HCP_YA_dktANDaseg_allConfounds/HCP_YA_dktANDaseg_allConfounds.csv")