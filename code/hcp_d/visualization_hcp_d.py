import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import warnings

from concat_aseg_dkt import concat_aseg_dkt

warnings.filterwarnings('ignore')


def plopiwaw_swarm(subgroup, hue_column):
    subgroup = subgroup[(subgroup['regressor'] == 'IDP (feature)') | (subgroup['regressor'] == 'IDP (dkt)') | (subgroup['regressor'] == 'IDP (aseg)')]
    sns.swarmplot(data=subgroup,
                  x='Yeo_name',
                  y="coef_",
                  size=15,
                  hue=hue_column,  # trick swarm to plot different colors for groups
                  legend=False,  # turn this back on if you want to get rid of xlabel altogether
                  ax=ax)
    plt.xticks(rotation=45, fontsize=50, ha='right')
    plt.yticks(fontsize=50)
    plt.xlabel("Yeo network", fontsize=50)
    plt.ylabel("Coefficient", fontsize=50)
    plt.ylim(-0.2, 0.2)
    # plt.ylim(-1, 1)
    plt.tight_layout()  # neccessary to get the x-axis labels to fit
    plt.axhline(y=0, color='lightgrey', linestyle='-', lw=8)


def map_yeo(results, dictionary):
    dataset_list = ['UKB', 'HCYA', 'HCA', 'HCD', 'ABCD', 'ANXPE']
    mappings = {dataset: pd.Series(dictionary['Yeo'].values, index=dictionary[dataset]).to_dict() for dataset in
                dataset_list}
    for dataset in dataset_list:
        results.loc[results['Dataset'] == "HCP_Development", 'Yeo'] = results['IDP_name'].map(mappings[dataset])
    yeo_mapping = {
        1: "Visual",
        2: "Somatomotor",
        3: "Dorsal Attention",
        4: "Ventral Attention",
        5: "Limbic",
        6: "Frontoparietal",
        7: "Default",
        8: "Subcortical"
    }
    results['Yeo_name'] = results['Yeo'].map(yeo_mapping)
    results = results.sort_values(by='Yeo')
    return results


# load dictionary 
dictionary = pd.read_csv("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/data/structural_dictionary.csv") # originally from "/scratch/l.lexi/WAPIAW2024/Source_Code/maps/structural_dictionary.csv" 

aseg_out_path = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/HCP_D/siteconf_aseg/linreg_siteconf_aseg.csv"
dkt_out_path = "/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/HCP_D/siteconf_dkt/linreg_siteconf_dkt.csv"

df = concat_aseg_dkt(aseg_out_path, dkt_out_path)


#df = pd.read_csv(regression_out_fpath, header=0, index_col=0)

# map yeo network
df = map_yeo(df, dictionary)


#
df['significance'] = df['p_fdr'].apply(lambda x: 'significant' if x <= 0.05 else 'not significant')
fig, ax = plt.subplots(figsize=(40, 30))
plopiwaw_swarm(df, 'significance')
#plt.show()
plt.savefig("/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/HCP_D/")

# To run the code in terminal:
# module load python
# source /export/anaconda/anaconda3/anaconda3-2020.07/bin/activate your_conda_env 
# python path_to_script
