import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import warnings
import sys

sys.path.append("/scratch/l.lexi/WAPIAW2024/Source_Code")

warnings.filterwarnings('ignore')


def plopiwaw_swarm(subgroup, hue_column):
    subgroup = subgroup[(subgroup['regressor'] == 'IDP (feature)')]
    sns.swarmplot(data=subgroup,
                  x='Yeo_name',
                  y="coef_",
                  size=15,
                  hue=hue_column,
                  legend=True,
                  ax=ax,
                  dodge=True,
                  #hue_order=['']
                  )
    plt.xticks(rotation=45, fontsize=50, ha='right')
    plt.yticks(fontsize=50)
    plt.xlabel("Yeo network", fontsize=50)
    plt.ylabel("Coefficient", fontsize=50)
    plt.setp(ax.get_legend().get_texts(), fontsize='30')
    plt.setp(ax.get_legend().get_title(), fontsize='30')
    plt.ylim(-.4, .4)
    plt.tight_layout()  # neccessary to get the x-axis labels to fit
    plt.axhline(y=0, color='lightgrey', linestyle='-', lw=8)


def map_yeo(results, dictionary, func_dictionary):
    dataset_list = ['UKB', 'HCP_YA', 'HCP_Aging', 'HCP_Development', 'ABCD', 'ANXPE']
    mappings = {dataset: pd.Series(dictionary['Yeo'].values, index=dictionary[dataset]).to_dict() for dataset in
                dataset_list}
    for dataset in dataset_list:
        results.loc[results['Dataset'] == dataset, 'Yeo'] = results['IDP_name'].map(mappings[dataset])

    # map func_idps
    func_map = dict(zip(func_dictionary['node'].values, func_dictionary['Yeo'].values))
    results.loc[results['Yeo'].isna(), 'Yeo'] = results['IDP_name'].map(func_map)

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
    idp_name_mapping = dict(zip(dictionary['UKB'].values, dictionary['HCP_YA'].values))
    results.loc[results['Dataset'] == 'UKB', 'IDP_name'] = results['IDP_name'].map(idp_name_mapping)
    return results

def set_confounds_type(confounds):
    if pd.isna(confounds):
        return 'None'
    else:
        return 'all'

def set_confounds_type_detail(confounds):
    if pd.isna(confounds):
        return 'None'
    else:
        return confounds


def get_all_val(indir):
    vals = pd.DataFrame()
    file_list = [file for file in os.listdir(indir) if file.endswith('.csv')]
    for i in file_list:
        valset = pd.read_csv(os.path.join(indir, i), index_col=0, header=0)
        vals = pd.concat([vals, valset], axis=0, ignore_index=True)
    vals['confounds_type'] = vals['confounds'].apply(set_confounds_type)
    vals['confounds_type_details'] = vals['confounds'].apply(set_confounds_type_detail)
    return vals

'''
# load dictionary
dictionary = pd.read_csv("/scratch/l.lexi/WAPIAW2024/Source_Code/maps/structural_dictionary.csv") 
func_dictionary = pd.read_csv("/scratch/l.lexi/WAPIAW2024/Source_Code/maps/functional_dictionary.csv")
df = get_all_val('/scratch/l.lexi/WAPIAW2024/data_for_figures')
df = map_yeo(df, dictionary,func_dictionary)

#
df['significance'] = df['p_fdr'].apply(lambda x: 'significant' if x <= 0.05 else 'not significant')

fig, ax = plt.subplots(figsize=(60, 40))
plopiwaw_swarm(df, 'Dataset')
plt.title("comparing across datasets", fontsize=40)
plt.show()

fig, ax = plt.subplots(figsize=(60, 40))
plopiwaw_swarm(df, 'confounds_type')
plt.title("comparing across controlling confounds", fontsize=40)
plt.show()

data = 'UKB' # options: ['UKB', 'HCP_YA', 'HCP_Development', 'HCP_Aging' ,'ABCD', 'ANXPE']
subgroup = df[(df['Dataset'] == data)]
fig, ax = plt.subplots(figsize=(60, 40))
plopiwaw_swarm(subgroup, 'phenotype')
plt.title(f"comparing across phenotypes for {data}", fontsize=40)
plt.savefig(f"comparing across phenotypes for {data}")

# To run the code in terminal:
# module load python
# source /export/anaconda/anaconda3/anaconda3-2020.07/bin/activate your_conda_env 
# python path_to_script
'''