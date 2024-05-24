import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as sci
import visualization
import pingouin as pg
import warnings
from statsmodels.stats.multicomp import MultiComparison

#are all column names the same? We don't have information on what the individual idp types are atm
#need a column that is IDP type- i.e. glasser, schaefer, dkt. Ty is adding this, everyone will need to rerun their data
dataset_list = ['UKB', 'HCP_YA', 'HCP_Aging', 'HCP_Development', 'ABCD', 'ANXPE']

#omnibus needs to have a merged dataset. 3 vars- network (7 yeo + subcort
#one for datasets
#outcome var is coefficient

#new approach is a function that can handle an arbitrary number of factors, presented in list
#exaple_function(data,factor_list[dataset, network, etc])

def one_way_anovas(data,outcome, factors,outname):
    '''
    runs one way anovas
    :param data: the dataset with all your variables of interest as a pandas dataframe.
    :param outcome: the predicted variable; usually coefficient in this study! this is the column name
    :param factors: a LIST of column names that you want to be the independent variables.
    :return: the anova output in csv form
    '''
    #read factor list
    #loop and do one way anovas

    for factor in factors:
        factor_levels = data[factor].unique()
        factor_coef = []
        for level in factor_levels:

            add = data.loc[data[factor] == level, outcome]
            factor_coef.append(add)

        a = sci.f_oneway(*factor_coef)
        a_frame = pd.DataFrame(a)
        a_frame['component'] = ['F statistic','pvalue']
        a_frame.columns = ['value', 'component']
        a_frame.to_csv(f"{outname}_{factor}_anova_output.csv")
    return


def concat_aseg_dkt(dkt_fpath, aseg_fpath):
    aseg_df = pd.read_csv(aseg_fpath, header=0, index_col=0)
    dkt_df = pd.read_csv(dkt_fpath, header=0, index_col=0)

    return pd.concat([dkt_df, aseg_df])


def concat_datasets(dataset_list=["ABCD", "ANXPE", "hcp_aging", "hcp_dev", "hcp_young", "UKB"], 
                    path_to_reg_outputs="/scratch/l.lexi/WAPIAW2024/data_for_figures", 
                    confounds="allconf"):
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
    
    dictionary = pd.read_csv("/scratch/l.lexi/WAPIAW2024/Source_Code/maps/structural_dictionary.csv") 
    func_dictionary = pd.read_csv("/scratch/l.lexi/WAPIAW2024/Source_Code/maps/functional_dictionary.csv")

    # we reset the `out` index because we need the indices for each row to be unique for map_yeo() to work
    out = visualization.map_yeo(out.reset_index(), dictionary, func_dictionary)

    return out


def welchWrapper(outdir, dataDir='/scratch/l.lexi/WAPIAW2024/data_for_figures', 
                 group="Dataset", 
                 dv='coef_', 
                 between='phenotype',
                 dataset_list=["ABCD", "ANXPE", "hcp_aging", "hcp_dev", "hcp_young", "UKB"], 
                 confounds="allconf"):
    '''
    This function will run Welch ANOVA(s)
    '''
    #go through data directory and concatenate into one dataset
    concatData = concat_datasets(dataset_list=dataset_list, path_to_reg_outputs=dataDir, confounds=confounds)
    #concatData = concat
    #runs welch one way anovas
    output = pd.DataFrame()
    
    if group == 'False':
        #if group is False, run welch on entire dataset
        welch = pg.welch_anova(dv=dv, between=between, data=concatData)
        welch = pd.DataFrame(welch)
        welch['Group'] = 'All Subjects'
        output.to_csv(f"{between}_welch.csv")
    else:
        #identify our groups- these will be unique values in the group column
        group_list = concatData[group].unique()

        #loop through our list of groups. Do welch on each subset of the data
        for i, sub_group in enumerate(group_list):
            #filter to one group at a time
            temp = concatData.loc[concatData[group] == sub_group]

            welch = pg.welch_anova(dv = dv,between=between,data=temp)
            welch = pd.DataFrame(welch)
            welch['Group'] = sub_group

            mc = MultiComparison(temp[dv], temp[between])
            posthoc_results = mc.tukeyhsd()
            post_tukey = pd.DataFrame(data=posthoc_results._results_table.data[1:], columns=posthoc_results._results_table[0])
            

            post_tukey.to_csv(os.path.join(outdir, f"{sub_group}_{between}_tukey.csv"))

            if i == 0:
                #initialize dataframe
                output = welch
                #post = post_tukey
            else:
                #if we've been through the loop before, add current welch to output
                output = pd.concat([output,welch],axis=0)
                #post = pd.concat([post,post_tukey],axis=0)

        #save results
        out_path = os.path.join(outdir, f"{between}_welch.csv")
        output.to_csv(out_path)

    return output


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Run One-Way Welch ANOVA")

    parser.add_argument(
        '-o',
        '--out_dir',
        type=str,
        help='Path to directory in which you want to output the Welch ANOVA output.')
    
    parser.add_argument(
        '-d',
        '--data_dir',
        type=str,
        default="/scratch/l.lexi/WAPIAW2024/data_for_figures",
        help="Path to directory containing all regression output files."
    )

    parser.add_argument(
        "-g",
        "--group",
        default="False",
        type=str,
        help="Set to the column name of the variable you'd like to group by when doing Welch ANOVAs. If you don't want to group by any variables, either (1) don't input anything for this argument OR (2) input the STRING `False`"
    )

    parser.add_argument(
        "-y",
        "--dependent_var",
        default="coef_",
        type=str,
        help="Set to the name of the column of the variable you'd like to be the dependent variable in the ANOVA."
    )

    parser.add_argument(
        "-b",
        "--between",
        default="phenotype",
        type=str,
        help="Set to name of the column of the variable you'd like to be the independent variable in the ANOVA"
    )

    parser.add_argument(
        "-l",
        "--dataset_list",
        default="ABCD,ANXPE,hcp_aging,hcp_dev,hcp_young,UKB",
        type=str,
        help="Comma-separated list of the datasets you'd like to include in the ANOVA. Must look like: `ABCD,ANXPE,hcp_aging,hcp_dev,hcp_young,UKB` **No spaces between!**"
    )

    parser.add_argument(
        "-c",
        "--confounds",
        default="allconf",
        type=str,
        help="if `both`, will concatenate tables for both confound and no-confound refression; if `allconf`, will concatenate tables for ONLY regressions using confounds; if `noconf`, will concatenate tables for ONLY regressions using NO confounds"
    )

    args = parser.parse_args()

    dataset_list = [item for item in args.dataset_list.split(",")]

    welchWrapper(outdir=args.out_dir, 
                 dataDir=args.data_dir,
                 group=args.group,
                 dv=args.dependent_var,
                 between=args.between,
                 dataset_list=dataset_list,
                 confounds=args.confounds)



    
