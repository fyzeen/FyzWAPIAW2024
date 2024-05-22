import numpy as np
import os
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

# Generate random p-values (replace this with your own p-values)
#np.random.seed(0)
#p_values = np.random.uniform(0, 1, size=1000)
def qq(phen='rds',study='UKB',showplot=True,indata='linreg_noconfs.csv',indir='/Users/petralenzini/work/WAPIAW/2024/QQ/',outdir='/Users/petralenzini/work/WAPIAW/2024/QQ/'):
    p_values=pd.read_csv(os.path.join(indir,indata), header=0)
    pvals=p_values.loc[(p_values.regressor.str.contains('IDP')) & (p_values.phenotype==phen)][['p_raw']]

    # Sort the p-values
    sorted_p_values = np.sort(pvals.p_raw)

    # Calculate the expected quantiles for a uniform distribution
    expected_quantiles = stats.uniform.ppf(np.linspace(0, 1, len(sorted_p_values)))

    # Create QQ plot
    plt.figure(figsize=(8, 6))
    plt.scatter(expected_quantiles, sorted_p_values, color='b', alpha=0.5)
    plt.plot([0, 1], [0, 1], color='r', linestyle='--')  # Diagonal line for reference
    plt.title('QQ Plot for Raw '+phen+' P-values in '+ study + ' '+indata)
    plt.xlabel('Expected Quantiles (Uniform Distribution)')
    plt.ylabel('Ordered P-values')
    plt.grid(True)
    plt.savefig(outdir + 'QQ Plot for Raw '+phen+' P-values in '+ study + '.png')
    if showplot:
        plt.show()

qq(phen='upps_neg',
    study='HCP-YA',
    showplot=False,
    indata='linreg_noconf_dkt.csv',
    indir='/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/HCP_D/noconf_dkt',
    outdir='/Users/fyzeen/FyzeenLocal/GitHub/FyzWAPIAW2024/regression_outputs/HCP_D/')
