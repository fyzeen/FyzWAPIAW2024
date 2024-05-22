#!/bin/bash -l
#SBATCH --job-name=fy_HCA
#SBATCH --account=janine_bijsterbosch
#SBATCH --output="/scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_oa/logs/fy_HCA.out%j"
#SBATCH --error="/scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_oa/logs/fy_HCA.err%j"
#SBATCH --time=00:20:00
#SBATCH --mem=5gb
#SBATCH --partition=tier2_cpu

# activate env
conda activate /home/ahmadf/miniconda/envs/WAPIAW2024

# constants (double-check the dataset path!)
base_dir="/scratch/l.lexi/WAPIAW2024"
exp_dir="${base_dir}/Experiments/HCP_team"
dataset_dir="${base_dir}/Data/HCP_Aging"

# path to regression code, #${exp_dir}/hcp_d/regression_hcd.py #"${base_dir}/Source_Code/regression.py"
reg_script="${base_dir}/Source_Code/regression.py"

# input parameters
IDP_filepath="${dataset_dir}/aseg_yeo.csv" #dkt_yeo.csv, aseg_yeo.csv, /func_idp/Netmats/
phenotype_filepath="${dataset_dir}/phenotypes.csv"
confound_filepath="${exp_dir}/hcp_oa/confounds_num.csv"
outdir="${exp_dir}/hcp_oa/Results"

# optional arguments (subject lists): if you use these, they will be paths to *your* files in *your* directories!
sublist_filepath="${dataset_dir}/subjects_list.csv"
ptylist_filepath="${exp_dir}/use_phenotypes.csv"
cnflist_filepath="${exp_dir}/hcp_oa/use_confounds.csv"

## sample python call with all arguments specified (except for IDP subset list; likely to be unusued), to be used as a template for later calls
# python ${reg_script} -x "${IDP_filepath}" -y "${phenotype_filepath}" -c "${confound_filepath}" -s "${sublist_filepath}" -p "${ptylist_filepath}" -C "${cnflist_filepath}" -o "${outdir}" -N "${outname}" -v

## actual python call, with only arguments of interest specified, 
echo Reg file found in: ${reg_script}
# # with no conf
# outname="linreg_noconf_dkt"
# python ${reg_script} -x "${IDP_filepath}" -y "${phenotype_filepath}" -s "${sublist_filepath}" -o "${outdir}" -N "${outname}" -v
# with conf
outname="linreg_allconf_noSite_aseg"
python ${reg_script} -x "${IDP_filepath}" -y "${phenotype_filepath}" -s "${sublist_filepath}" -c "${confound_filepath}" -C "${cnflist_filepath}" -o "${outdir}" -N "${outname}" -v

conda activate