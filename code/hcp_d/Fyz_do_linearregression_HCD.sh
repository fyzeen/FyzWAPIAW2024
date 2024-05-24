#!/bin/bash -l
#SBATCH --job-name=fy_HCD
#SBATCH --account=janine_bijsterbosch
#SBATCH --output="/scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_d/logs/fy_HCD.out%j"
#SBATCH --error="/scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_d/logs/fy_HCD.err%j"
#SBATCH --time=01:55:00
#SBATCH --mem=10gb
#SBATCH --partition=tier2_cpu

# activate env
conda activate /home/ahmadf/miniconda/envs/WAPIAW2024

# constants (double-check the dataset path!)
base_dir="/scratch/l.lexi/WAPIAW2024"
exp_dir="${base_dir}/Experiments/HCP_team"
dataset_dir="${base_dir}/Data/HCP_Development"

# path to regression code, #${exp_dir}/hcp_d/regression_hcd.py #"${base_dir}/Source_Code/regression.py"
reg_script="${base_dir}/Source_Code/regression.py"

# input parameters
IDP_filepath="${dataset_dir}/dkt_yeo.csv" #dkt_yeo.csv, aseg_yeo.csv, /func_idp/Netmats/
phenotype_filepath="${dataset_dir}/phenotypes.csv"
confound_filepath="${exp_dir}/hcp_d/confounds_num.csv"
outdir="${exp_dir}/hcp_d/Results"

# optional arguments (subject lists): if you use these, they will be paths to *your* files in *your* directories!
sublist_filepath="${dataset_dir}/subjects_list.csv"
ptylist_filepath="${exp_dir}/use_phenotypes.csv"
cnflist_filepath="${exp_dir}/hcp_d/use_confounds.csv"

## sample python call with all arguments specified (except for IDP subset list; likely to be unusued), to be used as a template for later calls
# python ${reg_script} -x "${IDP_filepath}" -y "${phenotype_filepath}" -c "${confound_filepath}" -s "${sublist_filepath}" -p "${ptylist_filepath}" -C "${cnflist_filepath}" -o "${outdir}" -N "${outname}" -v

## actual python call, with only arguments of interest specified, 
echo Reg file found in: ${reg_script}
# with no conf
#outname="linreg_noconf_dkt"
#python ${reg_script} -x "${IDP_filepath}" -y "${phenotype_filepath}" -s "${sublist_filepath}" -o "${outdir}" -N "${outname}" -v
# with conf
outname="linreg_siteconf_dkt" # change to dkt or aseg
python ${reg_script} -x "${IDP_filepath}" -y "${phenotype_filepath}" -s "${sublist_filepath}" -o "${outdir}" -N "${outname}" -v

conda activate