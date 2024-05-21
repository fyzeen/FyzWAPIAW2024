#!/bin/bash -l
#SBATCH --job-name=HCP_YA_dktYeo_allConfounds
#SBATCH --account=janine_bijsterbosch
#SBATCH --output="/scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_ya/logs/HCP_YA_dktYeo_allConfounds.out%j"
#SBATCH --error="/scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_ya/logs/HCP_YA_dktYeo_allConfounds.err%j"
#SBATCH --time=01:55:00
#SBATCH --mem=25gb
#SBATCH --partition=tier2_cpu

conda activate /home/ahmadf/miniconda/envs/WAPIAW2024

# constants (double-check the dataset path!)
base_dir="/scratch/l.lexi/WAPIAW2024"
exp_dir="${base_dir}/Experiments/HCP_team/hcp_ya"
dataset_dir="${base_dir}/Data/HCP_YA"

# path to regression code
reg_script="${base_dir}/Source_Code/regression.py"

# input parameters
IDP_filepath="${dataset_dir}/dkt_yeo.csv"
phenotype_filepath="${dataset_dir}/phenotypes.csv"
confound_filepath="${dataset_dir}/confounds_noFamilyID.csv"
outdir="${exp_dir}/results"
outname="HCP_YA_dktYeo_allConfounds"

# optional arguments (subject lists): if you use these, they will be paths to *your* files in *your* directories!
sublist_filepath="${exp_dir}/HCP_YA_CleanedSubjList.csv" # can also be .txt
ptylist_filepath="${exp_dir}/use_phenotypes.csv"
cnflist_filepath="${exp_dir}/use_confounds.csv"



## sample python call with all arguments specified (except for IDP subset list; likely to be unusued), to be used as a template for later calls
# python ${reg_script} -x "${IDP_filepath}" -y "${phenotype_filepath}" -c "${confound_filepath}" -s "${sublist_filepath}" -p "${ptylist_filepath}" -C "${cnflist_filepath}" -o "${outdir}" -N "${outname}" -v

## actual python call, with only arguments of interest specified
python ${reg_script} -x "${IDP_filepath}" -y "${phenotype_filepath}" -c "${confound_filepath}" -C "${cnflist_filepath}" -s "${sublist_filepath}" -o "${outdir}" -N "${outname}" -v
