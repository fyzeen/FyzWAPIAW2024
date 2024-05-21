#!/bin/bash -l
#SBATCH --job-name=make_confounds_numeric
#SBATCH --account=janine_bijsterbosch
#SBATCH --output="/scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_ya/logs/make_confounds_numeric.out%j"
#SBATCH --error="/scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_ya/logs/make_confounds_numeric.err%j"
#SBATCH --time=01:55:00
#SBATCH --mem=25gb
#SBATCH --partition=tier2_cpu

conda activate /home/ahmadf/miniconda/envs/WAPIAW2024

python /scratch/l.lexi/WAPIAW2024/Experiments/HCP_team/hcp_ya/code/make_confounds_numeric.py