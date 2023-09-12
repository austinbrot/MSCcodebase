#!/usr/bin/bash
#SBATCH --job-name=batch_msc_analyses
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err
#SBATCH --time=6:00:00
#SBATCH -p normal,russpold
#SBATCH -c 8
#SBATCH --mem=128GB

# Assumes MSC_codebase is in home directory. Change if not.
MSC_CODEBASE_PATH="${HOME}/MSCcodebase"

echo "Running batch MSC analysis from ${MSC_CODEBASE_PATH}"

ml matlab

matlab -r "addpath(genpath('${MSC_CODEBASE_PATH}')); batch_MSC_analyses; exit"

