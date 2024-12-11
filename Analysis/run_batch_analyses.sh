#!/usr/bin/bash
#SBATCH --job-name=batch_msc_analyses
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err
#SBATCH --time=23:30:00
#SBATCH -p bigmem,normal,hns,russpold
#SBATCH -c 14
#SBATCH --mem 360G
#SBATCH --mail-type END
#SBATCH --mail-user abrotman@stanford.edu

# Assumes MSC_codebase is in home directory. Change if not.
MSC_CODEBASE_PATH="${HOME}/MSCcodebase"

echo "Running batch MSC analysis from ${MSC_CODEBASE_PATH}"

ml biology workbench
ml matlab

matlab -batch "addpath(genpath('${MSC_CODEBASE_PATH}')); batch_MSC_analyses; exit"

