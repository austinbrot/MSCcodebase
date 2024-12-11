#!/usr/bin/bash
#SBATCH --job-name=msc_reliability
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err
#SBATCH --time=24:00:00
#SBATCH -p normal
#SBATCH -c 14
#SBATCH --mem 64G
#SBATCH --mail-type END
#SBATCH --mail-user abrotman@stanford.edu

# Assumes MSC_codebase is in home directory. Change if not.
MSC_CODEBASE_PATH="${HOME}/MSCcodebase"

echo "Running reliability MSC analysis from ${MSC_CODEBASE_PATH}"

ml biology workbench
ml matlab

matlab -batch "addpath(genpath('${MSC_CODEBASE_PATH}')); MSC_reliability; exit"

