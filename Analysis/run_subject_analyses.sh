#!/usr/bin/bash
#SBATCH --job-name=msc_subject_analysis
#SBATCH --output=log/%x_subject_%a.out # %a for array task ID
#SBATCH --error=log/%x_subject_%a.err  # %a for array task ID
#SBATCH --time=12:00:00 # Estimate 12 hours per subject, adjust as needed
#SBATCH -p normal,hns,russpold # Partitions for subject jobs
#SBATCH -c 12      # Number of CPUs, adjust based on script needs (e.g., Infomap)
#SBATCH --mem=96G  # Memory per subject job
#SBATCH --mail-type=FAIL,END # Notify on fail or end of entire array
#SBATCH --mail-user=abrotman@stanford.edu # Replace with your email

# Note: The --array directive will be added by the master script when submitting with sbatch

# Environment setup from master script
if [ -z "${HOME_DIR}" ]; then echo "HOME_DIR not set"; exit 1; fi
if [ -z "${MSC_CODEBASE_PATH}" ]; then echo "MSC_CODEBASE_PATH not set"; exit 1; fi
if [ -z "${LOG_DIR}" ]; then echo "LOG_DIR not set"; exit 1; fi
if [ -z "${SUBJECT_ID_LIST_CSV}" ]; then echo "SUBJECT_ID_LIST_CSV not set"; exit 1; fi
if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then echo "SLURM_ARRAY_TASK_ID not set. This script must be run as a SLURM job array task."; exit 1; fi

cd "${MSC_CODEBASE_PATH}/Analysis"

# Ensure log directory exists
mkdir -p "${LOG_DIR}"

# --- Determine current subject ID ---
# Convert CSV to bash array
IFS=',' read -r -a SUBJECT_ID_BASH_ARRAY <<< "${SUBJECT_ID_LIST_CSV}"

# SLURM_ARRAY_TASK_ID is 1-based, bash arrays are 0-based
ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

if [ ${ARRAY_INDEX} -lt 0 ] || [ ${ARRAY_INDEX} -ge ${#SUBJECT_ID_BASH_ARRAY[@]} ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is out of bounds for subject list of length ${#SUBJECT_ID_BASH_ARRAY[@]}."
    exit 1
fi

CURRENT_SUBJECT_ID=${SUBJECT_ID_BASH_ARRAY[ARRAY_INDEX]}

if [ -z "${CURRENT_SUBJECT_ID}" ]; then
    echo "Error: Failed to determine CURRENT_SUBJECT_ID for SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi

echo "Running MSC subject analysis for SLURM Array Task ID: ${SLURM_ARRAY_TASK_ID}, Subject ID: ${CURRENT_SUBJECT_ID}"

# Load modules
ml biology workbench
ml matlab

# Run MATLAB script, passing the current subject ID as an argument
matlab -nodisplay -nosplash -r "try; addpath(genpath(\'${MSC_CODEBASE_PATH}\')); batch_MSC_analyses_BIDS_mod(\'${CURRENT_SUBJECT_ID}\'); catch e; fprintf('MATLAB Error: %s\\n', e.message); for i=1:numel(e.stack), fprintf('File: %s, Name: %s, Line: %d\\n', e.stack(i).file, e.stack(i).name, e.stack(i).line); end; exit(1); end; exit(0);"

MATLAB_EXIT_CODE=$?
if [ ${MATLAB_EXIT_CODE} -ne 0 ]; then
    echo "MATLAB script batch_MSC_analyses_BIDS_mod failed for subject ${CURRENT_SUBJECT_ID} (Array Task ID ${SLURM_ARRAY_TASK_ID}) with exit code ${MATLAB_EXIT_CODE}."
    exit ${MATLAB_EXIT_CODE}
else
    echo "MATLAB script batch_MSC_analyses_BIDS_mod completed successfully for subject ${CURRENT_SUBJECT_ID} (Array Task ID ${SLURM_ARRAY_TASK_ID})."
fi

echo "Subject analysis job for ${CURRENT_SUBJECT_ID} (Array Task ID ${SLURM_ARRAY_TASK_ID}) finished." 