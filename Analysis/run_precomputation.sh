#!/usr/bin/bash
#SBATCH --job-name=msc_precompute_vertex_data
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err
#SBATCH --time=04:00:00 # Estimate 4 hours, adjust as needed
#SBATCH -p bigmem # Use a partition with large memory
#SBATCH -c 8      # Number of CPUs, adjust based on paircorr_mod parallelization or I/O needs
#SBATCH --mem=450G # Request substantial memory, adjust based on largest subject/session count
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=abrotman@stanford.edu # Replace with your email

# Environment setup
if [ -z "${HOME_DIR}" ]; then echo "HOME_DIR not set"; exit 1; fi
if [ -z "${MSC_CODEBASE_PATH}" ]; then echo "MSC_CODEBASE_PATH not set"; exit 1; fi
if [ -z "${LOG_DIR}" ]; then echo "LOG_DIR not set"; exit 1; fi
if [ -z "${SUBJECT_ID_LIST_CSV}" ]; then echo "SUBJECT_ID_LIST_CSV not set"; exit 1; fi

cd "${MSC_CODEBASE_PATH}/Analysis" # Ensure MATLAB runs from the correct directory

# Ensure log directory exists (it should be created by master script)
mkdir -p "${LOG_DIR}"

echo "Starting MATLAB precomputation script for subjects: ${SUBJECT_ID_LIST_CSV}"

# Load modules
ml biology workbench
ml matlab

# Run MATLAB script, passing the subject list CSV as an argument
# Note the single quotes around the MATLAB command and escaped single quotes for the argument string
matlab -nodisplay -nosplash -r "try; addpath(genpath(\'${MSC_CODEBASE_PATH}\')); precompute_vertexwise_data(\'${SUBJECT_ID_LIST_CSV}\'); catch e; fprintf('MATLAB Error: %s\\n', e.message); for i=1:numel(e.stack), fprintf('File: %s, Name: %s, Line: %d\\n', e.stack(i).file, e.stack(i).name, e.stack(i).line); end; exit(1); end; exit(0);"

MATLAB_EXIT_CODE=$?
if [ ${MATLAB_EXIT_CODE} -ne 0 ]; then
    echo "MATLAB script precompute_vertexwise_data failed with exit code ${MATLAB_EXIT_CODE}."
    exit ${MATLAB_EXIT_CODE}
else
    echo "MATLAB script precompute_vertexwise_data completed successfully."
fi

echo "Precomputation job finished." 