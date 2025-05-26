#!/usr/bin/bash

# Master script to run the two-stage MSC analysis pipeline

# --- Configuration ---
# USER: Define the list of subject IDs to process.
# These should be the full subject names, e.g., MSC01, MSC02, etc.
SUBJECT_ID_LIST=("MSC01") # Example: ("MSC01" "MSC02" "MSC03" "MSC04" "MSC05" "MSC06" "MSC07" "MSC08" "MSC09" "MSC10")
# Convert to a comma-separated string for MATLAB
SUBJECT_ID_LIST_CSV=$(IFS=,; echo "${SUBJECT_ID_LIST[*]}")

# Determine the total number of subjects from the list
TOTAL_SUBJECTS_TO_PROCESS=${#SUBJECT_ID_LIST[@]}

if [ "${TOTAL_SUBJECTS_TO_PROCESS}" -eq 0 ]; then
    echo "Error: SUBJECT_ID_LIST is empty. Please define subjects to process."
    exit 1
fi
echo "Processing subjects: ${SUBJECT_ID_LIST_CSV}"
echo "Total subjects: ${TOTAL_SUBJECTS_TO_PROCESS}"

# Paths - adjust if your structure is different
export HOME_DIR="${HOME}"
export MSC_CODEBASE_PATH="${HOME_DIR}/MSCcodebase" # Assumes MSCcodebase is in home directory
export LOG_DIR="${MSC_CODEBASE_PATH}/Analysis/log" # Centralized log directory

PRECOMP_SCRIPT_NAME="run_precomputation.sh"
SUBJECT_ANALYSIS_SCRIPT_NAME="run_subject_analyses.sh"

# --- Ensure scripts are executable ---
chmod +x "${MSC_CODEBASE_PATH}/Analysis/${PRECOMP_SCRIPT_NAME}"
chmod +x "${MSC_CODEBASE_PATH}/Analysis/${SUBJECT_ANALYSIS_SCRIPT_NAME}"

# --- Create Log Directory ---
echo "Creating log directory: ${LOG_DIR}"
mkdir -p "${LOG_DIR}"

# --- Step 1: Run Precomputation ---
echo "Submitting precomputation job..."
# Pass necessary environment variables, including the subject list CSV
PRECOMP_JOB_ID=$(sbatch \
    --export=ALL,HOME_DIR="${HOME_DIR}",MSC_CODEBASE_PATH="${MSC_CODEBASE_PATH}",LOG_DIR="${LOG_DIR}",SUBJECT_ID_LIST_CSV="${SUBJECT_ID_LIST_CSV}" \
    "${MSC_CODEBASE_PATH}/Analysis/${PRECOMP_SCRIPT_NAME}")

# Extract the job ID number
PRECOMP_JOB_ID=$(echo ${PRECOMP_JOB_ID} | awk '{print $NF}')

echo "Precomputation job submitted with ID: ${PRECOMP_JOB_ID}"
echo "You can monitor it with: squeue -u ${USER}"
echo "Or check logs in: ${LOG_DIR}"

# --- Step 2: Run Subject Analyses as a Job Array ---
if [ -z "${PRECOMP_JOB_ID}" ] || [ "${PRECOMP_JOB_ID}" == "Submitted" ] || ! [[ "${PRECOMP_JOB_ID}" =~ ^[0-9]+$ ]]; then
    echo "Error: Precomputation job ID not captured or invalid: '${PRECOMP_JOB_ID}'. Cannot submit dependent subject analyses."
    echo "Please check SLURM queue and logs for ${PRECOMP_SCRIPT_NAME}."
    exit 1
fi

echo ""
echo "Once precomputation job ${PRECOMP_JOB_ID} completes successfully,"
echo "the subject analysis job array will be submitted."

# Submit subject analysis job array, dependent on precomputation job
# The --array directive specifies tasks from 1 to TOTAL_SUBJECTS_TO_PROCESS
# Pass the SUBJECT_ID_LIST_CSV to the array job script as well
SUBJECT_ARRAY_JOB_ID=$(sbatch \
    --depend=afterok:${PRECOMP_JOB_ID} \
    --array=1-${TOTAL_SUBJECTS_TO_PROCESS} \
    --export=ALL,HOME_DIR="${HOME_DIR}",MSC_CODEBASE_PATH="${MSC_CODEBASE_PATH}",LOG_DIR="${LOG_DIR}",SUBJECT_ID_LIST_CSV="${SUBJECT_ID_LIST_CSV}" \
    "${MSC_CODEBASE_PATH}/Analysis/${SUBJECT_ANALYSIS_SCRIPT_NAME}")

# Extract the job ID number
SUBJECT_ARRAY_JOB_ID=$(echo ${SUBJECT_ARRAY_JOB_ID} | awk '{print $NF}')

echo "Subject analysis job array submitted with ID: ${SUBJECT_ARRAY_JOB_ID}"
echo "It will run for subjects 1 to ${TOTAL_SUBJECTS_TO_PROCESS} (index-wise) after job ${PRECOMP_JOB_ID} succeeds."
echo "Monitor with: squeue -u ${USER}"
echo "Individual subject logs will be in: ${LOG_DIR}"

echo ""
echo "Pipeline submission complete." 