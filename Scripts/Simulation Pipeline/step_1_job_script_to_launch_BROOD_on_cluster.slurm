#!/bin/bash
#
# brood.slurm
#
#SBATCH --job-name=sim_%A_%a
#SBATCH --output=logs/sim_%A_%a.out  # Using relative paths for logs
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --time=12-12:00:00 # Reconfigure the following arguments based on your cluster's use policies
#SBATCH --mem=50G
#SBATCH --cpus-per-task=6
#SBATCH --ntasks-per-node=2
#SBATCH --qos=normal
#SBATCH --array=1-300%30

################################################################################
# USER CONFIGURATION: Set these paths to match your local environment
################################################################################
# The root directory where your simulation outputs and tasks.txt will live
PROJECT_ROOT="/path/to/your/project_folder"

# The directory where your R scripts are located (one folder above /Simulation_Pipeline)
SCRIPT_DIR="/path/to/your/scripts_folder"
################################################################################

# Validation: Stop if the user hasn't updated the paths
if [[ "$PROJECT_ROOT" == "/path/to/your/"* ]] || [[ "$SCRIPT_DIR" == "/path/to/your/"* ]]; then
    echo "ERROR: Please update PROJECT_ROOT and SCRIPT_DIR in the SLURM script."
    exit 1
fi

task_file="${PROJECT_ROOT}/tasks.txt"
mkdir -p "${PROJECT_ROOT}/logs"

#### 2) Extract this task’s parameters ####
cd "$PROJECT_ROOT"

read seed gm n < <( sed -n "${SLURM_ARRAY_TASK_ID}p" "$task_file" )

echo "[$(date)] Task ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}: seed=$seed, model=$gm, size=$n"

#### 3) Run your R script ####
Rscript "${SCRIPT_DIR}/Simulation_Pipeline/step_2_simulation_worker.R" \
    --seed "$seed" --model "$gm" --size "$n" --root "$PROJECT_ROOT" --script_root "$SCRIPT_DIR"