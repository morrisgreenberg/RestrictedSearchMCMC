#!/bin/bash
#
# generate_tasks.sh
# Purpose: Generates the tasks.txt file containing parameter combinations
# for the brood.slurm array job.

################################################################################
# USER CONFIGURATION
################################################################################
# This should match the PROJECT_ROOT defined in your brood.slurm script
PROJECT_ROOT="/path/to/your/project_folder"
################################################################################

# Validation
if [[ "$PROJECT_ROOT" == "/path/to/your/"* ]]; then
    echo "ERROR: Please update PROJECT_ROOT in generate_tasks.sh to match your local setup."
    exit 1
fi

task_file="${PROJECT_ROOT}/tasks.txt"

# Create the project directory if it doesn't exist
mkdir -p "$PROJECT_ROOT"

# Define parameter ranges
graph_models=("ERs" "SBM_2" "hSBM_3")
nodes=(20 80 140 200)
seeds=({101..125})

# Clear the task file before starting
> "$task_file"

# Nested loops to generate all combinations (Cross-product)
for seed in "${seeds[@]}"; do
    for gm in "${graph_models[@]}"; do
        for n in "${nodes[@]}"; do
            echo "$seed $gm $n" >> "$task_file"
        done
    done
done

echo "Successfully generated $(wc -l < "$task_file") tasks."
echo "Task file location: $task_file"