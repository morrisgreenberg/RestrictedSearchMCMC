#!/usr/bin/env Rscript

# 1. Parse Arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 10) {
  stop("Usage: Rscript step_2_simulation_worker.R --seed [s] --model [m] --size [n] --root [path_project] --script_root [path_scripts]")
}

seed         <- as.numeric(args[2])
graph_model  <- args[4]
nodes        <- as.numeric(args[6])
proj_root    <- args[8]
script_root  <- args[10]


print(paste("Running simulation with Seed:", seed, "Model:", graph_model, "Nodes:", nodes))

# 2. Libraries
library(tidyverse)
library(BiDAG)
library(matrixStats)
library(Matrix)
library(gtools)
library(Rfast)
library(pcalg)
library(RBGL)
library(foreach)
library(doParallel)
library(igraph)
library(yardstick)
library(BDgraph)

# 3. Path Management

# Set up local paths relative to the project root
source(file.path(script_root, "BROOD_Functions.R"))
source(file.path(script_root, "Data_Generation_Functions.R"))
source(file.path(script_root, "Simulation_Pipeline", "step_3_simulation_loop.R"))

# Use an 'Output' folder inside the project directory for portability
sim_root_output <- file.path(proj_root, "Output")
if(!dir.exists(sim_root_output)){
  dir.create(sim_root_output, recursive=TRUE)
}

# 4. Simulation Parameters
endpoint <- ifelse(nodes < 100, 4, 1)
replicates <- 1:endpoint
pc_alg_thresh <- min(0.4, 20/nodes)
dataset_multiplier <- c(0.5, 2, 10)
data_model <- c("Gaussian", "FCM")
Bs <- 25000
save_waiting_times <- FALSE
thinning_rate_full <- 2500
burn_in_rate <- 0.1
adapt_space_prob <- 0.1
max_spars <- c(max(12, round(3 + 0.06 * nodes)), max(10, round(3 + 0.05 * nodes)))
max_change <- 1

# 5. Logging / Debugging
debug_filename <- paste0("debug_", seed, "_", graph_model, "_", nodes, ".txt")
debug_file     <- file.path(proj_root, "logs", debug_filename)

# 6. Parallel Setup
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
cl <- makeCluster(max(1, n_cores - 1), "PSOCK", outfile = debug_file)

# Export variables and functions to the cluster
clusterExport(cl, varlist=ls())
clusterEvalQ(cl, {
  library(tidyverse); library(BiDAG); library(matrixStats); library(gtools)
  library(Rfast); library(pcalg); library(RBGL); library(pROC)
  library(igraph); library(yardstick); library(Matrix); library(BDgraph)
})

registerDoParallel(cl)

# 7. Execution
print("Job about to start...")
output <- parallelized_job(seed, graph_model, nodes, sim_root_output,
                           pc_alg_thresh, dataset_multiplier,
                           data_model, replicates, 
                           thinning_rate_full, burn_in_rate, max_spars, 
                           max_change, adapt_space_prob, save_waiting_times,
                           Bs)
print("Job complete!")
stopCluster(cl)

# 8. Save Output
summary_dir <- file.path(sim_root_output, "summaries")
if(!dir.exists(summary_dir)) dir.create(summary_dir)

file_all <- paste0("gm_", graph_model, "_n_", nodes, "_seed_", seed, ".Rdata")
save(output, file = file.path(summary_dir, file_all))

print(paste("Saved output to:", file.path(summary_dir, file_all)))