# RestrictedSearchMCMC

Implementation of our Restricted Search MCMC methods for graph inference, Birth-death processes Restricted Over Order Distributions (BROOD).

## Standard Usage (Local R Session)

Most simply, to run BROOD for a specific graph problem, you should call the `graph_mcmc()` function in `./Scripts/BROOD_Functions.R`. Please first source the script to load all helper functions into R, and ensure you have the following package dependencies installed: `Matrix`, `gtools`, `Rfast`, `BiDAG`, `Rcpp` and `RcppArmadillo`. A working C++ compiler is also required, since `BROOD_Functions.R` compiles a small Rcpp routine (`Scripts/score_table.cpp`) on load for the BGe scoring path.

```R
source("./Scripts/BROOD_Functions.R")
```

**Note**: `score_table.cpp` must be present in `Scripts/`, alongside `BROOD_Functions.R`, and this script assumes it is being sourced with the repository root as the working directory (matching the usage shown above).

By default, BROOD uses BGe scoring (Heckerman and Geiger, 1995) as natively provided in `BiDAG`. To perform DAG-Wishart (Ben-David et al, 2015) scoring, please source the data generation script as well:

```R
source("./Scripts/Data_Generation_Functions.R")
```

**Note**: This allows a user to access `usrDAGcorescore` that points to our DAG-Wishart implementation, which overwrites the native `usrDAGcorescore` in the `BiDAG` R package (Suter and Kuipers). Its purpose is to act as a wrapper for creating the data structure used to perform restricted graph MCMC in `BiDAG` based on a user-specified scoring function.

### Example for 50-node Graph Inference

The example below illustrates how to run our work on a probabilistically generated graph via the Erdős–Rényi model with accompanying synthetic data generated from a Gaussian structural equation model. Please install `pcalg`, `Rfast`, and `yardstick` to generate the graph and see the AUC and F1 performance.

```R
library(pcalg)
library(Rfast)
library(yardstick)
source("./Scripts/BROOD_Functions.R")
source("./Scripts/Data_Generation_Functions.R")

set.seed(4321)
N <- 50
n_samples <- 100

# ---- 1. Generate a random ground-truth DAG ----
trueDAGedges_50 <- as(pcalg::randDAG(n = N, d = 4, wFUN = list(runif, min = 0.4, max = 2)), "matrix")

# spectral_rescale() rescales edge weights so the resulting SEM has
# a stable, well-conditioned covariance structure.
trueDAGedges_curr_50 <- spectral_rescale(trueDAGedges_50)

# ---- 2. Generate SEM data from the true DAG ----
data_50 <- rmvDAG(trueDAGedges_curr_50, n_samples, standardise = FALSE)

# ---- 3. Build the initial (restricted) search space via the PC algorithm ----
pc_fit_50 <- pcalg::skeleton(suffStat = list(C = cor(data_50), n = nrow(data_50)),
                             indepTest = pcalg::gaussCItest, alpha = 0.4,
                             labels = paste0("V", 1:ncol(data_50)), method = "stable",
                             m.max = 10)
space_PC_50 <- 1 * as(pc_fit_50@graph, "matrix")

# ---- 4. Get BGe score parameters ----
score_par_sim_50 <- BiDAG::scoreparameters("bge", data_50)

# ---- 5. Run BROOD ----

# space_move_prob=0 means the search space doesn't adapt; equivalent to hybrid order MCMC
no_adaptation_50 <- graph_mcmc(space_PC_50, score_par_sim_50, iter = 40000, space_move_prob=0,
                               thinning=20, max_sparsity = 15, temper=0.05, verbose = TRUE)

result_50 <- graph_mcmc(space_PC_50, score_par_sim_50, iter = 40000, space_move_prob=0.1,
                        thinning=20, max_sparsity = 15, temper=0.05, verbose = TRUE)

# ==============================================================================
# Useful outputs
# ==============================================================================

# ---- Edge-probability metrics against the known ground truth ----
trueDAG_50 <- (trueDAGedges_curr_50 > 0) * 1
t_vec <- as.numeric(trueDAG_50)
t_f <- factor(t_vec, levels = c(1, 0))

probs_noadapt_50 <- t(as.matrix(Reduce('+', no_adaptation_50$graphs))) / length(no_adaptation_50$graphs)
p_vec_noadapt <- as.numeric(probs_noadapt_50)
p_f_noadapt <- factor(1 * (p_vec_noadapt > 0.5), levels = c(1, 0))

probs_adapt_50 <- t(as.matrix(Reduce('+', result_50$graphs))) / length(result_50$graphs)
p_vec_adapt <- as.numeric(probs_adapt_50)
p_f_adapt <- factor(1 * (p_vec_adapt > 0.5), levels = c(1, 0))

cat("AUC:            ", Rfast::auc(t_vec, p_vec_noadapt), "\n")
cat("PR-AUC:         ", yardstick::pr_auc_vec(t_f, p_vec_noadapt), "\n")
cat("F1:             ", yardstick::f_meas_vec(t_f, p_f_noadapt), "\n")
cat("Mean true-edge probability: ", sum(p_vec_noadapt[t_vec == 1]) / max(1, sum(t_vec == 1)), "\n")
cat("Mean false-edge probability:", sum(p_vec_noadapt[t_vec == 0]) / max(1, sum(t_vec == 0)), "\n")

cat("AUC:            ", Rfast::auc(t_vec, p_vec_adapt), "\n")
cat("PR-AUC:         ", yardstick::pr_auc_vec(t_f, p_vec_adapt), "\n")
cat("F1:             ", yardstick::f_meas_vec(t_f, p_f_adapt), "\n")
cat("Mean true-edge probability: ", sum(p_vec_adapt[t_vec == 1]) / max(1, sum(t_vec == 1)), "\n")
cat("Mean false-edge probability:", sum(p_vec_adapt[t_vec == 0]) / max(1, sum(t_vec == 0)), "\n")


# ---- Search space improvement ----
# Starting search space compared to the true graph
table(trueDAG_50, t(as.matrix(result_50$spaces[[1]])))

# Updated search space by the end
table(trueDAG_50, t(as.matrix(result_50$spaces[[length(result_50$spaces)]])))

# ---- Search space size progression ----
# |E_H| at each thinned snapshot -- shows how much the search space grew
# beyond its PC-algorithm starting point over the course of the run.
space_sizes <- sapply(result_50$spaces, function(H) sum(as.matrix(H)))
iters_saved <- seq_along(space_sizes)  # one entry per thinned sample

plot(iters_saved, space_sizes, type = "l", lwd = 2, col = "steelblue",
     xlab = "Thinned sample index", ylab = expression("Search space size " * "|" * E[H] * "|"),
     main = "BROOD search space size over the course of the run")
abline(h = sum(space_PC_50), lty = 2, col = "grey40")
legend("bottomright", legend = "PC-algorithm starting size", lty = 2, col = "grey40", bty = "n")

```
Using our birth-death-based space changes improves ROC AUC because it adds and removes edges to the search space throughout the sampler, rather than just relying on the straight hybrid sampler which uses the starting search space throughout.

## Large-Scale Simulation Pipeline (SLURM Cluster)

We also provide a complete SLURM pipeline to execute the sampler on our simulation settings in the `./Scripts/Simulation_Pipeline` directory.

### 1. Configure Paths
Open `Scripts/Simulation_Pipeline/step_0_generate_tasks.slurm` and `step_1_job_script_to_launch_BROOD_on_cluster.slurm`. Update the `PROJECT_ROOT` and `SCRIPT_DIR` variables to match the absolute paths on your cluster.

### 2. Generate the Task Grid
Generate the flat-file database of task permutations. This will create `tasks.txt` in the root of the project:

```bash
bash Scripts/Simulation_Pipeline/step_0_generate_tasks.slurm
```

### 3. Launch Jobs on a Cluster

```bash
sbatch Scripts/Simulation_Pipeline/step_1_job_script_to_launch_BROOD_on_cluster.slurm
```

### 4. Output Structure
 
Once executed, the pipeline will automatically generate the following directories in your project root:

- `/logs/`: Contains `.out` and `.err` files from SLURM, and `debug_*.txt` R worker logs.

- `/Output/`: Contains the full `.Rdata` MCMC trace binaries.
- `/Output/summaries/`: Contains the lightweight, aggregated metrics (ROC, PR, F1, Time) ready for visualization.

## References

### BROOD Paper:
[1] Greenberg, M., Campbell, K., Craiu, R. Restricted Search Space Graph MCMC via Birth-Death Processes. arXiv:2604.10863 [stat].

### Other Related Work Directly Influencing This Directory:

[2] Heckerman, D., Geiger D. Learning Bayesian networks: a unification for discrete and Gaussian domains. In *Proceedings of the Eleventh conference on Uncertainty in artificial intelligence*, UAI’95, pages 274–284, San Francisco, CA, USA. Morgan Kaufmann Publishers Inc.

[3] Ben-David, E., Li, T., Massam, H., and Rajaratnam, B. (2015). High dimensional Bayesian inference for Gaussian directed acyclic graph models. arXiv:1109.4371 [math, stat].

[4] Kuipers, J., Moffa, G., and Heckerman, D. (2014). Addendum on the scoring of Gaussian directed acyclic graphical models. *The Annals of Statistics*, 42(4).

[5] Kuipers, J., Suter, P., and Moffa, G. (2021). Efficient sampling and structure learning of Bayesian networks. *Journal of Computational and Graphical Statistics*, 0(0):1–12.

[6] Suter, P., Kuipers, J., Moffa, G., & Beerenwinkel, N. (2023). Bayesian Structure Learning and Sampling of Bayesian Networks with the R Package BiDAG. *Journal of Statistical Software*, 105(9), 1–31. https://doi.org/10.18637/jss.v105.i09
