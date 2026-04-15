# RestrictedSearchMCMC

Implementation of our Restricted Search MCMC methods for graph inference, Birth-death processes Restricted Over Order Distributions (BROOD).

## Standard Usage (Local R Session)

Most simply, to run BROOD for a specific graph problem, you should call the `graph_mcmc()` function in `./Scripts/BROOD_Functions.R`. Please first source the script to load all helper functions into R, and ensure you have the following package dependencies installed: `Matrix`, `gtools`, `Rfast`, `BiDAG`.

```R
source("./Scripts/BROOD_Functions.R")
```

By default, BROOD uses BGe scoring (Heckerman and Geiger, 1995) as natively provided in `BiDAG`. To perform DAG-Wishart (Ben-David et al, 2015) scoring, please source the data generation script as well:

```R
source("./Scripts/Data_Generation_Functions.R")
```

**Note**: This allows a user to access `usrDAGcorescore` that points to our DAG-Wishart implementation, which overwrites the native `usrDAGcorescore` in the `BiDAG` R package (Suter and Kuipers). Its purpose is to act as a wrapper for creating the data structure used to perform restricted graph MCMC in `BiDAG` based on a user-specified scoring function.

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
