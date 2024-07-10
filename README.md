# RestrictedSearchMCMC

Implementation of our Restricted Search MCMC methods for graph inference.

## Current (Potential) To-Dos:
1. Add FCM data generation for simulations
2. Add function to assess score-equivalence between 2 graphs (for the purpose of comparing true graph vs. output). We may rely on the `bnlearn` package for this.

## Completed Tasks:

- 2023-09-06: Updated `generate_order_table` to be as/more efficient than equivalent BiDAG implementation
- 2023-09-19: Created new functions `banned_parents_mapping`, `plus_parents_mapping`, `score_plus_space`, `create_banned_plus_parent_table`, `sample_plus_graph`, and modified `sample_graph` to make a faster algorithm
- 2023-09-20: Updated the sampler to incorporate new scoring functions that were recently created, deleted `generate_order_score`, `generate_order_plus_score`, `score_full_space_order`, and `score_full_space_plus_order` as they became defunct after the new functions were implemented
- 2023-09-26: Updated the sampler to pass the plus score and banned score objects throughout as often as possible (instead of inefficiently rerunning them each time an expansion is proposed).
- 2023-10-29: Updated the sampler proposal bouncing to correct form for full support kernels at each transition (by creating `sample_graph_random`, `implement_order_random`, `sample_from_2_graphs`, and updated all other functions to accommodate for these)
- 2024-02-16: Updated contraction method to select edges based on maximal edge set size rather than user-inputted tolerance threshold
- 2024-03-30: Added `calculate_birth_rate`, `calculate_death_rate`, `sample_minus_graph`, and updated `banned_parents_mapping` to allow for minus set output.
- 2024-04-19: Updated `graph_mcmc` and `mcmc_sampler_step` to allow for the birth-death process updates.
- 2024-04-25: Added `rmvlogDAG` and `rmvlogexpDAG` for adding FMC data generation models, updated looping to account for equivalence.
- 2024-06-16: Added `g2Q` to turn stochastic block model output from the `igraph` package (using `sample_sbm`) DAGs. Pushed code with warm-up of 10*n built in.
- 2024-06-17: Updated `create_banned_parent_table` to use the fast zeta transform method explained in Viinikka et al. "Towards Scalable Bayesian Learning of Causal DAGs" (NeurIPS 2020).
- 2024-07-10: Completed coding overhaul to make plus scoring and banned tables more efficient, via new functions `score_plus_space_new`, `create_banned_plus_parent_table_new`, `logMinusExp`, and `bge_score_plus_parent`.
