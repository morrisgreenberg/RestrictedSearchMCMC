# RestrictedSearchMCMC

Implementation of our Restricted Search MCMC methods for graph inference.

## Current (Potential) To-Dos:
1. Consider different designs of adapting (and finding optimal parameters for mixing and total variation from full space sampling)
2. Test whether using `keras` for predicting weights can improve the algorithm
3. Consider sparsity updating

## Completed Tasks:

- 2023-09-06: Updated `generate_order_table` to be as/more efficient than equivalent BiDAG implementation
- 2023-09-19: Created new functions `banned_parents_mapping`, `plus_parents_mapping`, `score_plus_space`, `create_banned_plus_parent_table`, `sample_plus_graph`, and modified `sample_graph` to make a faster algorithm
- 2023-09-20: Updated the sampler to incorporate new scoring functions that were recently created, deleted `generate_order_score`, `generate_order_plus_score`, `score_full_space_order`, and `score_full_space_plus_order` as they became defunct after the new functions were implemented
- 2023-09-26: Updated the sampler to pass the plus score and banned score objects throughout as often as possible (instead of inefficiently rerunning them each time an expansion is proposed).
- 2023-10-29: Updated the sampler proposal bouncing to correct form for full support kernels at each transition (by creating `sample_graph_random`, `implement_order_random`, `sample_from_2_graphs`, and updated all other functions to accommodate for these)
