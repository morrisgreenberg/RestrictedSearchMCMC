# RestrictedSearchMCMC

Implementation of our Restricted Search MCMC methods for graph inference.

## Current (Potential) To-Dos:
1. Update the sampler to incorporate new scoring functions that were recently created
2. Test whether using `keras` for predicting weights can improve the algorithm
3. Consider different designs of adapting
4. Consider sparsity updating

## Completed Tasks:

- 2023-09-06: Updated `generate_order_table` to be as/more efficient than equivalent BiDAG implementation
- 2023-09-20: Created new functions `banned_parents_mapping`, `plus_parents_mapping`, `score_plus_space`, `create_banned_plus_parent_table`, `sample_plus_graph`, and modified `sample_graph` to make a faster algorithm
