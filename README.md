# RestrictedSearchMCMC

Implementation of our Restricted Search MCMC methods for graph inference.

## Current (Potential) To-Dos:
1. Test whether using `keras` for predicting weights can improve the algorithm
2. Consider different designs of adapting
3. Consider sparsity updating

## Completed Tasks:

- 2023-09-06: Updated `generate_order_table` to be as/more efficient than equivalent BiDAG implementation
- 2023-09-19: Created new functions `banned_parents_mapping`, `plus_parents_mapping`, `score_plus_space`, `create_banned_plus_parent_table`, `sample_plus_graph`, and modified `sample_graph` to make a faster algorithm
- 2023-09-20: Updated the sampler to incorporate new scoring functions that were recently created, deleted `generate_order_score`, `generate_order_plus_score`, `score_full_space_order`, and `score_full_space_plus_order` as they became defunct after the new functions were implemented
