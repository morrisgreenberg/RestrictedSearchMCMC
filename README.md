# RestrictedSearchMCMC

Implementation of our Restricted Search MCMC methods for graph inference.

## Current (Potential) To-Dos:
1. Update `generate_order_plus_table` to use hashing instead of Hadamard products.
2. Test whether using `keras` for predicting weights can improve the algorithm
3. Consider sparsity updating

## Completed Tasks:

- 2023-09-06: Updated `generate_order_table` to be as/more efficient than equivalent BiDAG implementation
