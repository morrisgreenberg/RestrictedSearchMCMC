# Restricted Search Space MCMC sampler:

# Restricted Search Space MCMC sampler:



#description: wrapper function to perform the MCMC chain
#parameters:  H_0 - starting search space
#             param - scoring parameter object, constructed from package BiDAG
#             alpha - numerator prior on the shrinkage weights
#             beta - denominator prior on the shrinkage weights
#             lambda - parameter for a Poisson sparsity prior 
#                      (currently on hold)
#             rho - exponent in shrinkage weights, between 0 and 1 
#                   (larger weight means quicker to approach a denominator of
#                    num of total steps instead of num of considered steps)
#             zeta - exponent for tolerance threshold, between 0 and 1
#                    (weight of 1 causes inv-linear decay, less is inv-sublinear)
#             start_epsilon - tolerance at step 1 of the chain
#             bounce - probability we bounce to a random graph in the space in the proposal
#             d - number of graphs to sample at expansion steps
#             thresh - threshold to add edges to the space from the d samples
#             B - number of steps of the chain
#             start_contract - step in the chain where contraction begins
#             bound_contract - bound of minimal number of steps an edge
#                              needs to be considered before removing it
#             blacklist - list of any parent structures that are not allowed
#             move_type - how to propose move types. Three current methods:
#                         a. relocate - always performs node relocation
#                         b. random - node relocation (NR) 1/3 of the time,
#                                     local transposition (LT) 1/3,
#                                     global swap (GS) 1/3
#                         c. Kuipers - NR 6/(t+7), LT t/(t+7), GS 1/(t+7)
#             verbose - prints the step in the chain number if TRUE, and when
#                       shrinking/expanding occurs
graph_mcmc <- function(H_0, param, alpha=1.25, beta=2, lambda=2.5, 
                       rho=1/1000, zeta=0.85, start_epsilon=0.1, 
                       bounce=0.000000005, d_expand=1, d_shrink=1, 
                       thresh=0.000000001, B=25000, 
                       start_contract=20, bound_contract=100,  
                       blacklist=NULL, move_type="relocate", verbose=TRUE){
  N <- nrow(H_0)
  prec_b <- 1:N
  K_b <- round(sqrt(N/2))
  G_b <- matrix(0, nrow=N, ncol=N)
  H_b <- H_0
  epsilon_b <- start_epsilon
  mappings <- parents_mapping(H_0)
  full_scores <- score_full_space(H_0, mappings,param, N)
  banned_scores <- create_banned_parent_table(H_0, mappings,full_scores, N)
  banned_mappings <- banned_parents_mapping(mappings, prec_b, TRUE)
  plus_mappings <- plus_parents_mapping(H_0, 1, mappings, blacklist)
  full_plus_scores <- score_plus_space(H_0, mappings, plus_mappings, param, N, 
                                       has_scores_orig = TRUE, H_scores=full_scores)
  banned_plus_scores <- create_banned_plus_parent_table(H_0, mappings, plus_mappings,
                                                        full_plus_scores$full_list)
  Gs <- array(dim=c(N, N, B))
  skels <- array(dim=c(N, N, B))
  precs <- matrix(nrow=B, ncol=N)
  Hs <- array(dim=c(N, N, B))
  Ks <- numeric(B)
  selected <- matrix(0, nrow=N, ncol=N)
  considered <- matrix(0, nrow=N, ncol=N)
  weights_matr <- array(dim=c(N, N, B))
  weight_vec <- numeric(B)
  
  for(b in 1:B){
    if(verbose){
      if(b %% 1 == 0){
        print(paste("b: ", b))
      }
    }
    sampler_step <- mcmc_sampler_step(prec_b, G_b, H_b, K_b, epsilon_b, 
                                      alpha, beta, selected, considered, 
                                      b, d_expand, d_shrink, banned_scores, 
                                      mappings, banned_mappings, plus_mappings, 
                                      full_scores, banned_plus_scores, 
                                      full_plus_scores, param, lambda, rho, 
                                      bounce, thresh, start_contract, move_type,
                                      bound_contract, blacklist, verbose)
    
    G_b <- sampler_step$G_t_plus1
    Gs[,,b] <- G_b
    skels[,,b] <- (G_b + t(G_b)>0)*1
    H_b <- sampler_step$H_t_plus1
    Hs[,,b] <- H_b
    prec_b <- sampler_step$prec_t_plus1
    precs[b,] <- prec_b
    selected <- sampler_step$s
    considered <- sampler_step$c
    K_b <- sampler_step$K_t_plus1
    Ks[b] <- K_b
    banned_scores <- sampler_step$banned_scores
    full_scores <- sampler_step$order_scores
    mappings <- sampler_step$par_mappings
    banned_plus_scores <- sampler_step$banned_plus_scores
    full_plus_scores <- sampler_step$order_plus_scores
    plus_mappings <- sampler_step$plus_par_mappings
    banned_mappings <- sampler_step$banned_par_mappings
    weights_matr[,,b] <- sampler_step$weights
    epsilon_b <- start_epsilon/(b^zeta)
    weight_vec[b] <- sampler_step$weight
  }
  return(list(orders=precs, graphs=Gs, skeletons=skels, spaces=Hs, 
              sparsity=Ks, weights=weights_matr, weight = weight_vec))
}

#description: 1 step of the MCMC sampler
#parameters:  prec_t - order at step t
#             G_t - graph at step t
#             H_t - search space at step t
#             K_t - sparsity at step t (currently on hold)
#             e_t - shrinkage threshold at step t
#             alpha - numerator prior on shrinkage weights
#             beta - denominator prior on shrinkage weights
#             selected - vector of num times each edge was selected in 1:t
#             considered - vector of num times each edge was considered in 1:t
#             t - step number in the chain
#             d - number of graphs to draw at expansion
#             space_banned_score_list - banned score list for scoring orders
#             map_pars - hash tables for order scoring
#             full_score_list - score lists for every valid parent set
#             param - score parameter object, constructed from package BiDAG
#             lamb - parameter for a Poisson sparsity prior (currently on hold)
#             rho - exponent in shrinkage weights, between 0 and 1 
#                   (larger weight means quicker to approach a denominator of
#                    num of total steps instead of num of considered steps)
#             bounce - probability we bounce to a random graph in the space in the proposal
#             thresh - threshold to add edges to the space from the d samples
#             start_contract - step in the chain where contraction begins
#             move_probs - how to propose move types. Three current methods:
#                          a. relocate - always performs node relocation
#                          b. random - node relocation (NR) 1/3 of the time,
#                                     local transposition (LT) 1/3,
#                                     global swap (GS) 1/3
#                         c. Kuipers - NR 6/(t+7), LT t/(t+7), GS 1/(t+7)
#             bound_contract - bound of minimal number of steps an edge
#                              needs to be considered before removing it
#             verbose - prints the step in the chain number if TRUE, and when
#                       shrinking/expanding occurs
mcmc_sampler_step <- function(prec_t, G_t, H_t, K_t, e_t, alpha, beta, 
                              selected, considered, t, d_expand, 
                              d_shrink, space_banned_score_list, 
                              map_pars, banned_pars, plus_pars,
                              full_score_list, plus_banned_list, 
                              plus_score_list, param, lamb, rho, bounce,
                              thresh, start_contract, move_probs, 
                              bound_contract, blacklist, verbose){
  N <- nrow(H_t)
  l <- 1
  #sampling whether we do any expansion/contraction steps
  birth_rates <- calculate_birth_rate(H_t, plus_banned_list, banned_pars)
  death_rates <- calculate_death_rate(H_t, full_score_list, prec_t, map_pars,
                                      banned_pars, space_banned_score_list)
  w_t <- 1/(sum(birth_rates)+sum(death_rates))
  is_adaption <- sample(c(T, F), size=1, prob=c(1/sqrt(t), 1-1/sqrt(t)))
  is_contraction <- F
  is_expansion <- F
  if(is_adaption){
    process <- sample(c("Birth", "Death"), size=1, prob=c(sum(birth_rates)*w_t, sum(death_rates)*w_t)) 
    if(process=="Birth"){
      is_expansion <- T
    }
    else{
      is_contraction <- T
    }
  }
  
  #sampling move type
  if(move_probs=="Kuipers"){
    move_type <- sample(c("global swap", "local transposition", "node relocation"), size=1,
                        prob=c(1/(t+7), t/(t+7), 6/(t+7)))
  }
  else if(move_probs=="relocate"){
    move_type <- "node relocation"
  }
  else{
    move_type <- sample(c("global swap", "local transposition", "node relocation"), size=1,
                        prob=c(1/3, 1/3, 1/3))
  }
  if(verbose){print(move_type)}
  if(is_adaption){
    if(is_expansion){
      G_set <- vector(mode="list", length=d_expand)
      #create a set of extra graphs to permanently add to the space
      for(i in 1:d_expand){
        temp <- sample_plus_graph(plus_score_list$full_list, prec_t, map_pars,
                                  plus_banned_list, banned_pars, plus_pars)
        G_set[[i]] <- temp
      }
      if(verbose){print("pre-expand")}
      space_output <- expand_search_space(H_t, G_set, thresh)
      H_t_plus1 <- space_output$H_new
      update_nodes <- space_output$updatenodes
      if(verbose){print("post-expand")}
    }
    else{
      G_set <- vector(mode="list", length=d_shrink)
      #create a set of extra graphs to permanently add to the space
      for(i in 1:d_shrink){
        temp <- sample_minus_graph(H_t, full_score_list, prec_t, map_pars,
                                   banned_pars, space_banned_score_list)
        G_set[[i]] <- temp
      }
      if(verbose){print("pre-shrink")}
      space_output <- shrink_search_space_v3(H_t, G_set, thresh)
      H_t_plus1 <- space_output$H_new
      update_nodes <- space_output$updatenodes
      if(verbose){print("post-shrink")}
    }
    
    new_mappings <- parents_mapping(H_t_plus1, N, update_nodes,
                                    TRUE, map_pars)
    new_full_scores <- score_full_space(H_t_plus1, new_mappings,
                                        param, N, update_nodes, 
                                        TRUE, full_score_list)
    
    new_banned_scores <- 
      create_banned_parent_table(H_t_plus1, new_mappings, new_full_scores, 
                                 N,update_nodes, TRUE,
                                 space_banned_score_list)
    new_plus_mappings <- plus_parents_mapping(H_t_plus1, l, new_mappings, blacklist)
    new_plus_scores <- score_plus_space(H_t_plus1, new_mappings, new_plus_mappings, 
                                        param, N, update_nodes, TRUE, new_full_scores,
                                        TRUE, plus_score_list)
    new_banned_plus_scores <- create_banned_plus_parent_table(H_t_plus1, new_mappings,
                                                              new_plus_mappings, 
                                                              new_plus_scores$full_list,
                                                              N, update_nodes, TRUE,
                                                              plus_banned_list)
    #considering the current space, and create new space scores for expansion
    # K_Ht <- max(rowSums(H_t))
    # l <- K_t - K_Ht
    # l <- 1
    # plus_pars <- plus_parents_mapping(H_t, l, map_pars)
    #sample new order
    prec_prime <- implement_order_v2(prec_t, move_type,
                                     new_banned_scores, new_mappings)
    prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                         new_banned_scores, new_mappings)
    banned_pars <- banned_parents_mapping(new_mappings, prec_t_plus1, TRUE)
    graph_t_plus1 <- sample_graph(new_full_scores, prec_t_plus1, new_mappings, 
                                  banned_pars)
    
    
    
  }
  
  else{
    #standard sampling of an order
    prec_prime <- implement_order_v2(prec_t, move_type, 
                                     space_banned_score_list, map_pars)
    prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                         space_banned_score_list, map_pars)
    banned_pars <- banned_parents_mapping(map_pars, prec_t_plus1, TRUE)
    
    #standard inventory of current search space sparsity
    # K_Ht <- max(rowSums(H_t))
    # l <- K_t - K_Ht
    
    #standard sampling of graph
    graph_t_plus1 <- sample_graph(full_score_list, prec_t_plus1, map_pars, banned_pars)
    
    #standard update of the search space
    H_t_plus1 <- H_t
    new_full_scores <- full_score_list
    new_banned_scores <- space_banned_score_list
    new_mappings <- map_pars
    new_plus_mappings <- plus_pars
    new_plus_scores <- plus_score_list
    new_banned_plus_scores <- plus_banned_list
  }
  
  weights_list <- create_weights(selected, considered, H_t, graph_t_plus1, t,
                                 bound_contract, alpha, beta, 
                                 rho)
  selected_new <- weights_list$selected_new
  considered_new <- weights_list$considered_new
  weights_new <- weights_list$weights_new
  K_t_plus1 <- K_t
  
  return(list(prec_t_plus1=prec_t_plus1, G_t_plus1=graph_t_plus1, 
              H_t_plus1=H_t_plus1, K_t_plus1=K_t_plus1, 
              s=selected_new, c=considered_new,
              order_scores = new_full_scores, 
              banned_scores = new_banned_scores,
              par_mappings = new_mappings,
              order_plus_scores = new_plus_scores,
              banned_plus_scores = new_banned_plus_scores,
              banned_par_mappings = banned_pars,
              plus_par_mappings = new_plus_mappings,
              weights = weights_new,
              weight = w_t))
  
}

#description: wrapper function to perform the MCMC chain
#parameters:  H_0 - starting search space
#             param - scoring parameter object, constructed from package BiDAG
#             alpha - numerator prior on the shrinkage weights
#             beta - denominator prior on the shrinkage weights
#             lambda - parameter for a Poisson sparsity prior 
#                      (currently on hold)
#             rho - exponent in shrinkage weights, between 0 and 1 
#                   (larger weight means quicker to approach a denominator of
#                    num of total steps instead of num of considered steps)
#             zeta - exponent for tolerance threshold, between 0 and 1
#                    (weight of 1 causes inv-linear decay, less is inv-sublinear)
#             start_epsilon - tolerance at step 1 of the chain
#             bounce - probability we bounce to a random graph in the space in the proposal
#             d - number of graphs to sample at expansion steps
#             thresh - threshold to add edges to the space from the d samples
#             B - number of steps of the chain
#             start_contract - step in the chain where contraction begins
#             bound_contract - bound of minimal number of steps an edge
#                              needs to be considered before removing it
#             blacklist - list of any parent structures that are not allowed
#             move_type - how to propose move types. Three current methods:
#                         a. relocate - always performs node relocation
#                         b. random - node relocation (NR) 1/3 of the time,
#                                     local transposition (LT) 1/3,
#                                     global swap (GS) 1/3
#                         c. Kuipers - NR 6/(t+7), LT t/(t+7), GS 1/(t+7)
#             verbose - prints the step in the chain number if TRUE, and when
#                       shrinking/expanding occurs
graph_mcmc_old <- function(H_0, param, alpha=1.25, beta=2, lambda=2.5, 
                           rho=1/1000, zeta=0.85, start_epsilon=0.1, 
                           bounce=0.000000005, d=1, thresh=0.000000001, B=25000, 
                           start_contract=20, bound_contract=100,  
                           blacklist=NULL, move_type="relocate", verbose=TRUE){
  N <- nrow(H_0)
  prec_b <- 1:N
  K_b <- round(sqrt(N/2))
  G_b <- matrix(0, nrow=N, ncol=N)
  H_b <- H_0
  epsilon_b <- start_epsilon
  mappings <- parents_mapping(H_0)
  full_scores <- score_full_space(H_0, mappings,param, N)
  banned_scores <- create_banned_parent_table(H_0, mappings,full_scores, N)
  plus_mappings <- plus_parents_mapping(H_0, 1, mappings, blacklist)
  full_plus_scores <- score_plus_space(H_0, mappings, plus_mappings, param, N, 
                                       has_scores_orig = TRUE, H_scores=full_scores)
  banned_plus_scores <- create_banned_plus_parent_table(H_0, mappings, plus_mappings,
                                                        full_plus_scores$full_list)
  Gs <- array(dim=c(N, N, B))
  precs <- matrix(nrow=B, ncol=N)
  Hs <- array(dim=c(N, N, B))
  Ks <- numeric(B)
  selected <- matrix(0, nrow=N, ncol=N)
  considered <- matrix(0, nrow=N, ncol=N)
  weights_matr <- array(dim=c(N, N, B))
  
  for(b in 1:B){
    if(verbose){
      if(b %% 1 == 0){
        print(paste("b: ", b))
      }
    }
    sampler_step <- mcmc_sampler_step(prec_b, G_b, H_b, K_b, epsilon_b, 
                                      alpha, beta, selected, considered, 
                                      b, d, banned_scores, mappings,
                                      full_scores, banned_plus_scores, 
                                      plus_mappings, full_plus_scores,
                                      param, lambda, rho, bounce, 
                                      thresh, start_contract, move_type,
                                      bound_contract, blacklist, verbose)
    
    G_b <- sampler_step$G_t_plus1
    Gs[,,b] <- G_b
    H_b <- sampler_step$H_t_plus1
    Hs[,,b] <- H_b
    prec_b <- sampler_step$prec_t_plus1
    precs[b,] <- prec_b
    selected <- sampler_step$s
    considered <- sampler_step$c
    K_b <- sampler_step$K_t_plus1
    Ks[b] <- K_b
    banned_scores <- sampler_step$banned_scores
    full_scores <- sampler_step$order_scores
    mappings <- sampler_step$par_mappings
    banned_plus_scores <- sampler_step$banned_plus_scores
    full_plus_scores <- sampler_step$order_plus_scores
    plus_mappings <- sampler_step$plus_par_mappings
    weights_matr[,,b] <- sampler_step$weights
    epsilon_b <- start_epsilon/(b^zeta)
  }
  return(list(orders=precs, graphs=Gs, spaces=Hs, sparsity=Ks,
              weights=weights_matr))
}


#description: 1 step of the MCMC sampler
#parameters:  prec_t - order at step t
#             G_t - graph at step t
#             H_t - search space at step t
#             K_t - sparsity at step t (currently on hold)
#             e_t - shrinkage threshold at step t
#             alpha - numerator prior on shrinkage weights
#             beta - denominator prior on shrinkage weights
#             selected - vector of num times each edge was selected in 1:t
#             considered - vector of num times each edge was considered in 1:t
#             t - step number in the chain
#             d - number of graphs to draw at expansion
#             space_banned_score_list - banned score list for scoring orders
#             map_pars - hash tables for order scoring
#             full_score_list - score lists for every valid parent set
#             param - score parameter object, constructed from package BiDAG
#             lamb - parameter for a Poisson sparsity prior (currently on hold)
#             rho - exponent in shrinkage weights, between 0 and 1 
#                   (larger weight means quicker to approach a denominator of
#                    num of total steps instead of num of considered steps)
#             bounce - probability we bounce to a random graph in the space in the proposal
#             thresh - threshold to add edges to the space from the d samples
#             start_contract - step in the chain where contraction begins
#             move_probs - how to propose move types. Three current methods:
#                          a. relocate - always performs node relocation
#                          b. random - node relocation (NR) 1/3 of the time,
#                                     local transposition (LT) 1/3,
#                                     global swap (GS) 1/3
#                         c. Kuipers - NR 6/(t+7), LT t/(t+7), GS 1/(t+7)
#             bound_contract - bound of minimal number of steps an edge
#                              needs to be considered before removing it
#             verbose - prints the step in the chain number if TRUE, and when
#                       shrinking/expanding occurs
mcmc_sampler_step_old <- function(prec_t, G_t, H_t, K_t, e_t, alpha, beta, 
                                  selected, considered, t, d, 
                                  space_banned_score_list, map_pars,
                                  full_score_list, plus_banned_list, 
                                  plus_pars, plus_score_list,
                                  param, lamb, rho, bounce,
                                  thresh, start_contract, move_probs, 
                                  bound_contract, blacklist, verbose){
  N <- nrow(H_t)
  l <- 1
  #sampling whether we do any expansion/contraction steps
  is_contraction <- sample(c(T, F), size=1, prob=c(1/sqrt(t), 1-1/sqrt(t)))
  is_expansion <- sample(c(T, F), size=1, prob=c(1/sqrt(t), 1-1/sqrt(t)))
  is_bounce <- sample(c(T, F), size=1, prob=c(bounce, 1-bounce))
  
  #sampling move type
  if(move_probs=="Kuipers"){
    move_type <- sample(c("global swap", "local transposition", "node relocation"), size=1,
                        prob=c(1/(t+7), t/(t+7), 6/(t+7)))
  }
  else if(move_probs=="relocate"){
    move_type <- "node relocation"
  }
  else{
    move_type <- sample(c("global swap", "local transposition", "node relocation"), size=1,
                        prob=c(1/3, 1/3, 1/3))
  }
  if(verbose){print(move_type)}
  if(is_expansion){
    #considering the current space, and create new space scores for expansion
    # K_Ht <- max(rowSums(H_t))
    # l <- K_t - K_Ht
    # l <- 1
    # plus_pars <- plus_parents_mapping(H_t, l, map_pars)
    #sample new order
    if(!is_bounce){
      prec_prime <- implement_order_v2(prec_t, move_type,
                                       space_banned_score_list, map_pars)
      prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                           space_banned_score_list, map_pars)
      banned_pars <- banned_parents_mapping(map_pars, prec_t_plus1, TRUE)
      graph_t_plus1 <- sample_graph(full_score_list, prec_t_plus1, map_pars, 
                                    banned_pars)
      G_set <- vector(mode="list", length=d)
      #create a set of extra graphs to permanently add to the space
      for(i in 1:d){
        temp <- sample_plus_graph(plus_score_list$full_list, prec_t_plus1, map_pars,
                                  plus_banned_list, banned_pars, plus_pars)
        G_set[[i]] <- temp
      }
      
    }
    else{
      if(verbose){print("bouncing proposed")}
      prec_prime <- implement_order_random(prec_t)
      graph_prime <- sample_graph_random(prec_prime)
      banned_pars_prime <- banned_parents_mapping(map_pars, prec_prime, TRUE)
      banned_pars <- banned_parents_mapping(map_pars, prec_t, TRUE)
      graph_output <- sample_from_2_graphs(G_t, graph_prime, full_score_list, map_pars, param)
      graph_t_plus1 <- graph_output$graph
      if(graph_output$choice == "t"){
        if(verbose){print("bouncing rejected")}
        prec_t_plus1 <- prec_t
        G_set <- vector(mode="list", length=d)
        #create a set of extra graphs to permanently add to the space
        for(i in 1:d){
          temp <- sample_plus_graph(plus_score_list$full_list, prec_t_plus1, map_pars,
                                    plus_banned_list, banned_pars, plus_pars)
          G_set[[i]] <- temp
        }
      }
      else{
        if(verbose){print("bouncing accepted")}
        prec_t_plus1 <- prec_prime
        G_set <- vector(mode="list", length=d+1)
        #create a set of extra graphs to permanently add to the space
        for(i in 1:d){
          temp <- sample_plus_graph(plus_score_list$full_list, prec_t_plus1, map_pars,
                                    plus_banned_list, banned_pars_prime, plus_pars)
          G_set[[i]] <- temp
        }
        G_set[[d+1]] <- graph_t_plus1
      }
    }
    
    
    if(verbose){print("pre-expand")}
    space_output <- expand_search_space(H_t, G_set, thresh)
    H_t_plus1 <- space_output$H_new
    update_nodes <- space_output$updatenodes
    new_mappings <- parents_mapping(H_t_plus1, N, update_nodes,
                                    TRUE, map_pars)
    new_full_scores <- score_full_space(H_t_plus1, new_mappings,
                                        param, N, update_nodes, 
                                        TRUE, full_score_list)
    if(verbose){print("post-expand")}
    new_banned_scores <- 
      create_banned_parent_table(H_t_plus1, new_mappings, new_full_scores, 
                                 N,update_nodes, TRUE,
                                 space_banned_score_list)
    new_plus_mappings <- plus_parents_mapping(H_t_plus1, l, new_mappings, blacklist)
    new_plus_scores <- score_plus_space(H_t_plus1, new_mappings, new_plus_mappings, 
                                        param, N, update_nodes, TRUE, new_full_scores,
                                        TRUE, plus_score_list)
    new_banned_plus_scores <- create_banned_plus_parent_table(H_t_plus1, new_mappings,
                                                              new_plus_mappings, 
                                                              new_plus_scores$full_list,
                                                              N, update_nodes, TRUE,
                                                              plus_banned_list)
    if(verbose){print("post-banned")}
  }
  
  else{
    #standard sampling of an order
    if(!is_bounce){
      prec_prime <- implement_order_v2(prec_t, move_type, 
                                       space_banned_score_list, map_pars)
      prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                           space_banned_score_list, map_pars)
      banned_pars <- banned_parents_mapping(map_pars, prec_t_plus1, TRUE)
      
      #standard inventory of current search space sparsity
      # K_Ht <- max(rowSums(H_t))
      # l <- K_t - K_Ht
      
      #standard sampling of graph
      graph_t_plus1 <- sample_graph(full_score_list, prec_t_plus1, map_pars, banned_pars)
      
      #standard update of the search space
      H_t_plus1 <- H_t
      if((t < start_contract ) | (!is_contraction)){
        new_full_scores <- full_score_list
        new_banned_scores <- space_banned_score_list
        new_mappings <- map_pars
        new_plus_mappings <- plus_pars
        new_plus_scores <- plus_score_list
        new_banned_plus_scores <- plus_banned_list
      }
    }
    else{
      if(verbose){print("bouncing proposed")}
      prec_prime <- implement_order_random(prec_t)
      graph_prime <- sample_graph_random(prec_prime)
      graph_output <- sample_from_2_graphs(G_t, graph_prime, full_score_list, map_pars, param)
      graph_t_plus1 <- graph_output$graph
      if(graph_output$choice == "t"){
        if(verbose){print("bouncing rejected")}
        prec_t_plus1 <- prec_t
        H_t_plus1 <- H_t
        if((t < start_contract ) | (!is_contraction)){
          new_full_scores <- full_score_list
          new_banned_scores <- space_banned_score_list
          new_mappings <- map_pars
          new_plus_mappings <- plus_pars
          new_plus_scores <- plus_score_list
          new_banned_plus_scores <- plus_banned_list
        }
      }
      else{
        if(verbose){print("bouncing accepted")}
        prec_t_plus1 <- prec_prime
        space_output <- expand_search_space(H_t, list(graph_t_plus1), thresh)
        H_t_plus1 <- space_output$H_new
        update_nodes <- space_output$updatenodes
        new_mappings <- parents_mapping(H_t_plus1, N, update_nodes,
                                        TRUE, map_pars)
        new_full_scores <- score_full_space(H_t_plus1, new_mappings,
                                            param, N, update_nodes, 
                                            TRUE, full_score_list)
        new_banned_scores <- 
          create_banned_parent_table(H_t_plus1, new_mappings, new_full_scores, 
                                     N,update_nodes, TRUE,
                                     space_banned_score_list)
        new_plus_mappings <- plus_parents_mapping(H_t_plus1, l, new_mappings, blacklist)
        new_plus_scores <- score_plus_space(H_t_plus1, new_mappings, new_plus_mappings, 
                                            param, N, update_nodes, TRUE, new_full_scores,
                                            TRUE, plus_score_list)
        new_banned_plus_scores <- create_banned_plus_parent_table(H_t_plus1, new_mappings,
                                                                  new_plus_mappings, 
                                                                  new_plus_scores$full_list,
                                                                  N, update_nodes, TRUE,
                                                                  plus_banned_list)
      }
    }
  }
  
  weights_list <- create_weights(selected, considered, H_t, graph_t_plus1, t,
                                 bound_contract, alpha, beta, 
                                 rho)
  selected_new <- weights_list$selected_new
  considered_new <- weights_list$considered_new
  weights_new <- weights_list$weights_new
  if(is_contraction & (t >= start_contract)){
    #perform contraction of the space
    n <- param$N
    set_size <- calculate_set_size(graph_t_plus1, n, N)
    if(is_expansion){
      if(verbose){print("pre-shrink")}
      space_output <- shrink_search_space_v2(H_t_plus1, weights_new, set_size)
      H_t_plus1 <- space_output$H_new
      update_nodes <- space_output$updatenodes
      new_mappings <- parents_mapping(H_t_plus1, N, update_nodes,
                                      TRUE, map_pars)
      new_full_scores <- score_full_space(H_t_plus1, new_mappings,
                                          param, N, update_nodes, 
                                          TRUE, full_score_list)
      if(verbose){print("post-shrink")}
      new_banned_scores <- 
        create_banned_parent_table(H_t_plus1, new_mappings,new_full_scores, 
                                   N,update_nodes, TRUE,
                                   space_banned_score_list)
      new_plus_mappings <- plus_parents_mapping(H_t_plus1, l, new_mappings, blacklist)
      new_plus_scores <- score_plus_space(H_t_plus1, new_mappings, new_plus_mappings, 
                                          param, N, update_nodes, TRUE, new_full_scores,
                                          TRUE, plus_score_list)
      new_banned_plus_scores <- create_banned_plus_parent_table(H_t_plus1, new_mappings,
                                                                new_plus_mappings, 
                                                                new_plus_scores$full_list,
                                                                N, update_nodes, TRUE,
                                                                plus_banned_list)
      if(verbose){print("post-banned")}
    }
    else{
      if(verbose){print("pre-shrink")}
      n <- param$N
      set_size <- calculate_set_size(graph_t_plus1, n, N)
      space_output <- shrink_search_space_v2(H_t, weights_new, set_size)
      H_t_plus1 <- space_output$H_new
      update_nodes <- space_output$updatenodes
      new_mappings <- parents_mapping(H_t_plus1, N, update_nodes,
                                      TRUE, map_pars)
      new_full_scores <- score_full_space(H_t_plus1, new_mappings,
                                          param, N, update_nodes, 
                                          TRUE, full_score_list)
      if(verbose){print("post-shrink")}
      new_banned_scores <- 
        create_banned_parent_table(H_t_plus1, new_mappings,
                                   new_full_scores, N,
                                   update_nodes, TRUE, 
                                   space_banned_score_list)
      new_plus_mappings <- plus_parents_mapping(H_t_plus1, l, new_mappings, blacklist)
      new_plus_scores <- score_plus_space(H_t_plus1, new_mappings, new_plus_mappings, 
                                          param, N, update_nodes, TRUE, new_full_scores,
                                          TRUE, plus_score_list)
      new_banned_plus_scores <- create_banned_plus_parent_table(H_t_plus1, new_mappings,
                                                                new_plus_mappings, 
                                                                new_plus_scores$full_list,
                                                                N, update_nodes, TRUE,
                                                                plus_banned_list)
      if(verbose){print("post-banned")}
    }
    
  }
  
  #sample sparsity
  # min_K <- min(K_t, K_Ht)
  # max_K <- max(K_t, K_Ht)
  # K_t_plus1_r <- exp(K_rel_prob_pois(min_K, abs(l), prec_t_plus1, graph_t_plus1, lamb))
  # r <- min(1, K_t_plus_r)
  # K_t_plus1 <- sample(min_K, max_K, 1, prob=c(1-r, r))
  K_t_plus1 <- K_t
  
  return(list(prec_t_plus1=prec_t_plus1, G_t_plus1=graph_t_plus1, 
              H_t_plus1=H_t_plus1, K_t_plus1=K_t_plus1, 
              s=selected_new, c=considered_new,
              order_scores = new_full_scores, 
              banned_scores = new_banned_scores,
              par_mappings = new_mappings,
              order_plus_scores = new_plus_scores,
              banned_plus_scores = new_banned_plus_scores,
              plus_par_mappings = new_plus_mappings,
              weights = weights_new))
  
}
create_banned_parent_table <- function(H, map_pars, score_list, 
                                       N=ncol(H), updatenodes=1:N,
                                       has_scores_orig=FALSE, 
                                       old_scores=NULL){
  poset_pt<-list(length=N)
  for(i in updatenodes){
    Ni_rows <- nrow(map_pars$idx_pset[[i]])
    Ni_cols <- ncol(map_pars$idx_pset[[i]])
    poset_pt[[i]]<-matrix(NA,nrow=Ni_rows,ncol=Ni_cols)
    offsets<-rep(1,Ni_rows)
    
    if(Ni_rows>1){
      for(j in Ni_rows:2){
        col_sel <- c(1:map_pars$numpars_vec[[i]][j])
        p_nodes <- map_pars$idx_pset[[i]][j, col_sel]
        children <-map_pars$maps[[i]]$backwards[map_pars$maps[[i]]$forward[j]-2^(p_nodes)/2]
        poset_pt[[i]][cbind(children, offsets[children])]<-j
        offsets[children]<-offsets[children]+1
      }
    }
  }
  orderscore<-vector(mode="list", length=N)
  revnumpar_vec<-lapply(map_pars$numpars_vec,rev)
  num_pars <- rowSums(H)
  for(i in 1:N){
    if(!(i %in% updatenodes) & has_scores_orig){
      orderscore[[i]] <- old_scores[[i]]
    }
    else{
      N_psets <- nrow(poset_pt[[i]])
      N_ps <- num_pars[i]
      binom_coef <- choose(N_ps, c(0:N_ps))
      P_local<-vector(mode="numeric",length=N_psets)
      P_local[N_psets] <- score_list[[i]][1,1]
      P_local[1] <- logSumExp(score_list[[i]][,1])
      cutoff <- 1
      if(N_psets > 2){
        for(l in 1:(N_ps-1)){
          cutoff <- cutoff+binom_coef[l]
          for(j in (N_psets-1):min(cutoff, N_psets-1)){
            poset_pnodes <- poset_pt[[i]][j, c(1:revnumpar_vec[[i]][j])]
            p_total <- logSumExp(P_local[poset_pnodes])-log(N_ps-revnumpar_vec[[i]][j]-l+1)
            # max_par <- max(P_local[poset_pnodes])
            # p_total_v2 <- log(sum(exp(P_local[poset_pnodes]-max_par)))+max_par-log(N_ps-revnumpar_vec[[i]][j]-l+1)
            conj_score <- score_list[[i]][map_pars$maps[[i]]$backwards[N_psets-map_pars$maps[[i]]$forward[j]+1],1]
            max_amt <- max(conj_score, p_total)
            P_local[j] <- log(exp(conj_score-max_amt)+exp(p_total-max_amt))+max_amt
          }
        }
      }
      orderscore[[i]] <- as.matrix(P_local)
    }
    
  }
  return(orderscore)
}

#description: scoring the plus space conditional on the current order
#parameters:  H - matrix of edges in H that are compatible with the order
#             map_pars - hash tables for order scoring
#             param - score parameter object, constructed from package BiDAG
#             N - number of nodes
#             updatenodes - nodes where scores need to be updated
#             has_scores_orig - logical if scores have already been recorded
#             H_scores - score_list
score_full_space <- function(H, map_pars, param,
                             N=ncol(H), updatenodes=1:N,
                             has_scores_orig=FALSE, H_scores=NULL){
  score_list <- vector(mode="list", length=N)
  for(i in 1:N){
    combos <- map_pars$par_pset[[i]]
    par_vec <- map_pars$numpars_vec[[i]]
    if(!(i %in% updatenodes) & has_scores_orig){
      score_list[[i]] <- H_scores[[i]]
    }
    else{
      n_parent_sets <- nrow(combos)
      score_matr <- matrix(0, nrow=n_parent_sets, ncol=1)
      for(j in 1:n_parent_sets){
        if(j==1){
          parent_group <- integer(0)
        }
        else{
          parent_group <- combos[j, 1:c(par_vec[j])] 
        }
        score_matr[j, 1] <- bge_score_node(i, parent_group, N, param)
      }
      score_list[[i]] <- score_matr
    }
  }
  return(score_list)
}
score_plus_space <- function(H, map_pars, plus_pars, param,
                             N=ncol(H), updatenodes=1:N,
                             has_scores_orig=FALSE, H_scores=NULL, 
                             has_plus_orig=FALSE, H_plus_scores=NULL){
  tot_scores <- vector(mode="list", length=N)
  score_list <- vector(mode="list", length=N)
  for(i in 1:N){
    if(!(i %in% updatenodes) & has_plus_orig){
      score_list[[i]] <- H_plus_scores$full_list[[i]]
      tot_scores[[i]] <- H_plus_scores$tot_scores_add[[i]]
    }
    else{
      combos <- map_pars$par_pset[[i]]
      par_vec <- map_pars$numpars_vec[[i]]
      plus_combos <- plus_pars$par_pset[[i]]
      nonpar_vec <- plus_pars$numpars_vec[[i]]
      n_parent_sets <- nrow(combos)
      n_nonparent_sets <- nrow(plus_combos)
      if(n_nonparent_sets==1){
        if(has_scores_orig){
          score_matr <- H_scores[[i]]
          score_sums <- logSumExp(score_matr)
        }
        else{
          score_matr <- matrix(0, nrow=n_parent_sets, ncol=1)
          for(j in 1:n_parent_sets){
            if(j==1){
              parent_group <- integer(0)
            }
            else{
              parent_group <- combos[j, 1:c(par_vec[j])]
            }
            score_matr[j,1] <- bge_score_node(i, parent_group, N, param)
            
          }
          score_sums <- logSumExp(score_matr)
        }
        tot_scores[[i]] <- score_sums
        score_list[[i]] <- score_matr
      }
      else{
        score_lists <- vector(mode="list", length=n_nonparent_sets)
        score_sums <- vector(mode="numeric", length=n_nonparent_sets)
        for(k in 1:n_nonparent_sets){
          if(k == 1 & has_scores_orig){
            score_matr <- matrix(0, nrow=n_parent_sets, ncol=1)
            score_matr[, 1] <- H_scores[[i]][, 1]
            score_sums[1] <- logSumExp(score_matr[,1])
            score_lists[[1]] <- score_matr
          }
          else{
            n_nonpar <- nonpar_vec[k]
            nonparent_group <- plus_combos[k, 1:n_nonpar]
            score_matr <- matrix(0, nrow=n_parent_sets, ncol=1)
            for(j in 1:n_parent_sets){
              if(j==1){
                if(k==1){
                  full_set <- integer(0)
                }
                else{
                  full_set <- nonparent_group
                }
              }
              else{
                if(k==1){
                  full_set <- combos[j, 1:c(par_vec[j])] 
                }
                else{
                  parent_group <- combos[j, 1:c(par_vec[j])] 
                  full_set <- c(parent_group, nonparent_group)
                }
              }
              score_matr[j, 1] <- bge_score_node(i, full_set, N, param)
            }
            score_lists[[k]] <- score_matr 
            score_sums[k] <- logSumExp(score_matr[,1])
          }
        }
        score_list[[i]] <- score_lists
        tot_scores[[i]] <- score_sums
      }
    }
  }
  return(list("full_list"=score_list, "tot_scores_add"=tot_scores))
}


create_banned_plus_parent_table <- function(H, map_pars, plus_pars, 
                                            score_plus_list, 
                                            N=ncol(H), updatenodes=1:N,
                                            has_scores_orig=FALSE, 
                                            old_scores=NULL){
  poset_pt<-list(length=N)
  for(i in updatenodes){
    Ni_rows <- nrow(map_pars$idx_pset[[i]])
    Ni_cols <- ncol(map_pars$idx_pset[[i]])
    poset_pt[[i]]<-matrix(NA,nrow=Ni_rows,ncol=Ni_cols)
    offsets<-rep(1,Ni_rows)
    
    if(Ni_rows>1){
      for(j in Ni_rows:2){
        col_sel <- c(1:map_pars$numpars_vec[[i]][j])
        p_nodes <- map_pars$idx_pset[[i]][j, col_sel]
        children <-map_pars$maps[[i]]$backwards[map_pars$maps[[i]]$forward[j]-2^(p_nodes)/2]
        poset_pt[[i]][cbind(children, offsets[children])]<-j
        offsets[children]<-offsets[children]+1
      }
    }
  }
  orderscore<-vector(mode="list", length=N)
  revnumpar_vec<-lapply(map_pars$numpars_vec,rev)
  num_pars <- rowSums(H)
  for(i in 1:N){
    if(!(i %in% updatenodes) & has_scores_orig){
      orderscore[[i]] <- old_scores[[i]]
    }
    else{
      N_psets <- nrow(poset_pt[[i]])
      N_ps <- num_pars[i]
      binom_coef <- choose(N_ps, c(0:N_ps))
      N_outside <- length(plus_pars$par_pset[[i]][,1])
      P_local<-matrix(nrow=N_psets, ncol=N_outside)
      
      for(j in 1:N_outside){
        if(N_outside==1){
          P_local[N_psets, j] <- score_plus_list[[i]][1,1]
          P_local[1, j] <- logSumExp(score_plus_list[[i]][,1])
        }
        else{
          P_local[N_psets, j] <- score_plus_list[[i]][[j]][1,1]
          P_local[1, j] <- logSumExp(score_plus_list[[i]][[j]][,1])
        }
        cutoff <- 1
        if(N_psets > 2){
          for(l in 1:(N_ps-1)){
            cutoff <- cutoff + binom_coef[l]
            for(k in (N_psets-1):min(cutoff, N_psets-1)){
              poset_pnodes <- poset_pt[[i]][k, c(1:revnumpar_vec[[i]][k])]
              p_total <- logSumExp(P_local[poset_pnodes, j])-log(N_ps-revnumpar_vec[[i]][k]-l+1)
              ind_curr <- map_pars$maps[[i]]$backwards[N_psets-map_pars$maps[[i]]$forward[k]+1]
              if(N_outside == 1){
                conj_score <- score_plus_list[[i]][ind_curr,1]
              }
              else{
                conj_score <- score_plus_list[[i]][[j]][ind_curr,1]
              }
              
              max_amt <- max(conj_score, p_total)
              P_local[k, j] <- log(exp(conj_score-max_amt)+exp(p_total-max_amt))+max_amt
            }
          }
        }
      }
      orderscore[[i]] <- P_local
    }
    
  }
  return(orderscore)
}
#description: wrapper function for hash tables 
#             H - current space
#             N - number of nodes
#             updatenodes - nodes where scores need to be updated
parents_mapping <- function(H, N=ncol(H), updatenodes=1:N, has_map=FALSE,
                            old_map=NULL){
  if(has_map){
    par_list <- create_parent_table(H, N, updatenodes, TRUE, 
                                    old_map$par_pset)
    idx_list <- create_parent_table_idx(H, N, updatenodes, TRUE, 
                                        old_map$idx_pset)
  }
  else{
    par_list <- create_parent_table(H)
    idx_list <- create_parent_table_idx(H)
  }
  mapi <- list()
  maps <- vector(mode="list", length=N)
  numpars_vec <- vector(mode="list", length=N)
  par_names <- vector(mode="list", length=N)
  for(i in 1:N){
    if(has_map & !(i %in% updatenodes)){
      numpars_vec[[i]] <- old_map$numpars_vec[[i]]
      maps[[i]] <- old_map$maps[[i]]
      par_names[[i]] <- old_map$par_names[[i]]
    }
    else{
      parent_nodes <- which(H[i,]==1)
      N_psets <- nrow(idx_list[[i]])
      n_psets <- log2(N_psets)
      numpars_vec[[i]] <- rep(c(0:n_psets),choose(n_psets,c(0:n_psets)))
      P_local <- numeric(N_psets)
      P_localinv <- numeric(N_psets)
      P_local[1]<-1
      P_localinv[1]<-1
      if(N_psets>1){
        for(j in 2:N_psets){
          p_nodes <- idx_list[[i]][j,]
          P_local[j] <- sum(2^p_nodes, na.rm = TRUE)/2+1
          P_localinv[P_local[j]] <- j
        }
      }
      mapi$forward<-P_local
      mapi$backwards<-P_localinv
      maps[[i]]<- mapi
      par_names[[i]] <- parent_nodes
    }
  }
  return(list("maps"=maps, "par_pset"=par_list, 
              "idx_pset"=idx_list, "numpars_vec"=numpars_vec,
              "par_names"=par_names))
}
#description: creates hash table by node indices
#parameters:  H - search space
#             updatenodes - nodes where scores need to be updated
create_parent_table_idx <- function(H, N=ncol(H), updatenodes=1:N, 
                                    has_map=FALSE,old_map=NULL){
  par_list <- vector(mode="list", length=N)
  for(i in 1:N){
    if(has_map & !(i %in% updatenodes)){
      par_list[[i]] <- old_map[[i]]
    }
    else{
      parent_nodes <- which(H[i,]==1)
      N_parents <- length(parent_nodes)
      if(N_parents==0){
        par_list[[i]] <- matrix(NA, nrow=1, ncol=1)
      }
      else{
        parent_matrix <- rep(NA, times=N_parents)
        for(r in 1:N_parents){
          new_row <- combinations(N_parents, r)
          if(r < N_parents){
            for(j in 1:(N_parents-r))
              new_row <- cbind(new_row, NA)
          }
          parent_matrix <- rbind(parent_matrix, new_row)
        }
        par_list[[i]] <- parent_matrix
      }
    }
  }
  return(par_list)
}
#description: creates hash table by node names
#parameters:  H - search space
#             updatenodes - nodes where scores need to be updated
create_parent_table <- function(H, N=ncol(H), updatenodes=1:N, 
                                has_map=FALSE,old_map=NULL){
  par_list <- vector(mode="list", length=N)
  for(i in 1:N){
    if(has_map & !(i %in% updatenodes)){
      par_list[[i]] <- old_map[[i]]
    }
    else{
      parent_nodes <- which(H[i,]==1)
      N_parents <- length(parent_nodes)
      if(N_parents == 0){
        par_list[[i]] <- matrix(NA, nrow=1, ncol=1)
      }
      else{
        parent_matrix <- rep(NA, times=N_parents)
        for(r in 1:N_parents){
          new_row <- combinations(N_parents, r, parent_nodes)
          if(r < N_parents){
            for(j in 1:(N_parents-r))
              new_row <- cbind(new_row, NA)
          }
          parent_matrix <- rbind(parent_matrix, new_row)
        }
        par_list[[i]] <- parent_matrix
      }
    }
  }
  return(par_list)
}




# these two functions are directly taking dettwobytwo, DAGcorescore 
# functions found here: 
# https://github.com/cran/BiDAG/blob/c4a79e73c901c598c82002d4420df6996bde393f/R/corescore.R
# The BGE score method can be found in: 
# "ADDENDUM ON THE SCORING OF GAUSSIAN DIRECTED ACYCLIC GRAPHICAL MODELS" by Kuipers, Moffa, Heckerman
dettwobytwo <- function(D) {
  D[1,1]*D[2,2]-D[1,2]*D[2,1]
}

bge_score_node <- function(j,parentnodes,N,param){
  TN<-param$TN
  awpN<-param$awpN
  scoreconstvec<-param$scoreconstvec
  
  lp<-length(parentnodes) #number of parents
  awpNd2<-(awpN-N+lp+1)/2
  A<-TN[j,j]
  switch(as.character(lp),
         "0"={# just a single term if no parents
           corescore <- scoreconstvec[lp+1] -awpNd2*log(A)
         },
         
         "1"={# no need for matrices
           D<-TN[parentnodes,parentnodes]
           logdetD<-log(D)
           B<-TN[j,parentnodes]
           logdetpart2<-log(A-B^2/D)
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             corescore <- corescore - param$logedgepmat[parentnodes, j]
           }
         },
         
         "2"={# can do matrix determinant and inverse explicitly
           # but this is numerically unstable for large matrices!
           # so we use the same approach as for 3 parents
           D<-TN[parentnodes,parentnodes]
           detD<-dettwobytwo(D)
           logdetD<-log(detD)
           B<-TN[j,parentnodes]
           #logdetpart2<-log(A-(D[2,2]*B[1]^2+D[1,1]*B[2]^2-2*D[1,2]*B[1]*B[2])/detD) #also using symmetry of D
           logdetpart2<-log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
           }
         },
         
         {# otherwise we use cholesky decomposition to perform both
           D<-as.matrix(TN[parentnodes,parentnodes])
           choltemp<-chol(D)
           logdetD<-2*log(prod(choltemp[(lp+1)*c(0:(lp-1))+1]))
           B<-TN[j,parentnodes]
           logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
           }
         })
  return(corescore)
}

#description: Updates the order from a move type (proportional to posterior)
#parameters:  prec_t - order at current step
#             move_type - move type for constructing proposals. Types of moves:
#                         1. "local transposition" - swaps adjacent nodes
#                         2. "global swap" - swaps 2 nodes
#                         3. "node relocation" - moves a single node
#             space_banned_score_list - banned score table for scoring orders
#             map_pars - hash table for quick scoring
implement_order_v1 <- function(prec_t, move_type, space_banned_score_list,
                               map_pars){
  N <- length(prec_t)
  #sampling a node which will be moved to create the proposal
  node1 <- sample(1:N, size=1)
  other_node_locs <- other_nodes <- 1:N
  other_nodes <- other_nodes[-node1]
  loc_n1 <- which(prec_t==node1)
  other_node_locs <- other_node_locs[-loc_n1]
  if(move_type=="local transposition"){
    if(loc_n1==1){
      if(N > 2){return(c(prec_t[2], prec_t[1], prec_t[3:N]))}
      if(N==2){return(c(prec_t[2], prec_t[1]))}
      return(prec_t)
    }
    if(loc_n1==N){
      if(N > 2){return(c(prec_t[1:(N-2)], prec_t[N], prec_t[N-1]))}
      return(c(prec_t[2], prec_t[1]))
    }
    
    prec_i_min1 <- prec_t
    prec_i_min1[loc_n1] <- prec_t[loc_n1-1]
    prec_i_min1[loc_n1-1] <- node1
    
    prec_i_plus1 <- prec_t
    prec_i_plus1[loc_n1] <- prec_t[loc_n1+1]
    prec_i_plus1[loc_n1+1] <- node1
    return(sample_from_2_orders(prec_i_min1, prec_i_plus1, 
                                space_banned_score_list, map_pars))
  }
  if(move_type=="global swap"){
    prec_list <- vector(mode="list", length=N-1)
    counter <- 0
    for(i in other_node_locs){
      counter <- counter+1
      prec_list[[counter]] <- prec_t
      prec_list[[counter]][loc_n1] <- prec_t[i]
      prec_list[[counter]][i] <- node1
    }
    return(sample_from_multiple_orders(prec_list, space_banned_score_list,
                                       map_pars))
  }
  if(move_type=="node relocation"){
    prec_list <- vector(mode="list", length=N-1)
    counter <- 0
    for(i in other_node_locs){
      counter <- counter+1
      prec_list[[counter]] <- prec_t
      prec_list[[counter]][i] <- node1
      if(i < loc_n1){
        prec_list[[counter]][(i+1):loc_n1] <- prec_t[i:(loc_n1-1)]
      }
      else{
        prec_list[[counter]][loc_n1:(i-1)] <- prec_t[(loc_n1+1):i]
      }
    }
    return(sample_from_multiple_orders(prec_list, space_banned_score_list,
                                       map_pars))
  }
}


#description: Updates the order from a move type (mimics BiDAG implementation)
#parameters:  prec_t - order at current step
#             move_type - move type for constructing proposals. Types of moves:
#                         1. "local transposition" - swaps adjacent nodes
#                         2. "global swap" - swaps 2 nodes
#                         3. "node relocation" - moves a single node
#             space_banned_score_list - banned score table for scoring orders
#             map_pars - hash table for quick scoring
implement_order_v2 <- function(prec_t, move_type,
                               space_banned_score_list,
                               map_pars){
  N <- length(prec_t)
  #sampling a node which will be moved to create the proposal
  node1 <- sample(1:N, size=1)
  other_node_locs <- other_nodes <- 1:N
  other_nodes <- other_nodes[-node1]
  loc_n1 <- which(prec_t==node1)
  other_node_locs <- other_node_locs[-loc_n1]
  if(move_type=="local transposition"){
    node1 <- sample(1:(N-1), size=1)
    changed_nodes <- c(node1, node1+1)
    prec_tplus1 <- prec_t
    prec_tplus1[changed_nodes] <- prec_t[rev(changed_nodes)]
    return(prec_tplus1)
  }
  if(move_type=="global swap"){
    changed_nodes <- sample(1:N,2,replace=FALSE)
    prec_tplus1 <- prec_t
    prec_tplus1[changed_nodes] <- prec_t[rev(changed_nodes)]
    return(prec_tplus1)
  }
  if(move_type=="node relocation"){
    prec_list <- vector(mode="list", length=N-1)
    counter <- 0
    for(i in other_node_locs){
      counter <- counter+1
      prec_list[[counter]] <- prec_t
      prec_list[[counter]][i] <- node1
      if(i < loc_n1){
        prec_list[[counter]][(i+1):loc_n1] <- prec_t[i:(loc_n1-1)]
      }
      else{
        prec_list[[counter]][loc_n1:(i-1)] <- prec_t[(loc_n1+1):i]
      }
    }
    return(sample_from_multiple_orders(prec_list, space_banned_score_list,
                                       map_pars))
  }
}
#description: Updates the order randomly
#parameters:  prec_t - order at current step
implement_order_random <- function(prec_t){
  return(sample(prec_t, size=length(prec_t)))
}

#description: sampling step from M-H to choose between current order and the 
#             proposal, via the banned parent table
#parameters:  prec_t - order at current step
#             prec_prime - proposed order
#             space_banned_score_list - banned score table for scoring orders
#             map_pars - hash table for quick scoring
sample_from_2_orders <- function(prec_t, prec_prime, 
                                 space_banned_score_list, map_pars){
  N <- length(space_banned_score_list)
  score_t <- 0
  score_prime <- 0
  banned_pars_t <- banned_parents_mapping(map_pars, prec_t)
  index_banned_t <- banned_pars_t$banned_row
  banned_pars_prime <- banned_parents_mapping(map_pars, prec_prime)
  index_banned_prime <- banned_pars_prime$banned_row
  for(i in 1:N){
    # access lookup table for current node
    lookup_table <- space_banned_score_list[[i]]
    score_t <- score_t + lookup_table[index_banned_t[i]]
    score_prime <- score_prime + lookup_table[index_banned_prime[i]]
  }
  r <- min(exp(score_prime-score_t), 1)
  prec_t1 <- sample(c("prime", "t"), size=1, prob=c(r, 1-r))
  if(prec_t1=="t"){return(prec_t)}
  return(prec_prime)
}
#description: extension of the above function for more than 2 orders. 
#             This is useful for implementing proposals that are proportional 
#             to the posterior
#parameters:  prec_list - list of orders to score
#             space_banned_score_list - banned score table for scoring orders
#             map_pars - hash table for quick scoring
sample_from_multiple_orders <- function(prec_list, space_banned_score_list,
                                        map_pars){
  N <- length(space_banned_score_list)
  num_orders <- length(prec_list)
  index_i_list <- vector(mode="list", length=num_orders)
  banned_pars <- vector(mode="list", length=num_orders)
  index_banned <- vector(mode="list", length=num_orders)
  score_vec <- rep(0, num_orders)
  for(j in 1:num_orders){
    prec_curr <- prec_list[[j]]
    banned_pars_j <- banned_parents_mapping(map_pars, prec_curr)
    index_banned[[j]] <- banned_pars_j$banned_row
    for(i in 1:N){
      lookup_table <- space_banned_score_list[[i]]
      score_vec[j] <- score_vec[j] + lookup_table[index_banned[[j]][i]]
    }
  }
  # because we operate on the log scale, we use a log-sum function from 
  # the matrixStats package
  score_all <- logSumExp(score_vec)
  score_vec <- exp(score_vec-score_all)
  prec_drawn <- sample(1:num_orders, size=1, prob=score_vec)
  return(prec_list[[prec_drawn]])
}
#description: keeps track of the row compatible with the current order 
#             for each node's banned parent table (and score table)
#parameters:  map_pars - hash table for quick scoring
#             prec - current_order
#             create_plus_sets - indicator for whether to return all valid rows
#                                in the score table
banned_parents_mapping <- function(map_pars, prec, create_plus_sets=FALSE,
                                   create_minus_sets=FALSE){
  N <- length(map_pars$par_names)
  if(create_plus_sets){valid_parset_maps <- vector(mode="list", length=N)}
  if(create_minus_sets){minus_lookup_rows <- vector(mode="list", length=N)}
  if(create_minus_sets){valid_parset_maps2 <- vector(mode="list", length=N)}
  banned_lookup_row <- vector(length=N)
  
  for(i in 1:N){
    # find where the current node is in the order
    index_i <- which(prec==i)
    # find the banned parent sets 
    # (a.k.a, nodes listed after the current node for each order)
    all_parents <- map_pars$par_names[[i]]
    banned_par_idx <- integer(0)
    if(index_i != 1){banned_par_idx <- which(all_parents %in% prec[1:(index_i-1)])}
    
    # access the relevant row from the banned parent table
    index_banned <- ifelse(length(banned_par_idx)==0 | is.na(banned_par_idx[1]),1,
                           map_pars$maps[[i]]$backwards[sum(2^banned_par_idx)/2+1])
    banned_lookup_row[i] <- index_banned
    if(create_minus_sets){
      N_parents <- length(all_parents)
      banned_row_minus <- vector(length=N_parents)
      for(j in 1:N_parents){
        if(j %in% banned_par_idx){
          banned_par_idx2 <- banned_par_idx
        }
        else{
          if(length(banned_par_idx)==0 | is.na(banned_par_idx[1])){
            banned_par_idx2 <- j
          }
          else{
            banned_par_idx2 <- sort(c(banned_par_idx, j))
          }
        }
        index_banned2 <- ifelse(length(banned_par_idx2)==0 | is.na(banned_par_idx2[1]),1,
                               map_pars$maps[[i]]$backwards[sum(2^banned_par_idx2)/2+1])
        banned_row_minus[j] <- index_banned2
      }
      minus_lookup_rows[[i]] <- banned_row_minus
    }
    if(create_plus_sets | create_minus_sets){
      #finding which indices in the score table are valid per the order
      if(index_banned==1){
        allowed_rows<-c(1:nrow(map_pars$par_pset[[i]]))
      }
      else{
        tablesize <- dim(map_pars$par_pset[[i]])
        if(tablesize[1]==1 || length(banned_par_idx)==tablesize[2]){
          allowed_rows <- c(1)
        }
        else{
          allowed_rows <- c(2:tablesize[1])
          banned_pars <- map_pars$par_names[[i]][banned_par_idx]
          for(j in 1:tablesize[2]){
            banned_rows <- which(map_pars$par_pset[[i]][allowed_rows, j] %in% banned_pars)
            if(length(banned_rows)>0){allowed_rows <- allowed_rows[-banned_rows]}
          }
          allowed_rows<-c(1,allowed_rows)
          if(create_minus_sets){
            ind_list <- vector(mode="list", length=tablesize[2])
            
            for(k in 1:tablesize[2]){
              allowed_rows_k <- allowed_rows
              banned_k <- map_pars$par_names[[i]][k]
              for(j in 1:tablesize[2]){
                banned_rows <- which(map_pars$par_pset[[i]][allowed_rows_k, j] == banned_k)
                if(length(banned_rows)>0){allowed_rows_k <- allowed_rows_k[-banned_rows]}
              }
              ind_list[[k]]<-allowed_rows_k
            }
            if(create_minus_sets){
              valid_parset_maps2[[i]] <- ind_list
            }
          }
          #column-wise search for banned parents to eliminate rows
          
        }
      }
      if(create_plus_sets){
        valid_parset_maps[[i]] <- allowed_rows
      }
    }
    
  }
  if(create_minus_sets){
    if(create_plus_sets){
      return(list("banned_row"=banned_lookup_row, "banned_minus_rows"=minus_lookup_rows, 
                  "valid_pset_rows"=valid_parset_maps, "valid_mset_rows"=valid_parset_maps2))
    }
    return(list("banned_row"=banned_lookup_row, "banned_minus_rows"=minus_lookup_rows,
           "valid_mset_rows"=valid_parset_maps2))
  }
  if(create_plus_sets){
    return(list("banned_row"=banned_lookup_row, "valid_pset_rows"=valid_parset_maps))
  }
  return(list("banned_row"=banned_lookup_row))
}

#description: samples a graph from an order score list
#parameters:  order_score_list - list of parent scores compatible with order
sample_graph <- function(score_list, prec, map_pars, banned_map_pars){
  N <- length(score_list)
  G <- matrix(0, ncol=N, nrow=N)
  for(i in 1:N){
    valid_rows <- banned_map_pars$valid_pset_rows[[i]]
    lookup_table <- score_list[[i]][valid_rows,]
    numpars_vec <- map_pars$numpars_vec[[i]][valid_rows]
    if(numpars_vec[length(numpars_vec)]==1){
      par_sets <- as.matrix(map_pars$par_pset[[i]][valid_rows,], ncol=1)
    }
    else{par_sets <- map_pars$par_pset[[i]][valid_rows,]}
    
    prob_vec <- as.numeric(exp(lookup_table-logSumExp(lookup_table)))
    skel_idx <- sample(1:length(prob_vec), size=1, prob=prob_vec)
    if(skel_idx != 1){
      parents <- par_sets[skel_idx,1:numpars_vec[skel_idx]]
      G[i, parents] <- 1
    }
  }
  return(G)
}
#description: samples a graph randomly
#parameters:  prec - current order
sample_graph_random <- function(prec){
  N <- length(prec)
  G <- matrix(0, nrow=N, ncol=N)
  for(i in 1:N){
    index_i <- which(prec==i)
    if(index_i != N){
      num_graphs_size <- choose(N, 0:(N-index_i))
      num_parents <- sample(0:(N-index_i), size=1, prob = num_graphs_size/sum(num_graphs_size))
      if(num_parents > 0){
        if(index_i == N-1 & num_parents==1){
          G[i, prec[N]] <- 1
        }
        else{
          allowed_parents <- sample(prec[(index_i+1):N], size=num_parents)
          G[i, allowed_parents] <- 1
        }
      }
    }
  }
  return(G)
}

sample_plus_graph <- function(score_plus_list, prec, map_pars,
                              banned_plus_list, banned_map_pars, plus_pars){
  N <- length(score_plus_list)
  G <- matrix(0, nrow=N, ncol=N)
  for(i in 1:N){
    N_plus_sets <- ncol(banned_plus_list[[i]])
    tot_scores <- banned_plus_list[[i]][banned_map_pars$banned_row[i],]
    prob_vec <- as.numeric(exp(tot_scores - logSumExp(tot_scores)))
    out_idx <- sample(1:N_plus_sets, size=1, prob=prob_vec)
    
    valid_rows <- banned_map_pars$valid_pset_rows[[i]]
    if(!is.null(dim(score_plus_list[[i]]))){
      lookup_table <- score_plus_list[[i]][valid_rows,]
    }
    else{
      lookup_table <- score_plus_list[[i]][[out_idx]][valid_rows,]
    }
    numpars_vec <- map_pars$numpars_vec[[i]][valid_rows]
    if(numpars_vec[length(numpars_vec)]==1){
      par_sets <- as.matrix(map_pars$par_pset[[i]][valid_rows,], ncol=1)
    }
    else{par_sets <- map_pars$par_pset[[i]][valid_rows,]}
    prob_vec2 <- as.numeric(exp(lookup_table-logSumExp(lookup_table)))
    skel_idx <- sample(1:length(prob_vec2), size=1, prob=prob_vec2)
    out_parents <- plus_pars$par_pset[[i]][out_idx,1:plus_pars$numpars_vec[[i]][out_idx]]
    if(skel_idx != 1){
      skel_parents <- par_sets[skel_idx,1:numpars_vec[skel_idx]]
      parents <- c(out_parents, skel_parents)
    }
    else{
      parents <- out_parents
    }
    G[i, parents] <- 1
  }
  return(G)
}


sample_minus_graph <- function(H, score_list, prec, map_pars,
                               banned_map_pars, banned_score_list){
  
  N <- length(banned_score_list)
  G <- H
  
  if(length(banned_map_pars$banned_minus_rows)==0){
    banned_map_pars <- banned_parents_mapping(map_pars, prec, create_minus_sets=TRUE)
  }

  for(i in 1:N){
    minus_scores <- banned_score_list[[i]][banned_map_pars$banned_minus_rows[[i]],1]
    orig_score <- banned_score_list[[i]][banned_map_pars$banned_row[i],1]
    if(sum(is.na(minus_scores))>0){
      tot_scores <- orig_score
    }
    else{
      tot_scores <- c(orig_score, minus_scores)
    }
    N_minus_sets <- length(tot_scores)
    prob_vec <- as.numeric(exp(tot_scores - logSumExp(tot_scores)))
    out_idx <- sample(1:N_minus_sets, size=1, prob=prob_vec)
    
    
    par_names <- map_pars$par_names[[i]]
    if(out_idx != 1){
      G[i, par_names[out_idx-1]] <- 0
    }
  }
  return(G)
}



#description: sampling step from M-H to choose between current order and the 
#             proposal, via the banned parent table
#parameters:  G_t - graph at current step
#             G_prime - proposed graph
#             score_list - score table
#             map_pars - hash table for quick scoring
#             param - scoring parameter object, constructed from package BiDAG
sample_from_2_graphs <- function(G_t, G_prime, score_list,
                                 map_pars, param){
  N <- nrow(G_t)
  score_t <- 0
  score_prime <- 0
  for(i in 1:N){
    pars_t <- which(G_t[i,]==1)
    pars_prime <- which(G_prime[i,]==1)
    if(length(pars_t)==0){
      score_t <- score_t + score_list[[i]][1,]
    } 
    else if(min(pars_t %in% map_pars$par_names[[i]])){
      valid_pars <- which(map_pars$par_names[[i]] %in% pars_t)
      valid_row <- map_pars$map[[i]]$backwards[sum(2^(valid_pars)/2)+1]
      score_t <- score_t + score_list[[i]][valid_row,]
    }
    else{
      score_t <- score_t + bge_score_node(i, pars_t, N, param)
    }
    if(length(pars_prime)==0){
      score_prime <- score_prime + score_list[[i]][1,]
    }
    else if(min(pars_prime %in% map_pars$par_names[[i]])){
      valid_pars <- which(map_pars$par_names[[i]] %in% pars_prime)
      valid_row <- map_pars$map[[i]]$backwards[sum(2^(valid_pars)/2)+1]
      score_prime <- score_prime + score_list[[i]][valid_row,]
    }
    else{
      score_prime <- score_prime + bge_score_node(i, pars_prime, N, param)
    }
  }
  r <- min(exp(score_prime-score_t), 1)
  prec_t1 <- sample(c("prime", "t"), size=1, prob=c(r, 1-r))
  if(prec_t1=="t"){return(list("graph"=G_t, "choice"="t"))}
  return(list("graph"=G_prime, "choice"="prime"))
}

plus_parents_mapping <- function(H, plus_amt, map_pars, blacklist=NULL){
  N <- ncol(H)
  plus_list <- vector(mode="list", length=N)
  numpars_list <- vector(mode="list", length=N)
  for(i in 1:N){
    outside_set <- setdiff(1:N, c(i, map_pars$par_names[[i]]))
    if(!is.null(blacklist)){
      outside_set <- setdiff(outside_set, blacklist[[i]])
    }
    N_nonparents <- length(outside_set)
    if(N_nonparents==0){
      plus_list[[i]] <- matrix(NA, nrow=1, ncol=1)
      numpars_list[[i]] <- 0
    }
    else{
      numpars_vec <- vector(mode="numeric", 
                            length=sum(choose(N_nonparents, 0:plus_amt)))
      numpars_vec[1] <- 0
      iter <- 2
      nonparent_matrix <- rep(NA, times=plus_amt)
      for(r in 1:plus_amt){
        num_add <- choose(N_nonparents, r)
        numpars_vec[c(iter:(iter+num_add-1))] <- r
        new_row <- combinations(N_nonparents, r, outside_set)
        if(r < plus_amt){
          for(j in 1:(plus_amt-r))
            new_row <- cbind(new_row, NA)
        }
        nonparent_matrix <- rbind(nonparent_matrix, new_row)
        iter <- iter + num_add
      }
      plus_list[[i]] <- nonparent_matrix
      numpars_list[[i]] <- numpars_vec
    }
  }
  return(list("par_pset"=plus_list, "numpars_vec"=numpars_list))
}
#description: creates the edge weights based on the chain
#parameters:  selected - N*N matrix of # times each edge is selected
#             considered - N*N of # times each edge is considered
#             H - search space (used for considered)
#             new_graph - edges that were just selected (used for selected)
#             t - step in the chain
#             bound_contract - minimal number of steps for contraction
#             alpha - numerator prior
#             beta - denominator prior
#             rho - weighted average exponent
create_weights <- function(selected, considered, H, new_graph,
                           t, bound_contract, alpha, beta, rho){
  selected_new <- selected + new_graph
  considered_new <- considered + H
  non_contract_index <- which(considered_new < bound_contract)
  weights_new <- (selected_new + alpha)/(1/t^(rho)*considered_new + 
                                           t*(1-1/t^(rho)) + beta)
  weights_new[non_contract_index] <- 1
  
  diag(weights_new) <- 0
  return(list("weights_new"=weights_new, "selected_new"=selected_new,
              "considered_new"=considered_new))
}


#description: expands the search space by adding any edges from a drawn set 
#             of graphs that are not currently in the search space
#parameters:  H - current search space
#             Gs - set of graphs to add to search space
#             thresh - threshold for percentage of graphs that contain an edge 
#                      in order to add it to the space
expand_search_space <- function(H, Gs, thresh=0.2){
  
  G_all <- matrix(0, nrow=nrow(H), ncol=ncol(H))
  
  for(i in 1:length(Gs)){
    #G_new <- BiDAG:::dagadj2cpadj(Gs[[i]])
    #G_all <- G_all + G_new
    G_all <- G_all + Gs[[i]]
  }
  G_all <- (G_all/(length(Gs))>thresh)*1
  H_new <- (H + G_all>=1)*1
  
  added_nodes <- H_new - H
  update_nodes <- c(1:nrow(H))[rowSums(added_nodes)>0]
  
  return(list(H_new=H_new, updatenodes = update_nodes))
}

#description: contracts the space based on the weights and threshold 
#parameters:  H_t - search space matrix
#             weight_matrix - set of edge weights matrix
#             threshold - tolerance for keeping an edge weight
shrink_search_space <- function(H_t, weight_matrix, threshold){
  
  H_new <- ((weight_matrix >= threshold) & (H_t))*1
  added_nodes <- H_new - H_t
  update_nodes <- c(1:nrow(H_t))[rowSums(added_nodes)>0]
  return(list(H_new=H_new, updatenodes=update_nodes))
  
}

#description: contracts the space based on the weights and threshold 
#parameters:  H_t - search space matrix
#             weight_matrix - set of edge weights matrix
#             set_size - number of edges for keeping in contracted model
shrink_search_space_v2 <- function(H_t, weight_matrix, set_size){
  
  if(sum(H_t) < set_size){
    return(list(H_new=H_t, updatenodes=integer(0)))
  }
  threshold <- weight_matrix[order(-weight_matrix)[set_size]]
  H_new <- ((weight_matrix >= threshold) & (H_t))*1
  # H_new <- BiDAG:::dagadj2cpadj(H_new)
  added_nodes <- H_new - H_t
  update_nodes <- c(1:nrow(H_t))[rowSums(added_nodes)>0]
  return(list(H_new=H_new, updatenodes=update_nodes))
}


#description: expands the search space by adding any edges from a drawn set 
#             of graphs that are not currently in the search space
#parameters:  H - current search space
#             Gs - set of graphs to keep in search space
#             thresh - threshold for percentage of graphs that contain an edge 
#                      in order to add it to the space
shrink_search_space_v3 <- function(H, Gs, thresh=0.2){
  
  G_all <- matrix(0, nrow=nrow(H), ncol=ncol(H))
  
  for(i in 1:length(Gs)){
    G_all <- G_all + Gs[[i]]
  }
  G_all <- (G_all/(length(Gs))>thresh)*1
  H_new <- (G_all>=1)*1
  
  removed_nodes <- H - H_new
  update_nodes <- c(1:nrow(H))[rowSums(removed_nodes)>0]
  
  return(list(H_new=H_new, updatenodes = update_nodes))
}

calculate_set_size <- function(G_t, n, p){
  q_est <- sum(G_t)*2/(p*(p-1))
  d_est <- max(rowSums(G_t))
  k_est <- (log(log(p))-log(n)+2*log(d_est*n)-2*log(-log(q_est)))/(log(-log(q_est))-log(d_est*n))
  k_est <- ifelse(k_est > 0, k_est, 0)
  h <- 1/8*d_est*(n/log(p))^((1+k_est)/(2+k_est))*3*log(p)/log(2)
  return(floor(h))
}


calculate_birth_rate <- function(H, banned_plus_list, banned_map_pars){
  
  N <- length(banned_plus_list)
  outside_matrix <- 1-H-diag(N)
  update_matrix <- outside_matrix
  for(i in 1:N){
    N_plus_sets <- ncol(banned_plus_list[[i]])
    tot_scores <- banned_plus_list[[i]][banned_map_pars$banned_row[i],]
    B_e <- as.numeric(exp(tot_scores - tot_scores[1]))
    B_e <- ifelse(B_e>1, 1, B_e)
    update_matrix[i,as.numeric(which(outside_matrix[i,]==1))] <- B_e[-1]
  }
  return(update_matrix)
}

calculate_death_rate <- function(H, score_list, prec, map_pars, banned_map_pars,
                                 banned_score_list){
  if(length(banned_map_pars$banned_minus_rows)==0){
    banned_map_pars <- banned_parents_mapping(map_pars, prec, create_minus_sets=TRUE)
  }
  
  N <- nrow(H)
  update_matrix <- H
  for(i in 1:N){
    tot_scores <- banned_score_list[[i]][banned_map_pars$banned_minus_rows[[i]],1]
    orig_score <- banned_score_list[[i]][banned_map_pars$banned_row[i],1]
    D_e <- as.numeric(exp(tot_scores - orig_score))
    D_e <- ifelse(D_e>1, 1, D_e)
    update_matrix[i, as.numeric(which(H[i,]==1))] <- D_e
  }
  return(update_matrix)
}


# creates power set of "s" items. From the following stackOverflow post: 
# https://stackoverflow.com/questions/18715580/algorithm-to-calculate-power-set-all-possible-subsets-of-a-set-in-r
powerset <- function(s){
  N <- length(s)
  if(N==0){return(list(numeric(0)))}
  l <- vector(mode="list",length=2^N) ; l[[1]]=numeric()
  counter <- 1L
  for(x in 1L:N){
    for(subset in 1L:counter){
      counter <- counter+1L
      l[[counter]] <- c(l[[subset]],s[x])
    }
  }
  return(l)
}

# finds the indices of a power set of size 2^(max_pwr) with number of 
# parents being less than or equal to plus_amt. Useful for updating the
# space, if we want to consider sets of size 1:plus_amt to add to the 
# current parent set of a node, we can just create a power-set 
# representation of all other nodes and select those that are in this
# index set
index_finder_plus <- function(plus_amt, max_pwr){
  options <- choose(max_pwr, 1:plus_amt)
  N_options <- sum(options)
  index_list <- rep(NA, N_options)
  power_values <- 2^(0:max_pwr)
  counter <- 1
  for(k in 1:plus_amt){
    bit_idx <- t(comb_n(max_pwr, k))
    for(j in 1:nrow(bit_idx)){
      index_j <- bit_idx[j,]
      curr_val <- sum(power_values[index_j])+1
      index_list[counter] <- curr_val
      counter <- counter + 1
    }
  }
  return(index_list)
}

are_equivalent <- function(g1, g2){
  if(sum(g1) != sum(g2)){
    return(FALSE)
  }
  cp1 <- BiDAG:::dagadj2cpadj(g1)
  cp2 <- BiDAG:::dagadj2cpadj(g2)
  compare <- abs(cp1 - cp2)
  return(max(compare)==0)
}
