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
#             d - number of graphs to sample at expansion steps
#             thresh - threshold to add edges to the space from the d samples
#             B - number of steps of the chain
#             start_contract - step in the chain where contraction begins
#             bound_contract - bound of minimal number of steps an edge
#                              needs to be considered before removing it
#             move_type - how to propose move types. Three current methods:
#                         a. relocate - always performs node relocation
#                         b. random - node relocation (NR) 1/3 of the time,
#                                     local transposition (LT) 1/3,
#                                     global swap (GS) 1/3
#                         c. Kuipers - NR 6/(t+7), LT t/(t+7), GS 1/(t+7)
#             verbose - prints the step in the chain number if TRUE, and when
#                       shrinking/expanding occurs
graph_mcmc <- function(H_0, param, alpha=1.25, beta=2, lambda=2.5, 
                       rho=1/1000, zeta=0.85, start_epsilon=0.1, d=5, 
                       thresh=0.2, B=25000, start_contract=20,
                       bound_contract=100,
                       move_type="relocate", verbose=TRUE){
  N <- nrow(H_0)
  prec_b <- 1:N
  K_b <- round(sqrt(N/2))
  H_b <- H_0
  epsilon_b <- start_epsilon
  mappings <- parents_mapping(H_0)
  full_scores <- score_full_space(H_0, mappings,param, N)
  banned_scores <- create_banned_parent_table(H_0, mappings,full_scores, N)
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
    sampler_step <- mcmc_sampler_step(prec_b, H_b, K_b, epsilon_b, 
                                      alpha, beta, selected, considered, 
                                      b, d, banned_scores, mappings,
                                      full_scores, param, lambda, rho,
                                      thresh, start_contract, move_type,
                                      bound_contract, verbose)
    
    Gs[,,b] <- sampler_step$G_t_plus1
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
    weights_matr[,,b] <- sampler_step$weights
    epsilon_b <- start_epsilon/(b^zeta)
  }
  return(list(orders=precs, graphs=Gs, spaces=Hs, sparsity=Ks,
              weights=weights_matr))
}


#description: 1 step of the MCMC sampler
#parameters:  prec_t - order at step t
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
mcmc_sampler_step <- function(prec_t, H_t, K_t, e_t, alpha, beta, 
                              selected, considered, t, d, 
                              space_banned_score_list, map_pars,
                              full_score_list, param, lamb, rho,
                              thresh, start_contract, move_probs, 
                              bound_contract, verbose){
  N <- nrow(H_t)
  #sampling whether we do any expansion/contraction steps
  is_contraction <- sample(c(T, F), size=1, prob=c(1/sqrt(t), 1-1/sqrt(t)))
  is_expansion <- sample(c(T, F), size=1, prob=c(1/sqrt(t), 1-1/sqrt(t)))
  
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
    l <- 1
    #sample new order
    prec_prime <- implement_order_v2(prec_t, move_type,
                                     space_banned_score_list, map_pars)
    prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                         space_banned_score_list, map_pars)
    
    #sample new graph
    order_score_list <- generate_order_table(prec_t_plus1, H_t, param, 
                                             map_pars, has_scores=TRUE,
                                             H_scores=full_score_list)
    graph_t_plus1 <- sample_graph(order_score_list)
    
    
    #create a set of extra graphs to permanently add to the space
    G_set <- vector(mode="list", length=d)
    order_score_list_new <- 
      generate_order_plus_table(prec_t_plus1, H_t, max(l, 1), 
                                param, map_pars, has_scores=TRUE, 
                                H_scores=full_score_list)
    for(i in 1:d){
      temp <- sample_graph(order_score_list_new$full_list)
      G_set[[i]] <- temp
    }
    if(verbose){print("pre-expand")}
    space_output <- expand_search_space(H_t, G_set, thresh)
    H_t_plus1 <- space_output$H_new
    update_nodes <- space_output$updatenodes
    new_mappings <- parents_mapping(H_t_plus1, N, update_nodes,
                                    TRUE, map_pars)
    N_curr <- ncol(H_t_plus1)
    new_full_scores <- score_full_space(H_t_plus1, new_mappings,
                                        param, N_curr, update_nodes, 
                                        TRUE, full_score_list)
    if(verbose){print("post-expand")}
    new_banned_scores <- 
      create_banned_parent_table(H_t_plus1, new_mappings, new_full_scores, 
                                 N_curr,update_nodes, TRUE,
                                 space_banned_score_list)
    if(verbose){print("post-banned")}
  }
  
  else{
    #standard sampling of an order
    prec_prime <- implement_order_v2(prec_t, move_type, 
                                     space_banned_score_list, map_pars)
    prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                         space_banned_score_list, map_pars)
    order_score_list <- generate_order_table(prec_t_plus1, H_t, param,
                                             map_pars, has_scores=TRUE,
                                             H_scores=full_score_list)
    
    #standard inventory of current search space sparsity
    # K_Ht <- max(rowSums(H_t))
    # l <- K_t - K_Ht
    
    #standard sampling of graph
    graph_t_plus1 <- sample_graph(order_score_list)
    
    #standard update of the search space
    H_t_plus1 <- H_t
    if((t < start_contract ) | (!is_contraction)){
      new_full_scores <- full_score_list
      new_banned_scores <- space_banned_score_list
      new_mappings <- map_pars
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
    if(is_expansion){
      if(verbose){print("pre-shrink")}
      space_output <- shrink_search_space(H_t_plus1, weights_new, e_t)
      H_t_plus1 <- space_output$H_new
      update_nodes <- space_output$updatenodes
      new_mappings <- parents_mapping(H_t_plus1, N, update_nodes,
                                      TRUE, map_pars)
      N_curr <- ncol(H_t_plus1)
      new_full_scores <- score_full_space(H_t_plus1, new_mappings,
                                          param, N_curr, update_nodes, 
                                          TRUE, full_score_list)
      if(verbose){print("post-shrink")}
      new_banned_scores <- 
        create_banned_parent_table(H_t_plus1, new_mappings,new_full_scores, 
                                   N_curr,update_nodes, TRUE,
                                   space_banned_score_list)
      if(verbose){print("post-banned")}
    }
    else{
      if(verbose){print("pre-shrink")}
      space_output <- shrink_search_space(H_t, weights_new, e_t)
      H_t_plus1 <- space_output$H_new
      update_nodes <- space_output$updatenodes
      new_mappings <- parents_mapping(H_t_plus1, N, update_nodes,
                                      TRUE, map_pars)
      N_curr <- ncol(H_t_plus1)
      new_full_scores <- score_full_space(H_t_plus1, new_mappings,
                                          param, N_curr, update_nodes, 
                                          TRUE, full_score_list)
      if(verbose){print("post-shrink")}
      new_banned_scores <- 
        create_banned_parent_table(H_t_plus1, new_mappings,
                                   new_full_scores, N_curr,
                                   update_nodes, TRUE, 
                                   space_banned_score_list)
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
            max_par <- max(P_local[poset_pnodes])
            p_total_v2 <- log(sum(exp(P_local[poset_pnodes]-max_par)))+max_par-log(N_ps-revnumpar_vec[[i]][j]-l+1)
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
#description: generates a table of the valid parent sets compatible
#             with the current order. 
#parameters:  prec - order
#             H - search space
#             param - score parameter object, constructed from package BiDAG
#             map_pars - hash tables for order scoring
#             updatenodes - nodes where scores need to be updated
#             has_scores - logical if scores have already been recorded
#             H_scores - score list
generate_order_table <- function(prec, H, param, map_pars,
                                 updatenodes= 1:length(prec),
                                 has_scores=FALSE, H_scores=NULL){
  N <- length(prec)
  # creating relevant space (using Hadamard product between search space 
  # matrix and the lower left-lefthand triangular matrix that 
  # satisfies the order)
  prec_matr <- lower.tri(matrix(0, nrow=N, ncol=N))*1
  prec_matr2 <- prec_matr
  prec_matr[rev(prec),rev(prec)] <- prec_matr2
  order_H <- H * prec_matr
  if(has_scores){
    order_scores <- vector(mode="list", length=N)
    for(i in 1:N){
      lookup_table <- H_scores[[i]]
      curr_pars <- which(order_H[i,]==1)
      incl_idx <- which(map_pars$par_names[[i]] %in% curr_pars)
      parent_sets <- powerset(incl_idx)
      index_keep <- sapply(1:length(parent_sets), function(j){
        curr_subset <- parent_sets[[j]]
        par_index <- ifelse(length(curr_subset)==0, 1, sum(2^curr_subset)/2+1)
        return(map_pars$maps[[i]]$backwards[par_index])})
      score_matr <- matrix(0, nrow=length(index_keep), ncol=1)
      score_matr[1:length(index_keep), 1] <- lookup_table[index_keep,]
      names(score_matr) <- names(lookup_table)[index_keep]
      order_scores[[i]] <- score_matr
    }
    return(order_scores)
  }
  return(score_full_space_order(order_H, map_pars, param, prec))
}

#description: generates a table of the valid plus parent sets compatible
#             with the current order. 
#parameters:  prec - order
#             H - search space
#             plus_amt - amount of additional parents to allow beyond sets
#                        in the current search space
#             param - score parameter object, constructed from package BiDAG
#             map_pars - hash tables for order scoring
#             updatenodes - nodes where scores need to be updated
#             has_scores - logical if scores have already been recorded
#             H_scores - score list
generate_order_plus_table <- function(prec, H, plus_amt, param, 
                                      map_pars, updatenodes=1:length(prec),
                                      has_scores=FALSE, H_scores=NULL){
  N <- length(prec)
  # creating relevant space (using Hadamard product between search space 
  # matrix and the lower left-lefthand triangular matrix that 
  # satisfies the order)
  prec_matr <- lower.tri(matrix(0, nrow=N, ncol=N))*1
  prec_matr2 <- prec_matr
  prec_matr[rev(prec),rev(prec)] <- prec_matr2
  order_H <- H * prec_matr
  if(has_scores){
    order_scores <- vector(mode="list", length=N)
    for(i in 1:N){
      lookup_table <- H_scores[[i]]
      curr_pars <- which(order_H[i,]==1)
      incl_idx <- which(map_pars$par_names[[i]] %in% curr_pars)
      parent_sets <- powerset(incl_idx)
      index_keep <- sapply(1:length(parent_sets), function(j){
        curr_subset <- parent_sets[[j]]
        par_index <- ifelse(length(curr_subset)==0, 1, sum(2^curr_subset)/2+1)
        return(map_pars$maps[[i]]$backwards[par_index])})
      score_matr <- matrix(0, nrow=length(index_keep), ncol=1)
      score_matr[1:length(index_keep), 1] <- lookup_table[index_keep,]
      names(score_matr) <- names(lookup_table)[index_keep]
      order_scores[[i]] <- score_matr
    }
    return(score_full_space_plus_order(order_H, map_pars, param, prec, 
                                       plus_amt,N, updatenodes, has_scores, 
                                       order_scores))
  }
  return(score_full_space_plus_order(order_H, map_pars, param, prec, plus_amt))
}
#description: scoring the space conditional on the current order
#parameters:  order_H - matrix of edges in H that are compatible with the order
#             map_pars - hash tables for order scoring
#             param - score parameter object, constructed from package BiDAG
#             prec - order
#             N - number of nodes
#             updatenodes - nodes where scores need to be updated
#             has_scores_orig - logical if scores have already been recorded
#             H_scores - score_list
score_full_space_order <- function(order_H, map_pars, param, prec,
                                   N=ncol(order_H), updatenodes=1:N,
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
      parent_indices <- which(order_H[i,]==1)
      n_obs_parent_sets <- 2^(length(parent_indices))
      score_matr <- matrix(0, nrow=n_obs_parent_sets, ncol=1)
      iter <- 1
      for(j in 1:n_parent_sets){
        pass_over <- FALSE
        if(j==1){
          parent_group <- integer(0)
        }
        else{
          parent_group <- combos[j, 1:c(par_vec[j])] 
          if(length(which(!(parent_group %in% parent_indices)))>0){
            pass_over <- TRUE
          }
        }
        if(!pass_over){
          score_matr[iter, 1] <- bge_score_node(i, parent_group, N, param)
          names(score_matr)[iter] <- stri_c(as.character(parent_group), collapse=",")
          iter <- iter+1
        }
        
      }
      score_list[[i]] <- score_matr
    }
  }
  return(score_list)
}
#description: scoring the plus space conditional on the current order
#parameters:  order_H - matrix of edges in H that are compatible with the order
#             map_pars - hash tables for order scoring
#             param - score parameter object, constructed from package BiDAG
#             prec - order
#             plus_amt - amount of additional parents to allow beyond sets
#                        in the current search space
#             N - number of nodes
#             updatenodes - nodes where scores need to be updated
#             has_scores_orig - logical if scores have already been recorded
#             H_scores - score_list
score_full_space_plus_order <- function(order_H, map_pars, param, prec, 
                                        plus_amt, N=ncol(order_H), 
                                        updatenodes=1:N, 
                                        has_scores_orig=FALSE, H_scores=NULL){
  tot_scores <- vector(mode="list", length=N)
  score_list <- vector(mode="list", length=N)
  score_list_v2 <- vector(mode="list", length=N)
  plus_set_list <- vector(mode="list", length=N)
  for(i in 1:N){
    if(!(i %in% updatenodes) & has_scores_orig){
      score_list[[i]] <- H_scores[[i]]
      score_list_v2[[i]] <- H_scores[[i]]
      tot_scores[[i]] <- logSumExp(H_scores[[i]])
      plus_set_list[[i]] <- numeric(0)
    }
    else{
      combos <- map_pars$par_pset[[i]]
      par_vec <- map_pars$numpars_vec[[i]]
      parent_indices <- which(order_H[i,]==1)
      nonparent_indices <- which(order_H[i,]==0)
      loc_i <- which(prec==i)
      nonparent_indices <- setdiff(nonparent_indices, prec[1:loc_i])
      plus_sets <- powerset(nonparent_indices)
      N_nonparents <- length(nonparent_indices)
      size_matches <- index_finder_plus(plus_amt, N_nonparents)
      plus_sets_2 <- plus_sets[size_matches]
      n_parent_sets <- nrow(combos)
      n_obs_parent_sets <- 2^(length(parent_indices))
      n_plus_sets <- length(plus_sets_2)
      if(length(nonparent_indices) == 0){
        if(has_scores_orig){
          score_matr <- H_scores[[i]]
          score_sums <- logSumExp(score_matr)
        }
        else{
          score_matr <- matrix(0, nrow=n_obs_parent_sets, ncol=1)
          iter <- 1
          for(j in 1:n_parent_sets){
            pass_over <- FALSE
            if(j==1){
              parent_group <- integer(0)
            }
            else{
              parent_group <- combos[j, 1:c(par_vec[j])] 
              if(length(which(!(parent_group %in% parent_indices)))>0){
                pass_over <- TRUE
              }
            }
            if(!pass_over){
              score_matr[iter,1] <- bge_score_node(i, parent_group, N, param)
              names(score_matr)[iter] <- stri_c(as.character(parent_group), collapse=",")
              iter <- iter+1
            }
            
            
          }
          score_sums <- logSumExp(score_matr)
        }
        tot_scores[[i]] <- score_sums
        score_list[[i]] <- score_matr
        score_list_v2[[i]] <- score_matr
        plus_set_list[[i]] <- numeric(0)
      }
      else{
        score_matr <- matrix(0, nrow=n_obs_parent_sets, 
                             ncol=n_plus_sets+1)
        score_matr_v2 <- matrix(0, nrow=n_obs_parent_sets, ncol=1)
        score_sums <- vector(mode="numeric", length=n_plus_sets+1)
        if(has_scores_orig){
          score_matr[, 1] <- H_scores[[i]][, 1]
          n_H_parset <- nrow(H_scores[[i]])
          names(score_matr)[1:n_H_parset]<-names(H_scores[[i]])
          score_matr_v2[, 1] <- score_matr[, 1]
          names(score_matr_v2)[1:n_H_parset]<-names(H_scores[[i]])
          score_sums[1] <- logSumExp(score_matr[,1])
        }
        iter <- 1
        for(j in 1:n_parent_sets){
          pass_over <- FALSE
          if(j==1){
            parent_group <- integer(0)
          }
          else{
            parent_group <- combos[j, 1:c(par_vec[j])]
            if(length(which(!(parent_group %in% parent_indices)))>0){
              pass_over <- TRUE
            }
          }
          if(!pass_over){
            if(!has_scores_orig){
              score_matr[iter, 1] <- bge_score_node(i, parent_group, N, param)
              names(score_matr)[iter] <- stri_c(as.character(parent_group), 
                                                collapse=",")
              score_matr_v2[iter, 1] <- score_matr[iter, 1]
              names(score_matr_v2)[iter] <- stri_c(as.character(parent_group), 
                                                   collapse=",")
              score_sums[1] <- logSumExp(score_matr[,1])
            }
            for(k in 1:n_plus_sets){
              plus_group <- plus_sets_2[[k]]
              parent_plus_group <- sort(c(parent_group, plus_group))
              score_matr[iter, k+1] <- 
                bge_score_node(i, c(parent_plus_group), N, param)
              names(score_matr)[(n_obs_parent_sets)*k+iter] <- stri_c(as.character(parent_plus_group), 
                                                                      collapse=",")
            }
            score_matr_v2[iter, 1] <- logSumExp(score_matr[iter,])
            names(score_matr_v2)[iter] <- stri_c(as.character(parent_group), 
                                                 collapse=",")
            iter <- iter+1
          }
          
          
        }
        tot_scores[[i]] <- score_sums
        score_list[[i]] <- score_matr
        score_list_v2[[i]] <- score_matr_v2
        plus_set_list[[i]] <- plus_sets_2
      }
    }
    
  }
  return(list("full_list"=score_list, "condensed_list"=score_list_v2,
              "tot_scores_add"=tot_scores, "plus_sets"=plus_set_list))
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
        names(score_matr)[j] <- stri_c(as.character(parent_group), collapse=",")
      }
      score_list[[i]] <- score_matr
    }
  }
  return(score_list)
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
  for(i in 1:N){
    # access lookup table for current node
    lookup_table <- space_banned_score_list[[i]]
    # find where the current node is in each order
    index_i_t <- which(prec_t==i)
    index_i_prime <- which(prec_prime==i)
    # find the banned parent sets 
    # (a.k.a, nodes listed after the current node for each order)
    all_parents <- map_pars$par_names[[i]]
    banned_pars_t <- ifelse(index_i_t != 1,
                            which(all_parents %in% prec_t[1:(index_i_t-1)]),
                            integer(0))
    banned_pars_prime <- ifelse(index_i_prime != 1,
                                which(all_parents %in% 
                                        prec_prime[1:(index_i_prime-1)]), 
                                integer(0))
    # access the relevant row from the banned parent table, and update 
    # score since we are on the log-scale, we add each node's relevant 
    # score to the full order score
    index_banned_t <- ifelse(is.na(banned_pars_t),1,
                             map_pars$maps[[i]]$backwards[sum(2^banned_pars_t)/2+1])
    index_banned_prime <- ifelse(is.na(banned_pars_prime), 1,
                                 map_pars$maps[[i]]$backwards[sum(2^banned_pars_prime)/2+1])
    score_t <- score_t + lookup_table[index_banned_t]
    score_prime <- score_prime + lookup_table[index_banned_prime]
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
  for(i in 1:N){
    lookup_table <- space_banned_score_list[[i]]
    all_parents <- map_pars$par_names[[i]]
    for(j in 1:num_orders){
      prec_curr <- prec_list[[j]]
      index_i_j <-  which(prec_curr==i)
      index_i_list[[j]] <- index_i_j
      banned_pars[[j]] <- ifelse(index_i_j != 1,
                                 which(all_parents %in% prec_curr[1:(index_i_j-1)]),
                                 integer(0))
      index_banned[[j]] <- ifelse(is.na(banned_pars[[j]]), 1,
                                  map_pars$maps[[i]]$backwards[sum(2^banned_pars[[j]])/2+1])
      score_vec[j] <- score_vec[j] + lookup_table[index_banned[[j]]]
    }
  }
  # because we operate on the log scale, we use a log-sum function from 
  # the matrixStats package
  score_all <- logSumExp(score_vec)
  score_vec <- exp(score_vec-score_all)
  prec_drawn <- sample(1:num_orders, size=1, prob=score_vec)
  return(prec_list[[prec_drawn]])
}
#description: samples a graph from an order score list
#parameters:  order_score_list - list of parent scores compatible with order
sample_graph <- function(order_score_list){
  N <- length(order_score_list)
  G <- matrix(0, ncol=N, nrow=N)
  for(i in 1:N){
    lookup_table <- order_score_list[[i]]
    prob_vec <- as.numeric(exp(lookup_table-logSumExp(lookup_table)))
    par_matr_index <- sample(1:length(prob_vec), size=1, prob = prob_vec)
    parents <- names(lookup_table)[par_matr_index]
    index_set_parents <- as.numeric(unlist(strsplit(parents, ",")))
    G[i, index_set_parents] <- 1
  }
  return(G)
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