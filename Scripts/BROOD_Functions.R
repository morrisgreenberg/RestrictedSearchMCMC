library(Matrix)
library(gtools)
library(Rfast)
library(BiDAG)

# ==============================================================================
# 1. MAIN WRAPPER: graph_mcmc
# ==============================================================================
#' BROOD MCMC Sampler
#' @description Performs MCMC over both the order space and the search space (H).
#' @param H_0 Initial adjacency matrix representing the starting search space.
#' @param param Scoring parameters (BiDAG object).
#' @param d_expand number of graphs to sample at expansion steps
#' @param d_shrink number of graphs to sample at contraction steps
#' @param thresh threshold to add edges to the space from the d samples
#' @param c_star scaling constant for death rates (smaller means more inclusive search space)
#' @param space_move_prob Probability of attempting a search space expansion/contraction.
#' @param iter number of steps of the chain
#' @param warm_up amount of steps before adaptation starts (defaults to 10*N)
#' @param max_sparsity limit on the number of parents a node can have
#' @param blacklist list of any parent structures that are not allowed
#' @param move_type how to propose move types. Three current methods:
#                           a. relocate - always performs node relocation
#                           b. random - node relocation (NR) 1/3 of the time,
#                                       local transposition (LT) 1/3,
#                                       global swap (GS) 1/3
#                           c. Kuipers - GS 6/(N+7), LT t/(N+7), NR 1/(N+7)
#' @param sample_parameters if we also want to sample the covariance parameters for the Gaussian DAG Model (for DAG Wishart score)
#' @param plus1 If TRUE, uses the "plus1" look-ahead score for orders.
#' @param score_type Either "bge" score (Heckerman and Geiger, 1995), or "dag-wishart" score (Ben-David et al, 2015)
#' @param rounded If birth and death rates are rounded at 1, as is done in Mohammadi and Wit (2015)
#' @param max_change Maximum number of edges to add/remove in a single expansion/contraction
#' @param save_all_weights If TRUE, saves all birth-death weights, or just the subset for the thinned samples
#' @param sparse If TRUE, uses the Matrix::Matrix(sparse=TRUE) type matrices, or standard R matrices in saved output
#' @param verbose Outputs progress as the Markov Chain progresses
graph_mcmc <- function(H_0, param, d_expand=1, d_shrink=1, 
                       thresh=1e-9, c_star=1, space_move_prob=0.1,
                       iter=25000, warm_up=NULL, max_sparsity=18, blacklist=NULL, 
                       thinning=250, move_type="relocate", sample_parameters=FALSE,
                       plus1=FALSE, score_type="bge", rounded=FALSE, max_change=1,
                       save_all_weights=TRUE, sparse=TRUE, verbose=TRUE) {
  
  if(!(score_type %in% c("bge", "dag_wishart"))) stop("Please use a valid score function: 'bge' or 'dag_wishart' ")
  
  N <- nrow(H_0)
  max_change <- if(is.null(max_change)) 1 else max_change
  if(max_change < 1 | max_change > N) stop("Please make max_change between 1:N")
  
  warm_up <- if(is.null(warm_up)) 10 * N else warm_up
  total_iter <- iter + warm_up + 1
  
  # -- Initialization --
  prec_b <- 1:N
  K_b <- round(sqrt(N/2))
  H_b <- H_0
  
  # Initial Mappings
  mappings <- parents_mapping(H_0)
  banned_mappings <- banned_parents_mapping(mappings, prec_b, TRUE, TRUE)
  plus_mappings <- plus_parents_mapping(H_0, 1, mappings, blacklist)
  
  # Pre-calculate Scores
  score_object <- score_plus_space(H_0, mappings, plus_mappings, param, score_type, N)
  banned_object <- create_banned_plus_parent_table(H_0, mappings, plus_mappings, score_object$full_list)
  
  full_scores <- score_object$curr_scores
  banned_scores <- banned_object$curr_scores
  full_plus_scores <- score_object$full_list
  banned_plus_scores <- banned_object$plus_scores
  
  order_score <- calculate_order_score(banned_mappings, banned_scores, plus1, 
                                       if(plus1) banned_plus_scores else NULL)
  
  # -- Storage Allocation --
  thinned_samples <- seq(warm_up + 1, total_iter, thinning)
  B_saved <- length(thinned_samples)
  
  # Initialize output containers
  Gs <- if(sparse) vector(mode="list", B_saved) else array(NA, c(N, N, B_saved))
  G1s <- if(!plus1) (if(sparse) vector(mode="list", B_saved) else array(NA, c(N, N, B_saved))) else NULL
  Hs <- if(sparse) vector(mode="list", B_saved) else array(NA, c(N, N, B_saved))
  skels <- if(sparse) vector(mode="list", B_saved) else array(NA, c(N, N, B_saved))
  precs <- matrix(nrow=B_saved, ncol=N)
  Ks <- numeric(B_saved)
  weight_vec <- numeric(B_saved)
  update_vec <- character(B_saved)
  Ls <- if(sample_parameters) (if(sparse) vector(mode="list", B_saved) else array(dim=c(N, N, B_saved))) else NULL
  Ds <- if(sample_parameters) matrix(nrow=B_saved, ncol=N) else NULL
  
  # -- MCMC Loop --
  for(b in 1:total_iter) {
    if(verbose && b %% 500 == 0) message(sprintf("Iteration: %d", b))
    
    save_curr_weight <- save_all_weights || (b %in% thinned_samples)
    
    step <- mcmc_sampler_step(prec_b, H_b, K_b, b, d_expand, d_shrink, 
                              banned_scores, mappings, banned_mappings, plus_mappings, 
                              full_scores, banned_plus_scores, full_plus_scores, 
                              order_score, param, c_star, space_move_prob, 
                              score_type, sample_parameters, rounded, max_change, 
                              thresh, move_type, warm_up, max_sparsity, plus1, 
                              blacklist, save_curr_weight, verbose=FALSE)
    
    # Update state
    prec_b <- step$prec_t_plus1
    K_b <- step$K_t_plus1
    banned_scores <- step$banned_scores
    full_scores <- step$order_scores
    mappings <- step$par_mappings
    banned_plus_scores <- step$banned_plus_scores
    full_plus_scores <- step$order_plus_scores
    order_score <- step$curr_order_score
    plus_mappings <- step$plus_par_mappings
    banned_mappings <- step$banned_par_mappings
    D_b <- if(sample_parameters) step$D_t_plus1 else NULL
    curr_weight <- if(save_curr_weight) step$weight else NULL
    curr_update <- step$update_type
    
    H_b <- if(sparse) Matrix::Matrix(step$H_t_plus1, sparse=TRUE) else step$H_t_plus1
    G_b <- if(sparse) Matrix::Matrix(step$G_t_plus1, sparse=TRUE) else step$G_t_plus1
    G1_b <- if(!plus1) (if(sparse) Matrix::Matrix(step$G1_t_plus1, sparse=TRUE) else step$G1_t_plus1) else NULL
    L_b <- if(sample_parameters) (if(sparse) Matrix::Matrix(step$L_t_plus1, sparse=TRUE) else step$L_t_plus1) else NULL
    
    
    # Store thinned samples
    if(b %in% thinned_samples) {
      idx <- which(thinned_samples == b)
      precs[idx,] <- prec_b
      Ks[idx] <- K_b
      weight_vec[idx] <- curr_weight
      update_vec[idx] <- curr_update
      if(sample_parameters) Ds[idx,] <- D_b
      
      if(!sparse){
        Gs[,,idx] <- G_b
        skels[,,idx] <- (G_b + t(G_b)>0)*1
        Hs[,,idx] <- H_b
        if(!plus1) G1s[,,idx] <- G1_b
        if(sample_parameters) Ls[,,idx] <- L_b
      }
      else{
        Gs[[idx]] <- G_b
        skels[[idx]] <- (G_b + t(G_b)>0)*1
        Hs[[idx]] <- H_b
        if(!plus1) G1s[[idx]] <- G1_b
        if(sample_parameters) Ls[[idx]] <- L_b
      }
      
    }
  }
  
  out_list <- list(orders=precs, graphs=Gs, skeletons=skels, spaces=Hs,
                   sparsity=Ks, weights=weight_vec, update=update_vec)
  if(!plus1) out_list$plusgraphs <- G1s
  if(sample_parameters) {
    out_list$D_vecs <- Ds
    out_list$L_matrs <- Ls
  }
  return(out_list)
}

# ==============================================================================
# 2. INTERNAL FUNCTIONS
# ==============================================================================

#description: 1 step of the MCMC sampler
#' @param prec_t order at step t
#' @param H_t search space at step t
#' @param K_t sparsity at step t (currently on hold)
#' @param t step number in the chain
#' @param d_expand number of graphs to draw at expansion
#' @param d_shrink number of graphs to draw at shrink
#' @param space_banned_score_list banned score list for scoring orders
#' @param map_pars hash tables for order scoring
#' @param full_score_list score lists for every valid parent set
#' @param plus_banned_list (+1) banned score list for scoring orders
#' @param plus_score_list (+1) score lists for every valid parent set
#' @param order_score nodewise score values for the current order
#' @param param score parameter object, constructed from package BiDAG
#' @param c_star scaling constant to offset death rates if wanting a larger search space
#' @param prob_adapt probability of transitioning to another space
#' @param thresh threshold to add edges to the space from the d samples
#' @param move_probs how to propose move types. Three current methods:
#                          a. relocate - always performs node relocation
#                          b. random - node relocation (NR) 1/3 of the time,
#                                     local transposition (LT) 1/3,
#                                     global swap (GS) 1/3
#                         c. Kuipers - NR 6/(t+7), LT t/(t+7), GS 1/(t+7)
#' @param warm_up amount of steps before adaptation starts (defaults to 10*N)
#' @param max_sparsity limit on the number of parents a node can have
#' @param plus1 indicates whether plus score or standard score is used
#' @param blacklist list of any parent structures that are not allowed
#' @param save_waiting_times whether to save the waiting times for each step
#' @param verbose prints the step in the chain number if TRUE, and when shrinking/expanding occurs
mcmc_sampler_step <- function(prec_t, H_t, K_t, t, d_expand, 
                              d_shrink, space_banned_score_list, 
                              map_pars, banned_pars, plus_pars,
                              full_score_list, plus_banned_list, 
                              plus_score_list, order_score,
                              param, c_star, prob_adapt, score_type,
                              to_sample_params, rounded, max_change, thresh, 
                              move_probs, warm_up, max_sparsity, plus1, 
                              blacklist, save_waiting_times, verbose){
  N <- nrow(H_t)
  l <- 1
  
  is_adaption <- F
  if(t>warm_up){
    is_adaption <- sample(c(T, F), 1, prob=c(prob_adapt, 1-prob_adapt))
  }
  is_contraction <- F
  is_expansion <- F
  
  if(is_adaption | save_waiting_times){
    birth_rates <- calculate_log_birth_rate(H_t, plus_banned_list, max_sparsity,
                                            TRUE, prec_t, plus_pars, banned_pars, rounded)
    death_rates<- calculate_log_death_rate(H_t, full_score_list, map_pars,
                                           space_banned_score_list, TRUE, prec_t, banned_pars,
                                           c_star, rounded)
    logsum_br <- logSumExp(birth_rates)
    logsum_dr <- logSumExp(death_rates)
    w_t <- -logSumExp(c(logsum_br, logsum_dr))
  }
  
  type_of_update <- "displacement"
  
  if(is_adaption){
    
    #sampling whether we do any expansion/contraction steps
    process <- sample(c("Birth", "Death"), size=1, 
                      prob=exp(c(logsum_br+w_t, logsum_dr+w_t)))
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
                        prob=c(6/(N+7), N/(N+7), 1/(N+7)))
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
      type_of_update <- "birth"
      if(max_change == 1){
        edge <- sample(1:N^2, 1, prob=exp(as.numeric(birth_rates-logsum_br)))
      }
      else{
        num_nodes <- sample(1:max_change, size=1)
        node_birth_rates_all <- rowLogSumExps(birth_rates)
        log_prob_vec <- as.numeric(node_birth_rates_all-logsum_br)
        if(length(which(exp(log_prob_vec)>0))<num_nodes){
          num_nodes <- length(which(exp(log_prob_vec)>0))
        }
        if(num_nodes == 0){
          edge <- integer(0)
        }
        else{
          row_choices <- sample(1:N, size=num_nodes, 
                                prob=exp(log_prob_vec))
          col_choices <- sapply(row_choices, function(i){
            sample(1:N, size=1, prob=exp(as.numeric(birth_rates[i,]-node_birth_rates_all[i])))
          })
          edge <- row_choices + (col_choices-1)*N
        }
      }
      
      if(verbose){print("pre-expand")}
      if(length(edge)>0){
        space_output <- expand_search_space(H_t, edge, max_sparsity)
        H_proposal <- space_output$H_new
        update_nodes <- space_output$updatenodes
      }
      else{
        H_proposal <- H_t
        update_nodes <- integer(0)
      }
      if(verbose){print("post-expand")}
    }
    else{
      type_of_update <- "death"
      if(max_change==1){
        edge <- sample(1:N^2, 1, prob=exp(as.numeric(death_rates-logsum_dr)))
      }
      else{
        num_nodes <- sample(1:max_change, size=1)
        node_death_rates_all <- rowLogSumExps(death_rates)
        log_prob_vec <- as.numeric(node_death_rates_all-logsum_dr)
        if(length(which(exp(log_prob_vec)>0))<num_nodes){
          num_nodes <- length(which(exp(log_prob_vec)>0))
        }
        if(num_nodes == 0){
          edge <- NULL
        }
        else{
          row_choices <- sample(1:N, size=num_nodes, 
                                prob=exp(log_prob_vec))
          col_choices <- sapply(row_choices, function(i){
            sample(1:N, size=1, prob=exp(as.numeric(death_rates[i,]-node_death_rates_all[i])))
          })
          edge <- row_choices + (col_choices-1)*N
        }
      }
      
      if(verbose){print("pre-shrink")}
      if(length(edge)>0){
        space_output <- shrink_search_space(H_t, edge)
        H_proposal <- space_output$H_new
        update_nodes <- space_output$updatenodes
      }
      else{
        H_proposal <- H_t
        update_nodes <- integer(0)
      }
      
      if(verbose){print("post-shrink")}
    }
    
    if(length(update_nodes)>0){
      proposed_mappings <- parents_mapping(H_proposal, N, update_nodes,TRUE, map_pars)
      proposed_plus_mappings <- plus_parents_mapping(H_proposal, l, proposed_mappings, blacklist)
      proposed_banned_mappings <- banned_parents_mapping(proposed_mappings, prec_t, TRUE, TRUE)
      if(is_expansion){
        score_object <- score_plus_space(H_proposal, proposed_mappings, proposed_plus_mappings,
                                         param, score_type, N, update_nodes, 
                                         has_scores_orig=TRUE, has_plus_orig=TRUE, 
                                         full_score_list, plus_score_list, map_pars, 
                                         plus_pars, H_t, is_shrink=FALSE)
        
        banned_object <- create_banned_plus_parent_table(H_proposal, proposed_mappings,
                                                         proposed_plus_mappings, 
                                                         score_object$full_list,
                                                         N, update_nodes, 
                                                         has_scores_orig=TRUE, 
                                                         space_banned_score_list,
                                                         plus_banned_list, map_pars,
                                                         H_t, is_shrink=FALSE)
      }
      else{
        score_object <- score_plus_space(H_proposal, proposed_mappings, proposed_plus_mappings,
                                         param, score_type, N, update_nodes, 
                                         has_scores_orig=TRUE, has_plus_orig=TRUE, 
                                         full_score_list, plus_score_list, map_pars, 
                                         plus_pars, H_t, is_shrink=TRUE)
        banned_object <- create_banned_plus_parent_table(H_proposal, proposed_mappings,
                                                         proposed_plus_mappings, 
                                                         score_object$full_list,
                                                         N, update_nodes, 
                                                         has_scores_orig=TRUE,
                                                         space_banned_score_list,
                                                         plus_banned_list, map_pars,
                                                         H_t, is_shrink=TRUE)
      }
      proposed_full_scores <- score_object$curr_scores
      proposed_plus_scores <- score_object$full_list
      
      proposed_banned_scores <- banned_object$curr_scores
      proposed_banned_plus_scores <- banned_object$plus_scores
      
      # M-H ratios are the ratios of waiting times of new vs. old
      birth_rates_proposed <- calculate_log_birth_rate(H_proposal, 
                                                       proposed_banned_plus_scores,
                                                       max_sparsity, TRUE, prec_t,
                                                       proposed_plus_mappings,
                                                       proposed_banned_mappings)
      death_rates_proposed <- calculate_log_death_rate(H_proposal, proposed_full_scores, 
                                                       proposed_mappings, proposed_banned_scores,
                                                       TRUE, prec_t, proposed_banned_mappings,
                                                       c_star)
      
      w_proposed <- -logSumExp(c(logSumExp(birth_rates_proposed),logSumExp(death_rates_proposed)))
      r_t <- w_proposed-w_t
      acpt_rate_t <- min(1, exp(r_t))
      
      # perform M-H accept/reject step
      u_t <- runif(1)
      to_update <- u_t <= acpt_rate_t
      if(verbose){
        print(paste("acceptance rate", acpt_rate_t, sep=": "))
        print(paste("w_proposed", w_proposed, sep=": "))
        print(paste("w_t", w_t, sep=": "))
        print(paste("u_t", u_t, sep=": "))
        print(paste("update?", to_update, sep=": "))
      }
    }
    else{
      to_update <- FALSE
    }
    
    if(to_update){
      new_mappings <- proposed_mappings
      new_plus_mappings <- proposed_plus_mappings
      new_banned_mappings <- proposed_banned_mappings
      new_full_scores <- proposed_full_scores
      new_plus_scores <- proposed_plus_scores
      new_banned_scores <- proposed_banned_scores
      new_banned_plus_scores <- proposed_banned_plus_scores
      if(plus1){
        new_order_score <- calculate_order_score(new_banned_mappings, new_banned_scores, TRUE,
                                                 new_banned_plus_scores)
      }
      else{
        new_order_score <- calculate_order_score(new_banned_mappings, new_banned_scores)
      }
      H_t_plus1 <- H_proposal
      w_t_plus1 <- w_proposed
    }
    else{
      new_mappings <- map_pars
      new_plus_mappings <- plus_pars
      new_banned_mappings <- banned_pars
      new_full_scores <- full_score_list
      new_plus_scores <- plus_score_list
      new_banned_scores <- space_banned_score_list
      new_banned_plus_scores <- plus_banned_list
      new_order_score <- order_score
      H_t_plus1 <- H_t
      w_t_plus1 <- w_t
      type_of_update <- "rejection"
    }
    #sample new order
    if(!plus1){
      update_order_obj <- implement_order(prec_t, move_type,
                                          new_banned_scores, new_mappings,
                                          new_banned_mappings, new_order_score, H=H_t_plus1)
      prec_prime <- update_order_obj$order
      order_score_prime <- update_order_obj$score 
      if(move_type != "node relocation"){
        prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                             new_banned_scores, new_mappings)
        if(sum(prec_t_plus1!=prec_t)>0){
          order_score <- order_score_prime
        }
      }
      else{
        prec_t_plus1 <- prec_prime
        order_score <- order_score_prime
      }
    }
    else{
      update_order_obj <- implement_order(prec_t, move_type,
                                          new_banned_scores, new_mappings,
                                          new_banned_mappings, new_order_score,
                                          TRUE, new_banned_plus_scores, H_t_plus1)
      prec_prime <- update_order_obj$order
      order_score_prime <- update_order_obj$score
      if(move_type != "node relocation"){
        prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                             new_banned_scores, new_mappings, TRUE, 
                                             new_banned_plus_scores)
        if(sum(prec_t_plus1!=prec_t)>0){
          order_score <- order_score_prime
        }
      }
      else{
        prec_t_plus1 <- prec_prime
        order_score <- order_score_prime
      }
      
    }
    banned_pars <- new_banned_mappings
    if(plus1){
      graph_t_plus1 <- sample_graph(new_full_scores, prec_t_plus1, new_mappings, 
                                    banned_pars, plus_1 = TRUE, new_plus_scores,
                                    new_banned_plus_scores, new_plus_mappings)
    }
    else{
      graph_t_plus1 <- sample_graph(new_full_scores, prec_t_plus1, new_mappings, 
                                    banned_pars, plus_1 = FALSE)
      graph_plus1_t_plus1 <- sample_graph(new_full_scores, prec_t_plus1, new_mappings, 
                                          banned_pars, plus_1 = TRUE, new_plus_scores,
                                          new_banned_plus_scores, new_plus_mappings)
    }
    
  }
  
  else{
    if(!plus1){
      #standard sampling of an order
      update_order_obj <- implement_order(prec_t, move_type, 
                                          space_banned_score_list, map_pars,
                                          banned_pars, order_score, H=H_t)
      prec_prime <- update_order_obj$order
      order_score_prime <- update_order_obj$score
      if(move_type != "node relocation"){
        prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                             space_banned_score_list, map_pars)
        if(sum(prec_t_plus1!=prec_t)>0){
          order_score <- order_score_prime
        }
      }
      else{
        prec_t_plus1 <- prec_prime
        order_score <- order_score_prime
      }
      
    }
    else{
      #standard sampling of an order
      update_order_obj <- implement_order(prec_t, move_type, 
                                          space_banned_score_list, map_pars,
                                          banned_pars, order_score,
                                          TRUE, plus_banned_list, H_t)
      prec_prime <- update_order_obj$order
      order_score_prime <- update_order_obj$score
      if(move_type != "node relocation"){
        prec_t_plus1 <- sample_from_2_orders(prec_t, prec_prime, 
                                             space_banned_score_list, map_pars,
                                             TRUE, plus_banned_list)
        if(sum(prec_t_plus1!=prec_t)>0){
          order_score <- order_score_prime
        }
      }
      else{
        prec_t_plus1 <- prec_prime
        order_score <- order_score_prime
      }
    }
    if(sum(prec_t_plus1 != prec_t)> 0){
      
      if(move_type != "node relocation"){
        banned_pars <- update_order_obj$banned_pars
      } else {
        banned_pars <- banned_parents_mapping(map_pars, prec_t_plus1, TRUE, TRUE)
      }
    }
    
    if(plus1){
      graph_t_plus1 <- sample_graph(full_score_list, prec_t_plus1, map_pars, banned_pars,
                                    plus_1=TRUE, plus_score_list, plus_banned_list,
                                    plus_pars)
    }
    else{
      graph_t_plus1 <- sample_graph(full_score_list, prec_t_plus1, map_pars, 
                                    banned_pars, plus_1 = FALSE)
      graph_plus1_t_plus1 <- sample_graph(full_score_list, prec_t_plus1, map_pars, banned_pars,
                                          plus_1=TRUE, plus_score_list, plus_banned_list,
                                          plus_pars)
    }
    
    
    #standard update of the search space
    H_t_plus1 <- H_t
    new_full_scores <- full_score_list
    new_banned_scores <- space_banned_score_list
    new_mappings <- map_pars
    new_plus_mappings <- plus_pars
    new_plus_scores <- plus_score_list
    new_banned_plus_scores <- plus_banned_list
    if(save_waiting_times){
      w_t_plus1 <- w_t
    }
  }
  
  K_t_plus1 <- K_t
  
  if(to_sample_params){
    pars_t_plus1 <- sample_DL_parameters(graph_t_plus1, prec_t_plus1, param)
    L_matr <- pars_t_plus1$L
    D_vec <- pars_t_plus1$D
  }
  
  
  output_list <- list(prec_t_plus1=prec_t_plus1, G_t_plus1=graph_t_plus1,
                      H_t_plus1=H_t_plus1, K_t_plus1=K_t_plus1,
                      order_scores = new_full_scores, 
                      banned_scores = new_banned_scores,
                      par_mappings = new_mappings,
                      order_plus_scores = new_plus_scores,
                      banned_plus_scores = new_banned_plus_scores,
                      curr_order_score = order_score,
                      banned_par_mappings = banned_pars,
                      plus_par_mappings = new_plus_mappings,
                      update_type = type_of_update)
  
  if(save_waiting_times){
    output_list$weight <- w_t_plus1
  }
  if(to_sample_params){
    output_list$L_t_plus1 <- L_matr
    output_list$D_t_plus1 <- D_vec
  }
  if(!plus1){
    output_list$G1_t_plus1 <- graph_plus1_t_plus1
  }
  
  return(output_list)
  
}

score_plus_space <- function(H, map_pars, plus_pars, param,
                             score_type, N=ncol(H), updatenodes=1:N,
                             has_scores_orig=FALSE, H_scores=NULL, 
                             has_plus_orig=FALSE, H_plus_scores=NULL,
                             map_pars_orig=NULL, plus_pars_orig=NULL,
                             H_orig=NULL, is_shrink=FALSE){
  score_list <- vector(mode="list", length=N)
  space_scores <- vector(mode="list", length=N)
  for(i in 1:N){
    if(!(i %in% updatenodes) & has_plus_orig){
      space_scores[[i]] <- H_scores[[i]]
      score_list[[i]] <- H_plus_scores[[i]]
    }
    else{
      if(is_shrink & has_plus_orig){
        combos_orig <- map_pars_orig$par_pset[[i]]
        removed_pars <- which(H_orig[i,]-H[i,]==1)
        N_removed <- length(removed_pars)
        if(N_removed==1){
          removed_idx <- integer(0)
          for(j in 1:ncol(combos_orig)){
            removed_idx <- c(removed_idx, which(combos_orig[,j] %in% removed_pars))
          }
          if(length(removed_idx)>0){
            plus_scores <- H_plus_scores[[i]]
            score_matr <- matrix(nrow=nrow(plus_scores)-length(removed_idx),
                                 ncol=ncol(plus_scores)+length(removed_pars))
            new_cols <- which(plus_pars$par_pset[[i]][,1] %in% removed_pars)
            score_matr[,-new_cols] <- plus_scores[-removed_idx,]
            score_matr[,new_cols] <- plus_scores[sort(removed_idx),1]
            space_scores[[i]] <- matrix(score_matr[,1], ncol=1)
            score_list[[i]] <- score_matr
          }
          else{
            space_scores[[i]] <- H_scores[[i]]
            score_list[[i]] <- H_plus_scores[[i]]
          }
        }
        else if(N_removed>1){
          removed_idx <- vector(mode="list", length=N_removed)
          removed_idx_all <- integer(0)
          for(k in 1:N_removed){
            removed_idx_vec <- integer(0)
            for(j in 1:ncol(combos_orig)){
              removed_idx_vec <- c(removed_idx_vec, which(combos_orig[,j] %in% removed_pars[k]))
            }
            removed_idx[[k]] <- removed_idx_vec
            removed_idx_all <- c(removed_idx_all, removed_idx_vec)
          }
          removed_idx_unique <- unique(removed_idx_all)
          removed_idx_multiple <- as.numeric(names(which(table(removed_idx_all)>1)))
          
          plus_scores <- H_plus_scores[[i]]
          score_matr <- matrix(nrow=nrow(plus_scores)-length(removed_idx_unique),
                               ncol=ncol(plus_scores)+length(removed_pars))
          new_cols <- which(plus_pars$par_pset[[i]][,1] %in% removed_pars)
          score_matr[,-new_cols] <- plus_scores[-removed_idx_unique,]
          for(k in 1:N_removed){
            removed_idx_k <- setdiff(removed_idx[[k]], removed_idx_multiple)
            score_matr[,new_cols[k]] <-plus_scores[sort(removed_idx_k),1]
          }
          space_scores[[i]] <- matrix(score_matr[,1], ncol=1)
          score_list[[i]] <- score_matr
        }
        else{
          space_scores[[i]] <- H_scores[[i]]
          score_list[[i]] <- H_plus_scores[[i]]
        }
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
              if(score_type == "bge"){
                score_matr[j,1] <- bge_score_node(i, parent_group, N, param)
              }
              else if(score_type == "dag_wishart"){
                score_matr[j,1] <- dagwishart_score_node(i, parent_group, N, param)
              }
              else{
                stop("Please use a valid score type: 'bge' or 'dag_wishart' ")
              }
              
            }
          }
          score_list[[i]] <- score_matr
          space_scores[[i]] <- matrix(score_matr[,1], ncol=1)
        }
        else{
          score_matr <- matrix(nrow=n_parent_sets, ncol=n_nonparent_sets)
          for(k in 1:n_parent_sets){
            if(k==1){
              if(score_type == "bge"){
                score_matr[k,] <- bge_score_plus_parent(i, NULL, plus_combos[-1,1], N, param)
              }
              else if(score_type == "dag_wishart"){
                score_matr[k,] <- dagwishart_score_plus_parent(i, NULL, plus_combos[-1,1], N, 
                                                               param)
              }
              else{
                stop("Please use a valid score type: 'bge' or 'dag_wishart' ")
              }
            }
            else{
              if(score_type == "bge"){
                score_matr[k,] <- bge_score_plus_parent(i, combos[k,1:c(par_vec[k])], 
                                                        plus_combos[-1,1], N, param)
              }
              else if(score_type == "dag_wishart"){
                score_matr[k,] <- dagwishart_score_plus_parent(i, combos[k,1:c(par_vec[k])], 
                                                               plus_combos[-1,1], N, param)
              }
              else{
                stop("Please use a valid score type: 'bge' or 'dag_wishart' ")
              }
            }
            
          }
          space_scores[[i]] <- matrix(score_matr[,1], ncol=1)
          score_list[[i]] <- score_matr
        }
      }
    }
  }
  return(list("curr_scores"=space_scores, "full_list"=score_list))
}

logMinusExp <-function(lx, ly){
  value <- lx-ly
  if(min(value)<=0){
    min_idx <- which.min(value)
    stop(paste("Computing log of a negative number. First number less than second at index",
               min_idx, sep=" "))
  }
  value_return <- lx
  idx_small <- (ly != -Inf) & (value <= log(2))
  idx_large <- (ly != -Inf) & (value > log(2))
  
  value_return[idx_small] <- lx[idx_small] + log(-expm1(-value[idx_small]))
  value_return[idx_large] <- lx[idx_large] + log1p(-exp(-value[idx_large]))
  return(value_return)
}

create_banned_plus_parent_table <- function(H, map_pars, plus_pars,
                                            score_plus_list,
                                            N=ncol(H), updatenodes=1:N,
                                            has_scores_orig=FALSE, orig_scores=NULL,
                                            orig_scores_plus=NULL, orig_map_pars=NULL,
                                            H_orig=NULL, is_shrink=FALSE){
  orderscore_plus<-vector(mode="list", length=N)
  orderscore_curr<-vector(mode="list", length=N)
  for(i in 1:N){
    if(!(i %in% updatenodes) & has_scores_orig){
      orderscore_plus[[i]] <- orig_scores_plus[[i]]
      orderscore_curr[[i]] <- orig_scores[[i]]
    }
    else{
      if(is_shrink & has_scores_orig){
        combos_orig <- orig_map_pars$par_pset[[i]]
        N_pars_orig <- length(orig_map_pars$par_names[[i]])
        removed_pars <- which(H_orig[i,]-H[i,]==1)
        N_removed <- length(removed_pars)
        if(N_removed==1){
          new_cols <- which(plus_pars$par_pset[[i]][,1] %in% removed_pars)
          keep_idx <- integer(0)
          if(N_pars_orig > 0){
            for(colsy in 1:N_pars_orig){
              keep_idx <- c(keep_idx, which(combos_orig[,colsy] %in% removed_pars))
            }
            if(length(keep_idx)>0){
              plus_scores <- orig_scores_plus[[i]]
              curr_scores <- orig_scores[[i]]
              banned_matr <- matrix(nrow=length(keep_idx),
                                    ncol=ncol(plus_scores)+N_removed)
              banned_matr[,-new_cols] <- plus_scores[sort(keep_idx),]
              banned_matr[,new_cols] <- logMinusExp(curr_scores[-keep_idx,], 
                                                    banned_matr[,1])
              
              orderscore_plus[[i]] <- banned_matr
              orderscore_curr[[i]] <- matrix(banned_matr[,1], ncol=1)
            }
            else{
              orderscore_plus[[i]] <- orig_scores_plus[[i]]
              orderscore_curr[[i]] <- orig_scores[[i]]
            }
          }
          else{
            orderscore_plus[[i]] <- orig_scores_plus[[i]]
            orderscore_curr[[i]] <- orig_scores[[i]] 
          }
        }
        else if(N_removed > 1){
          new_cols <- which(plus_pars$par_pset[[i]][,1] %in% removed_pars)
          keep_idx_list <- vector(mode="list", length=N_removed)
          keep_idx_all <- integer(0)
          if(N_pars_orig > 0){
            for(k in 1:N_removed){
              keep_idx <- integer(0)
              for(colsy in 1:N_pars_orig){
                keep_idx <- c(keep_idx, which(combos_orig[,colsy]==removed_pars[[k]]))
              }
              keep_idx_list[[k]] <- keep_idx
              keep_idx_all <- c(keep_idx_all, keep_idx)
            }
            keep_idx_every <- as.numeric(names(which(table(keep_idx_all)==N_removed)))
            keep_idx_off1 <- as.numeric(names(which(table(keep_idx_all)==(N_removed-1))))
          }
          else{
            orderscore_plus[[i]] <- orig_scores_plus[[i]]
            orderscore_curr[[i]] <- orig_scores[[i]]
          }
          if(length(keep_idx_every) > 0){
            plus_scores <- orig_scores_plus[[i]]
            curr_scores <- orig_scores[[i]]
            banned_matr <- matrix(nrow=length(keep_idx_every),
                                  ncol=ncol(plus_scores)+N_removed)
            banned_matr[,-new_cols] <- plus_scores[sort(keep_idx_every),]
            for(k in 1:N_removed){
              removed_idx <- setdiff(keep_idx_off1, keep_idx_list[[k]])
              banned_matr[,new_cols[k]] <- logMinusExp(curr_scores[sort(removed_idx),],
                                                       banned_matr[,1])
            }
            orderscore_plus[[i]] <- banned_matr
            orderscore_curr[[i]] <- matrix(banned_matr[,1], ncol=1)
          }
          else{
            orderscore_plus[[i]] <- orig_scores_plus[[i]]
            orderscore_curr[[i]] <- orig_scores[[i]]
          }
        }
        else{
          orderscore_plus[[i]] <- orig_scores_plus[[i]]
          orderscore_curr[[i]] <- orig_scores[[i]]
        }
        
      }
      else{
        N_outside <- length(plus_pars$par_pset[[i]][,1])
        N_pars <- length(map_pars$par_names[[i]])
        zeta_matr <- score_plus_list[[i]]
        N_scores <- nrow(zeta_matr)
        if(N_pars > 0){
          for(t in 1:N_pars){
            index_t <- numeric(0)
            for(colsy in 1:N_pars){
              index_t <- c(index_t, which(t==map_pars$idx_pset[[i]][,colsy]))
            }
            for(rowsy in index_t){
              mapped_val <- map_pars$maps[[i]]$forward[rowsy]-2^(t-1)
              min_rowsy <- map_pars$maps[[i]]$backwards[mapped_val]
              if(ncol(zeta_matr)>1){
                zeta_matr[rowsy,] <- colLogSumExps(zeta_matr[c(rowsy, min_rowsy),])
              }
              else{
                zeta_matr[rowsy,] <- logSumExp(zeta_matr[c(rowsy, min_rowsy),])
              }
            }
          }
          orderscore_plus[[i]] <- zeta_matr[N_scores:1,]
          orderscore_curr[[i]] <- matrix(zeta_matr[N_scores:1,1], ncol=1)
        }
        else{
          orderscore_plus[[i]] <- zeta_matr
          orderscore_curr[[i]] <- matrix(zeta_matr[,1], ncol=1)
        }
      }
    }
    
  }
  return(list("plus_scores"=orderscore_plus, "curr_scores"=orderscore_curr))
}

#description: wrapper function for hash tables 
#' @param H           current space
#' @param N           number of nodes
#' @param updatenodes nodes where scores need to be updated
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
#' @param H           search space
#' @param updatenodes nodes where scores need to be updated
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
        total_rows <- 2^N_parents
        parent_matrix <- matrix(NA, nrow=total_rows, ncol=N_parents)
        
        current_row <- 2
        for(r in 1:N_parents){
          new_rows <- gtools::combinations(N_parents, r)
          n_combos <- nrow(new_rows)
          
          parent_matrix[current_row:(current_row + n_combos - 1), 1:r] <- new_rows
          current_row <- current_row + n_combos
        }
        par_list[[i]] <- parent_matrix
      }
    }
  }
  return(par_list)
}
#description: creates hash table by node names
#' @param H search space
#' @param updatenodes nodes where scores need to be updated
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
        total_rows <- 2^N_parents
        parent_matrix <- matrix(NA, nrow=total_rows, ncol=N_parents)
        
        current_row <- 2
        for(r in 1:N_parents){
          new_rows <- gtools::combinations(N_parents, r, parent_nodes)
          n_combos <- nrow(new_rows)
          
          parent_matrix[current_row:(current_row + n_combos - 1), 1:r] <- new_rows
          current_row <- current_row + n_combos
        }
        par_list[[i]] <- parent_matrix
      }
    }
  }
  return(par_list)
}

# these two functions are directly adapting dettwobytwo, DAGcorescore
#' @author 	Polina Suter [aut, cre], Jack Kuipers [aut] (originally)
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
           D<-TN[parentnodes,parentnodes]
           detD<-dettwobytwo(D)
           logdetD<-log(detD)
           B<-TN[j,parentnodes]
           logdetpart2<-log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
           }
         },
         
         {# otherwise we use cholesky decomposition to perform both
           D<-as.matrix(TN[parentnodes,parentnodes])
           choltemp<-chol(D)
           logdetD<-2*sum(log(choltemp[(lp+1)*c(0:(lp-1))+1]))
           B<-TN[j,parentnodes]
           logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
           }
         })
  return(corescore)
}

#an efficient vectorized "plus-1" version of the BGe scoring function
bge_score_plus_parent <- function(j, parentnodes, plus_parentnodes, N, param){
  TN<-param$TN
  awpN<-param$awpN
  scoreconstvec<-param$scoreconstvec
  
  lp<-length(parentnodes) #number of parents
  lpp <- length(plus_parentnodes)
  corescore_vec <- numeric(lpp+1)
  awpNd2<-(awpN-N+lp+1)/2
  A<-TN[j,j]
  
  TN_noplus <- TN[parentnodes, parentnodes]
  TN_j_noplus <- TN[j, parentnodes]
  TN_j_plus <- TN[j, plus_parentnodes]
  TN_plus_off <- TN[parentnodes, plus_parentnodes]
  TN_plus_on <- TN[plus_parentnodes, plus_parentnodes]
  TN_plus_diag <- TN_plus_on[(lpp+1)*c(0:(lpp-1))+1]
  
  switch(as.character(lp),
         "0"={# just a single term if no parents
           corescore_vec[1] <- scoreconstvec[lp+1] -awpNd2*log(A)
           D<-TN_plus_diag
           logdetD<-log(D)
           B<-TN_j_plus
           logdetpart2<-log(A-B^2/D)
           corescores <- scoreconstvec[lp+2]-(awpNd2+1/2)*logdetpart2 - logdetD/2
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             corescores <- corescores - as.numeric(param$logedgepmat[parentnodes, j])
           }
           corescore_vec[2:(lpp+1)] <- corescores
         },
         
         "1"={# no need for matrices
           D<-TN_noplus
           logdetD<-log(D)
           B<-TN_j_noplus
           logdetpart2<-log(A-B^2/D)
           
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           
           logdetD_plus <- log(D*TN_plus_diag - as.numeric(TN_plus_off^2))
           logdetpart2_plus <- log((D-B^2/A)*(TN_plus_diag-TN_j_plus^2/A)-
                                     (TN_plus_off-TN_j_plus*B/A)^2)+log(A)-logdetD_plus
           corescores <- scoreconstvec[lp+2]-(awpNd2+1/2)*logdetpart2_plus - logdetD_plus/2
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalization
             penalty <- param$logedgepmat[parentnodes, j]
             plus_penalties <- as.numeric(param$logedgepmat[plus_parentnodes, j])
             corescore <- corescore - penalty
             corescores <- corescores - penalty - plus_penalties
           }
           corescore_vec[1] <- corescore
           corescore_vec[2:(lpp+1)] <- corescores
         },
         
         {# otherwise we use cholesky decomposition to perform both
           D <- TN_noplus
           choltemp<-chol(D)
           logdetD<-2*sum(log(choltemp[(lp+1)*c(0:(lp-1))+1]))
           B<-TN_j_noplus
           c_noplus <- backsolve(choltemp,B,transpose=TRUE)
           logdetpart2<-log(A-sum(c_noplus^2))
           
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           
           choltemp_new_12 <- backsolve(choltemp, TN_plus_off, transpose=TRUE)
           
           if(lpp > 1){
             choltemp_new_22 <- sqrt(TN_plus_diag - colsums(choltemp_new_12^2))
           }
           else{
             choltemp_new_22 <- sqrt(TN_plus_diag - sum(choltemp_new_12^2))
           }
           c_plusses <- as.numeric((TN_j_plus-t(choltemp_new_12) %*% c_noplus) / choltemp_new_22)
           logdetD_plus <- logdetD + 2*log(choltemp_new_22)
           logdetpart2_plus <- log(A-sum(c_noplus^2)-c_plusses^2)
           corescores <- scoreconstvec[lp+2]-(awpNd2+1/2)*logdetpart2_plus - logdetD_plus/2
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             penalty <- sum(param$logedgepmat[parentnodes, j])
             plus_penalties <- as.numeric(param$logedgepmat[plus_parentnodes, j])
             corescore <- corescore - penalty
             corescores <- corescores - penalty - plus_penalties
           }
           corescore_vec[1] <- corescore
           corescore_vec[2:(lpp+1)] <- corescores
           
         })
  
  
  return(corescore_vec)
}

dagwishart_score_node <- function(j,parentnodes,N,param){
  UN<-param$UN
  U0 <- param$U0
  
  awpN_new <- param$alpha_post[j]
  aw_new <- awpN_new-N
  
  scoreconstvec<-param$scoreconstlist
  
  lp<-length(parentnodes) #number of parents
  
  logedgepvec <- param$logedgepvec
  
  awpd2_new <- (aw_new - lp)/2 - 1
  awpNd2_new <- (awpN_new - lp)/2 - 1
  
  A<-UN[j,j]
  A0 <- U0[j,j]
  switch(as.character(lp),
         "0"={# just a single term if no parents
           corescore <- scoreconstvec[[lp+1]][j] -awpNd2_new*log(A)+awpd2_new*log(A0)
         },
         
         "1"={# no need for matrices
           D<-UN[parentnodes,parentnodes]
           D0 <- U0[parentnodes, parentnodes]
           logdetD<-log(D)
           logdetD0 <- log(D0)
           B<-UN[j,parentnodes]
           B0 <- U0[j, parentnodes]
           logdetpart2<-log(A-B^2/D)   
           logdetpart2_0 <- log(A0-B0^2/D0)
           corescore <- scoreconstvec[[lp+1]][j]-awpNd2_new*logdetpart2 - logdetD/2 + 
             awpd2_new*logdetpart2_0 + logdetD0/2
           if (!is.null(param$logedgepvec)) { # if there is an additional edge penalisation
             corescore <- corescore + param$logedgepvec[lp+1]
           }
         },
         
         "2"={# can do matrix determinant and inverse explicitly
           D<-UN[parentnodes,parentnodes]
           D0 <- U0[parentnodes, parentnodes]
           detD<-dettwobytwo(D)
           detD0 <- dettwobytwo(D0)
           logdetD<-log(detD)
           logdetD0 <- log(detD0)
           B<-UN[j,parentnodes]
           B0 <- U0[j, parentnodes]
           logdetpart2<-log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
           logdetpart2_0<-log(dettwobytwo(D0-(B0)%*%t(B0)/A0))+log(A0)-logdetD0
           corescore <- scoreconstvec[[lp+1]][j]-awpNd2_new*logdetpart2 - logdetD/2 +
             awpd2_new*logdetpart2_0 + logdetD0/2
           if (!is.null(param$logedgepvec)) { # if there is an additional edge penalisation
             corescore <- corescore + logedgepvec[lp+1]
           }
         },
         
         {# otherwise we use cholesky decomposition to perform both
           D<-as.matrix(UN[parentnodes,parentnodes])
           D0 <- as.matrix(U0[parentnodes,parentnodes])
           choltemp<-chol(D)
           choltemp0<-chol(D0)
           logdetD<-2*sum(log(choltemp[(lp+1)*c(0:(lp-1))+1]))
           logdetD0<-2*sum(log(choltemp0[(lp+1)*c(0:(lp-1))+1]))
           B<-UN[j,parentnodes]
           B0<-U0[j,parentnodes]
           logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
           logdetpart2_0<-log(A0-sum(backsolve(choltemp0,B0,transpose=TRUE)^2))
           corescore <- scoreconstvec[[lp+1]][j]-awpNd2_new*logdetpart2 - logdetD/2 +
             awpd2_new*logdetpart2_0 + logdetD0/2
           if (!is.null(param$logedgepvec)) { # if there is an additional edge penalisation
             corescore <- corescore + logedgepvec[lp+1]
           }
         })
  return(corescore)
}

dagwishart_score_plus_parent <- function(j, parentnodes, plus_parentnodes, N, param){
  UN<-param$UN
  U0 <- param$U0
  
  awpN_new <- param$alpha_post[j]
  aw_new <- awpN_new-N
  
  scoreconstvec<-param$scoreconstlist
  
  logedgepvec <- param$logedgepvec
  
  lp<-length(parentnodes) #number of parents
  lpp <- length(plus_parentnodes)
  corescore_vec <- numeric(lpp+1)
  
  awpd2_new <- (awpN_new - N - lp)/2 - 1
  awpNd2_new <- (awpN_new - lp)/2 - 1
  A<-UN[j,j]
  A0 <- U0[j,j]
  
  UN_noplus <- UN[parentnodes, parentnodes]
  UN_j_noplus <- UN[j, parentnodes]
  UN_j_plus <- UN[j, plus_parentnodes]
  UN_plus_off <- UN[parentnodes, plus_parentnodes]
  UN_plus_on <- UN[plus_parentnodes, plus_parentnodes]
  UN_plus_diag <- UN_plus_on[(lpp+1)*c(0:(lpp-1))+1]
  
  U0_noplus <- U0[parentnodes, parentnodes]
  U0_j_noplus <- U0[j, parentnodes]
  U0_j_plus <- U0[j, plus_parentnodes]
  U0_plus_off <- U0[parentnodes, plus_parentnodes]
  U0_plus_on <- U0[plus_parentnodes, plus_parentnodes]
  U0_plus_diag <- U0_plus_on[(lpp+1)*c(0:(lpp-1))+1]
  
  switch(as.character(lp),
         "0"={# just a single term if no parents
           corescore_vec[1] <- scoreconstvec[[lp+1]][j] -awpNd2_new*log(A)+awpd2_new*log(A0)
           D<-UN_plus_diag
           logdetD<-log(D)
           B<-UN_j_plus
           logdetpart2<-log(A-B^2/D)
           
           D0<-U0_plus_diag
           logdetD0<-log(D0)
           B0<-U0_j_plus
           logdetpart2_0<-log(A0-B0^2/D0)
           corescores <- scoreconstvec[[lp+2]][j]-(awpNd2_new+1/2)*logdetpart2 - logdetD/2 +
             awpd2_new*logdetpart2_0 + logdetD0/2
           if (!is.null(param$logedgepvec)) { # if there is an additional edge penalisation
             corescores <- corescores + as.numeric(param$logedgepvec[lp+1])
           }
           corescore_vec[2:(lpp+1)] <- corescores
         },
         
         "1"={# no need for matrices
           D<-UN_noplus
           logdetD<-log(D)
           B<-UN_j_noplus
           logdetpart2<-log(A-B^2/D)
           
           D0<-U0_noplus
           logdetD0<-log(D0)
           B0<-U0_j_noplus
           logdetpart2_0<-log(A0-B0^2/D0)
           
           corescore <- scoreconstvec[[lp+1]][j]-awpNd2_new*logdetpart2 - logdetD/2 +
             awpd2_new*logdetpart2_0 + logdetD0/2
           
           logdetD_plus <- log(D*UN_plus_diag - as.numeric(UN_plus_off^2))
           logdetpart2_plus <- log((D-B^2/A)*(UN_plus_diag-UN_j_plus^2/A)-
                                     (UN_plus_off-UN_j_plus*B/A)^2)+log(A)-logdetD_plus
           
           logdetD0_plus <- log(D0*U0_plus_diag - as.numeric(U0_plus_off^2))
           logdetpart2_0_plus <- log((D0-B0^2/A0)*(U0_plus_diag-U0_j_plus^2/A0)-
                                       (U0_plus_off-U0_j_plus*B0/A0)^2)+log(A0)-logdetD0_plus
           
           corescores <- scoreconstvec[[lp+2]][j]-(awpNd2_new+1/2)*logdetpart2_plus - logdetD_plus/2 +
             awpd2_new*logdetpart2_0_plus + logdetD0_plus/2
           if (!is.null(param$logedgepvec)) { # if there is an additional edge penalization
             penalty <- param$logedgepvec[lp+1]
             plus_penalties <- as.numeric(param$logedgepvec[lp+2])
             corescore <- corescore + penalty
             corescores <- corescores + plus_penalties
           }
           corescore_vec[1] <- corescore
           corescore_vec[2:(lpp+1)] <- corescores
         },
         
         {# otherwise we use cholesky decomposition to perform both
           D <- UN_noplus
           choltemp<-chol(D)
           logdetD<-2*sum(log(choltemp[(lp+1)*c(0:(lp-1))+1]))
           B<-UN_j_noplus
           c_noplus <- backsolve(choltemp,B,transpose=TRUE)
           logdetpart2<-log(A-sum(c_noplus^2))
           
           D0 <- U0_noplus
           choltemp0<-chol(D0)
           logdetD0<-2*sum(log(choltemp0[(lp+1)*c(0:(lp-1))+1]))
           B0<-U0_j_noplus
           c_noplus0 <- backsolve(choltemp0,B0,transpose=TRUE)
           logdetpart2_0<-log(A0-sum(c_noplus0^2))
           
           corescore <- scoreconstvec[[lp+1]][j]-awpNd2_new*logdetpart2 - logdetD/2 +
             awpd2_new*logdetpart2_0 + logdetD0/2
           
           choltemp_new_12 <- backsolve(choltemp, UN_plus_off, transpose=TRUE)
           choltemp_new_22 <- sqrt(UN_plus_diag - colsums(choltemp_new_12^2))
           c_plusses <- as.numeric((UN_j_plus-t(choltemp_new_12) %*% c_noplus) / choltemp_new_22)
           logdetD_plus <- logdetD + 2*log(choltemp_new_22)
           logdetpart2_plus <- log(A-sum(c_noplus^2)-c_plusses^2)
           
           choltemp0_new_12 <- backsolve(choltemp0, U0_plus_off, transpose=TRUE)
           choltemp0_new_22 <- sqrt(U0_plus_diag - colsums(choltemp0_new_12^2))
           c_plusses0 <- as.numeric((U0_j_plus-t(choltemp0_new_12) %*% c_noplus0) / choltemp0_new_22)
           logdetD0_plus <- logdetD0 + 2*log(choltemp0_new_22)
           logdetpart2_0_plus <- log(A0-sum(c_noplus0^2)-c_plusses0^2)
           
           corescores <- scoreconstvec[[lp+2]][j]-(awpNd2_new+1/2)*logdetpart2_plus - logdetD_plus/2 +
             awpd2_new*logdetpart2_0_plus + logdetD0_plus/2
           if (!is.null(param$logedgepvec)) { # if there is an additional edge penalisation
             penalty <- sum(param$logedgepvec[lp+1])
             plus_penalties <- as.numeric(param$logedgepvec[lp+2])
             corescore <- corescore + penalty
             corescores <- corescores + plus_penalties
           }
           corescore_vec[1] <- corescore
           corescore_vec[2:(lpp+1)] <- corescores
           
         })
  
  return(corescore_vec)
}

#description: Updates the order from a move type (mimics BiDAG implementation)
#' @param prec_t                  order at current step
#' @param move_type               move type for constructing proposals. Types of moves:
#                                   1. "local transposition" - swaps adjacent nodes
#                                   2. "global swap" - swaps 2 nodes
#                                   3. "node relocation" - moves a single node
#' @param space_banned_score_list banned score table for scoring orders
#' @param map_pars                hash table for quick scoring
#' @param banned_pars             hash table for quick banned scoring
#' @param order_score             nodewise scores for current order
#' @param plus1                   indicator for whether to consider (+1) scores for order moves
#' @param space_banned_plus_list  banned plus score table for (+1) scoring orders
#' @param H                       search space
implement_order <- function(prec_t, move_type,
                            space_banned_score_list, map_pars,
                            banned_pars, order_score,
                            plus_1=FALSE, space_banned_plus_list,
                            H=NULL){
  N <- length(prec_t)
  
  if(move_type %in% c("local transposition", "global swap")){
    if(move_type=="local transposition"){
      node1 <- sample(1:(N-1), size=1)
      changed_nodes <- c(node1, node1+1)
    } else {
      changed_nodes <- sample(1:N, 2, replace=FALSE)
    }
    moved_ids <- prec_t[changed_nodes]           # actual node identities that moved
    prec_tplus1 <- prec_t
    prec_tplus1[changed_nodes] <- prec_t[rev(changed_nodes)]
    
    update_nodes <- nodes_affected_by_order_move(moved_ids, H)
    
    banned_pars_new <- banned_parents_mapping(map_pars, prec_tplus1,
                                              create_plus_sets = plus_1,
                                              create_minus_sets = FALSE,
                                              updatenodes = update_nodes,
                                              old_map = banned_pars)
    
    curr_score <- calculate_order_score(banned_pars_new, space_banned_score_list, plus_1,
                                        space_banned_plus_list)
    
    return(list(order=prec_tplus1, score=curr_score, banned_pars=banned_pars_new))
  }
  
  if(move_type=="node relocation"){
    #sampling a node which will be moved to create the proposal
    node1 <- sample(1:N, size=1)
    loc_n1 <- which(prec_t==node1)
    if(!plus_1){
      return(sample_from_node_relocation(prec_t, loc_n1, banned_pars, order_score,
                                         space_banned_score_list, map_pars))
    }
    return(sample_from_node_relocation(prec_t, loc_n1, banned_pars, order_score,
                                       space_banned_score_list, map_pars, TRUE,
                                       space_banned_plus_list))
  }
}
# Helper function for expanding moved node(s) into the full set of nodes needing an update
#' @param moved_nodes node ids whose position in the order changed
#' @param H           current search space adjacency matrix
#'                     (H[i, j] == 1  <=>  j is currently a candidate parent of i)
nodes_affected_by_order_move <- function(moved_nodes, H){
  affected <- moved_nodes
  for(m in moved_nodes){
    # anyone who currently lists m as a candidate parent needs a refresh,
    # since whether m precedes or follows them in the order may have flipped
    affected <- c(affected, which(H[, m] == 1))
  }
  return(unique(affected))
}

#description: sampling step from M-H to choose between current order and the 
#             proposal, via the banned parent table
#' @param prec_t                  order at current step
#' @param prec_prime              proposed order
#' @param space_banned_score_list banned score table for scoring orders
#' @param map_pars                hash table for quick scoring
sample_from_2_orders <- function(prec_t, prec_prime, 
                                 space_banned_score_list, map_pars,
                                 plus_1 = FALSE, space_banned_plus_list=NULL){
  N <- length(space_banned_score_list)
  score_t <- 0
  score_prime <- 0
  banned_pars_t <- banned_parents_mapping(map_pars, prec_t)
  index_banned_t <- banned_pars_t$banned_row
  banned_pars_prime <- banned_parents_mapping(map_pars, prec_prime)
  index_banned_prime <- banned_pars_prime$banned_row
  for(i in 1:N){
    if(!plus_1){
      # access lookup table for current node
      lookup_table <- space_banned_score_list[[i]]
      score_t <- score_t + lookup_table[index_banned_t[i]]
      score_prime <- score_prime + lookup_table[index_banned_prime[i]]
    }
    else{
      # access lookup table for current node
      lookup_table <- space_banned_plus_list[[i]]
      index_allowed_nonpar_t <- banned_pars_t$allowed_nonpar_idx[[i]]
      index_allowed_nonpar_prime <- banned_pars_prime$allowed_nonpar_idx[[i]]
      if(length(index_allowed_nonpar_t)>0){
        score_t <- score_t + logSumExp(lookup_table[index_banned_t[i], 
                                                    c(1, index_allowed_nonpar_t+1)])
      }
      else{
        score_t <- score_t + lookup_table[index_banned_t[i], 1]
      }
      if(length(index_allowed_nonpar_prime)>1){
        score_prime <- score_prime + 
          logSumExp(lookup_table[index_banned_prime[i], c(1, index_allowed_nonpar_prime+1)])
      }
      else{
        score_prime <- score_prime + lookup_table[index_banned_prime[i], 1]
      }
    }
  }
  r <- min(exp(score_prime-score_t), 1)
  prec_t1 <- sample(c("prime", "t"), size=1, prob=c(r, 1-r))
  if(prec_t1=="t"){return(prec_t)}
  return(prec_prime)
}


#description: extension of the above function for more than 2 orders (in node relocation). 
#             This is useful for implementing proposals that are proportional 
#             to the posterior
#' @param prec_orig               original order
#' @param location                location of the node to be relocated
#' @param order_score             score values for the current order
#' @param space_banned_score_list banned score table for scoring orders
#' @param map_pars                hash table for quick scoring
#' @param plus_1                  whether to score based on the +1 list or not
#' @param space_banned_plus_list  (+1) banned score table for scoring (+1) orders
sample_from_node_relocation <- function(prec_orig, location,
                                        banned_pars, order_score,
                                        space_banned_score_list,
                                        map_pars,plus_1 = FALSE, 
                                        space_banned_plus_list=NULL){
  relocated_node <- prec_orig[location]
  N <- length(prec_orig)
  score_mat <- matrix(0, nrow=N, ncol=N)
  score_mat[location,] <- order_score
  
  if(location > 1){
    for(i in (location-1):1){
      new_node <- prec_orig[i]
      prec_new <- prec_orig
      prec_new[i] <- relocated_node
      prec_new[(i+1):location] <- prec_orig[i:(location-1)]
      if(location < N){
        score_mat[i, (location+1):N] <- score_mat[location, (location+1):N]
      }
      if(i > 1){
        score_mat[i, 1:(i-1)] <- score_mat[location, 1:(i-1)]
      }
      if(i < location-1){
        score_mat[i, (i+2):(location)] <- score_mat[i+1, (i+2):(location)]
      }
      allparents_i <- map_pars$par_names[[relocated_node]]
      allparents_iplus1 <- map_pars$par_names[[new_node]]
      nonparents_i <- which(!(1:N %in% c(allparents_i, relocated_node)))
      nonparents_iplus1 <- which(!(1:N %in% c(allparents_iplus1, new_node)))
      bannedrow_i <- 1
      bannedrow_iplus1 <- 1
      
      if(i > 1){
        bannedvalues_i <- which(allparents_i %in% prec_new[1:(i-1)])
        if(length(bannedvalues_i)>0 & !is.na(bannedvalues_i[1])){
          bannedrow_i <- map_pars$maps[[relocated_node]]$backwards[sum(2^bannedvalues_i)/2+1]
        }
      }
      bannedvalues_iplus1 <- which(allparents_iplus1 %in% prec_new[1:i])
      if(length(bannedvalues_iplus1)>0 & !is.na(bannedvalues_iplus1[1])){
        bannedrow_iplus1 <- map_pars$maps[[new_node]]$backwards[sum(2^bannedvalues_iplus1)/2+1]
      }
      if(!plus_1){
        score_mat[i, i] <- space_banned_score_list[[relocated_node]][bannedrow_i,1]
        score_mat[i, i+1] <- space_banned_score_list[[new_node]][bannedrow_iplus1,1]
      }
      else{
        allowedcol_i <- integer(0)
        allowedcol_iplus1 <- integer(0)
        if(i < N){
          allowedcol_i <- which(nonparents_i %in% prec_new[(i+1):N])
          if(i < N-1){
            allowedcol_iplus1 <- which(nonparents_iplus1 %in% prec_new[(i+2):N])
          }
        }
        score_mat[i, i] <- logSumExp(
          space_banned_plus_list[[relocated_node]][bannedrow_i, c(1, allowedcol_i+1)])
        score_mat[i, i+1] <- logSumExp(
          space_banned_plus_list[[new_node]][bannedrow_iplus1, c(1, allowedcol_iplus1+1)])
      }
    }
    
  }
  if(location < N){
    for(i in (location+1):N){
      new_node <- prec_orig[i]
      prec_new <- prec_orig
      prec_new[i] <- relocated_node
      prec_new[location:(i-1)] <- prec_orig[(location+1):i]
      
      score_mat[i, 1:(location-1)] <- score_mat[location, 1:(location-1)]
      if(i < N){
        score_mat[i, (i+1):N] <- score_mat[location,(i+1):N]
      }
      if(i > location+1){
        score_mat[i, location:(i-2)] <- score_mat[i-1, location:(i-2)]
      }
      allparents_i <- map_pars$par_names[[relocated_node]]
      allparents_iminus1 <- map_pars$par_names[[new_node]]
      nonparents_i <- which(!(1:N %in% c(allparents_i, relocated_node)))
      nonparents_iminus1 <- which(!(1:N %in% c(allparents_iminus1, new_node)))
      bannedrow_i <- 1
      bannedrow_iminus1 <- 1
      bannedvalues_i <- which(allparents_i %in% prec_new[1:(i-1)])
      bannedvalues_iminus1 <- which(allparents_iminus1 %in% prec_new[1:(i-2)])
      if(length(bannedvalues_i)>0 & !is.na(bannedvalues_i[1])){
        bannedrow_i <- map_pars$maps[[relocated_node]]$backwards[sum(2^bannedvalues_i)/2+1]
      }
      if((i-1)> 1 & length(bannedvalues_iminus1)>0 & !is.na(bannedvalues_iminus1[1])){
        bannedrow_iminus1 <- map_pars$maps[[new_node]]$backwards[sum(2^bannedvalues_iminus1)/2+1]
      }
      if(!plus_1){
        score_mat[i, i] <- space_banned_score_list[[relocated_node]][bannedrow_i,1]
        score_mat[i, i-1] <- space_banned_score_list[[new_node]][bannedrow_iminus1,1]
      }
      else{
        allowedcol_i <- integer(0)
        allowedcol_iminus1 <- integer(0)
        if(i < N){
          allowedcol_i <- which(nonparents_i %in% prec_new[(i+1):N])
          if(i < N-1){
            allowedcol_iminus1 <- which(nonparents_iminus1 %in% prec_new[i:N])
          }
        }
        score_mat[i, i] <- logSumExp(
          space_banned_plus_list[[relocated_node]][bannedrow_i, c(1, allowedcol_i+1)])
        score_mat[i, i-1] <- logSumExp(
          space_banned_plus_list[[new_node]][bannedrow_iminus1, c(1, allowedcol_iminus1+1)])
      }
    }
  }
  score_vec <- rowsums(score_mat)
  score_all <- logSumExp(score_vec)
  score_vec <- exp(score_vec-score_all)
  new_spot <- sample(1:N, size=1, prob=score_vec)
  prec_new <- prec_orig
  if(new_spot==location){
    prec_new <- prec_orig
  }
  else if(new_spot>location){
    prec_new[new_spot] <- relocated_node
    prec_new[location:(new_spot-1)] <- prec_orig[(location+1):new_spot]
    
  }
  else{
    prec_new[new_spot] <- relocated_node
    prec_new[(new_spot+1):location] <- prec_orig[new_spot:(location-1)]
  }
  return(list(order=prec_new, score=score_mat[new_spot,]))
}


calculate_order_score <- function(banned_pars,space_banned_score_list,
                                  plus_1 = FALSE, space_banned_plus_list=NULL){
  N <- length(banned_pars$banned_row)
  score_vec <- vector(mode="numeric", length=N)
  for(i in 1:N){
    if(!plus_1){
      score_vec[i] <- space_banned_score_list[[i]][banned_pars$banned_row[i],1]
    }
    else{
      score_vec[i] <- logSumExp(
        space_banned_plus_list[[i]][banned_pars$banned_row[i], 
                                    c(1, banned_pars$allowed_nonpar_idx[[i]]+1)])
    }
  }
  return(score_vec)
}

#description: keeps track of the row compatible with the current order 
#             for each node's banned parent table (and score table)
#' @param map_pars          hash table for quick scoring
#' @param prec              current_order
#' @param create_plus_sets  logical; if TRUE, also builds `valid_pset_rows`: for
#'                          each node, the row indices into that node's "+1"
#'                          parent-subset score table (map_pars$par_pset) that
#'                          remain valid under the current order (excluding any
#'                          parent subset that would require a candidate parent
#'                          to sit after the node in the order). Needed when
#'                          scoring single-edge-addition ("+1"/birth) proposals;
#'                          leave FALSE to skip that work when not needed.
#' @param create_minus_sets logical; if TRUE, also builds `banned_minus_rows`:
#'                          for each node, one banned-row lookup index per
#'                          *currently included* candidate parent, giving the
#'                          banned-row index that would result if that single
#'                          parent were removed. Needed when scoring single-edge
#'                          -deletion ("-1"/death) proposals; leave FALSE to
#'                          skip that work when not needed.
#' @param updatenodes       NULL (default) Otherwise, only these node ids
#'                          are recomputed; everything else is copied from `old_map`.
#' @param old_map           the previously-computed banned_parents_mapping()
#'                          output to update incrementally. Required whenever
#'                          `updatenodes` is supplied.
banned_parents_mapping <- function(map_pars, prec, create_plus_sets=FALSE,
                                   create_minus_sets=FALSE,
                                   updatenodes=NULL, old_map=NULL){
  N <- length(map_pars$par_names)
  
  full_recompute <- is.null(updatenodes) || is.null(old_map)
  if(full_recompute) updatenodes <- 1:N
  
  # seed containers: reuse `old_map` for incremental updates, fresh otherwise
  if(!full_recompute){
    banned_lookup_row <- old_map$banned_row
    allowed_nonpar    <- old_map$allowed_nonpar_idx
    if(create_minus_sets) minus_lookup_rows <- old_map$banned_minus_rows
    if(create_plus_sets)  valid_parset_maps  <- old_map$valid_pset_rows
  } else {
    banned_lookup_row <- vector(length=N)
    allowed_nonpar    <- vector(mode="list", length=N)
    if(create_minus_sets) minus_lookup_rows <- vector(mode="list", length=N)
    if(create_plus_sets)  valid_parset_maps  <- vector(mode="list", length=N)
  }
  
  for(i in updatenodes){
    # find where the current node is in the order
    index_i <- which(prec==i)
    # find the banned parent sets 
    # (a.k.a, nodes listed after the current node for each order)
    all_parents <- map_pars$par_names[[i]]
    non_parents <- which(!(1:N %in% c(all_parents, i)))
    banned_par_idx <- integer(0)
    allowed_nonpar_idx <- integer(0)
    if(index_i != 1){
      banned_par_idx <- which(all_parents %in% prec[1:(index_i-1)])
    }
    if(index_i != N){
      allowed_nonpar_idx <- which(non_parents %in% prec[(index_i+1):N])
    }
    allowed_nonpar[[i]] <- allowed_nonpar_idx
    
    # access the relevant row from the banned parent table
    index_banned <- ifelse(length(banned_par_idx)==0 | is.na(banned_par_idx[1]),1,
                           map_pars$maps[[i]]$backwards[sum(2^banned_par_idx)/2+1])
    banned_lookup_row[i] <- index_banned
    
    if(create_minus_sets){
      N_parents <- length(all_parents)
      banned_row_minus <- integer(N_parents)
      if(N_parents > 0){
        powers_of_2 <- 2^(0:(N_parents - 1))
        banned_mask <- logical(N_parents)
        banned_mask[banned_par_idx] <- TRUE
        banned_sum <- sum(powers_of_2[banned_mask])
        for(j in 1:N_parents){
          if(banned_mask[j]){
            new_sum <- banned_sum
          } else {
            new_mask <- banned_mask
            new_mask[j] <- TRUE
            new_sum <- sum(powers_of_2[new_mask])
          }
          index_banned2 <- if(new_sum==0) 1 else map_pars$maps[[i]]$backwards[new_sum+1]
          banned_row_minus[j] <- index_banned2
        }
      } else {
        banned_row_minus <- 1
      }
      minus_lookup_rows[[i]] <- banned_row_minus
    }
    
    if(create_plus_sets){
      #finding which indices in the score table are valid per the order
      if(index_banned==1){
        allowed_rows <- c(1:nrow(map_pars$par_pset[[i]]))
      } else {
        tablesize <- dim(map_pars$par_pset[[i]])
        if(tablesize[1]==1 || length(banned_par_idx)==tablesize[2]){
          allowed_rows <- c(1)
        } else {
          allowed_rows <- c(2:tablesize[1])
          banned_pars_i <- map_pars$par_names[[i]][banned_par_idx]
          for(j in 1:tablesize[2]){
            banned_rows <- which(map_pars$par_pset[[i]][allowed_rows, j] %in% banned_pars_i)
            if(length(banned_rows)>0){allowed_rows <- allowed_rows[-banned_rows]}
          }
          allowed_rows <- c(1, allowed_rows)
        }
      }
      valid_parset_maps[[i]] <- allowed_rows
    }
  }
  
  if(create_minus_sets){
    if(create_plus_sets){
      return(list("banned_row"=banned_lookup_row, "banned_minus_rows"=minus_lookup_rows,
                  "valid_pset_rows"=valid_parset_maps,"allowed_nonpar_idx"=allowed_nonpar))
    }
    return(list("banned_row"=banned_lookup_row, "banned_minus_rows"=minus_lookup_rows,
                "allowed_nonpar_idx"=allowed_nonpar))
  }
  if(create_plus_sets){
    return(list("banned_row"=banned_lookup_row, "valid_pset_rows"=valid_parset_maps,
                "allowed_nonpar_idx"=allowed_nonpar))
  }
  return(list("banned_row"=banned_lookup_row, "allowed_nonpar_idx"=allowed_nonpar))
}


#description: samples a graph from an order score list
#' @param score_list list of parent scores compatible with order
sample_graph <- function(score_list, prec, map_pars, banned_map_pars,
                         plus_1 = FALSE, score_plus_list = NULL, banned_plus_list = NULL,
                         plus_pars = NULL){
  if(!plus_1){
    return(sample_graph_noplus(score_list, prec, map_pars, banned_map_pars))
  }
  else{
    return(sample_plus_graph(score_plus_list, prec, map_pars,
                             banned_plus_list, banned_map_pars, plus_pars))
  }
}

#description: samples a graph from an order score list
#' @param score_list list of parent scores compatible with order
sample_graph_noplus <- function(score_list, prec, map_pars, banned_map_pars){
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
    
    prob_vec <- exp(lookup_table-logSumExp(lookup_table))
    skel_idx <- sample.int(length(prob_vec), size=1, prob=prob_vec)
    if(skel_idx != 1){
      parents <- par_sets[skel_idx,1:numpars_vec[skel_idx]]
      G[i, parents] <- 1
    }
  }
  return(G)
}

sample_plus_graph <- function(score_plus_list, prec, map_pars,
                              banned_plus_list, banned_map_pars, plus_pars){
  N <- length(score_plus_list)
  G <- matrix(0, nrow=N, ncol=N)
  for(i in 1:N){
    N_plus_sets <- length(banned_map_pars$allowed_nonpar_idx[[i]])+1
    if(N_plus_sets==1){
      out_idx <- 1
      out_parents <- integer(0)
      out_col <- 1
    }
    else{
      allowed_columns <- c(1, banned_map_pars$allowed_nonpar_idx[[i]]+1)
      tot_scores <- banned_plus_list[[i]][banned_map_pars$banned_row[i], allowed_columns]
      prob_vec <- exp(tot_scores - logSumExp(tot_scores))
      out_idx <- sample.int(N_plus_sets, size=1, prob=prob_vec)
      out_parents <- plus_pars$par_pset[[i]][allowed_columns[out_idx], 1:plus_pars$numpars_vec[[i]][allowed_columns[out_idx]]]
      out_col <- allowed_columns[out_idx]
    }
    valid_rows <- banned_map_pars$valid_pset_rows[[i]]
    
    lookup_table <- score_plus_list[[i]][valid_rows, out_col]
    numpars_vec <- map_pars$numpars_vec[[i]][valid_rows]
    if(numpars_vec[length(numpars_vec)]==1){
      par_sets <- as.matrix(map_pars$par_pset[[i]][valid_rows,], ncol=1)
    }
    else{par_sets <- map_pars$par_pset[[i]][valid_rows,]}
    prob_vec2 <- as.numeric(exp(lookup_table-logSumExp(lookup_table)))
    skel_idx <- sample.int(length(prob_vec2), size=1, prob=prob_vec2)
    if(skel_idx != 1){
      skel_parents <- par_sets[skel_idx,1:numpars_vec[skel_idx]]
      parents <- c(out_parents, skel_parents)
      G[i, parents] <- 1
    }
    else if(length(out_parents)> 0){
      parents <- out_parents
      G[i, parents] <- 1
    }
    
  }
  return(G)
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


#description: expands the search space by adding any edges from a drawn set 
#             of graphs that are not currently in the search space
#' @param H               current search space
#' @param edges           set of edges to add to search space
#' @param sparsity_limit  cap on the number of edges allowed
expand_search_space <- function(H, edges, sparsity_limit=18){
  
  H_new <- H
  H_new[edges] <- 1
  
  sparsity_count <- rowSums(H_new)
  sparsity_reached <- which(sparsity_count > sparsity_limit)
  if(length(sparsity_reached)>0){
    H_new[sparsity_reached,] <- H[sparsity_reached,]
  }
  
  added_nodes <- H_new - H
  update_nodes <- c(1:nrow(H))[rowSums(added_nodes)>0]
  
  return(list(H_new=H_new, updatenodes = update_nodes))
}

#description: contracts the search space by adding any edges from a drawn set 
#             of graphs that are not currently in the search space
#' @param H     current search space
#' @param edges set of edges to add to search space
shrink_search_space <- function(H, edges){
  
  H_new <- H
  H_new[edges] <- 0
  
  removed_nodes <- H - H_new
  update_nodes <- c(1:nrow(H))[rowSums(removed_nodes)>0]
  
  return(list(H_new=H_new, updatenodes = update_nodes))
}

calculate_log_birth_rate <- function(H, banned_plus_list, sparsity, has_order=FALSE, 
                                     prec=NULL, plus_map_pars=NULL, banned_map_pars=NULL,
                                     rounded=FALSE){
  N <- length(banned_plus_list)
  outside_matrix <- 1-H-diag(N)
  update_matrix <- outside_matrix
  update_matrix[outside_matrix==0] <- -1e12
  curr_setsize <- rowsums(H)
  birth_allowed <- curr_setsize < sparsity
  for(i in 1:N){
    if(birth_allowed[i]){
      N_plus_sets <- ncol(banned_plus_list[[i]])
      if(has_order){
        tot_scores <- banned_plus_list[[i]][banned_map_pars$banned_row[i],]
        index_i <- which(prec==i)
        allowed_indices <- banned_map_pars$allowed_nonpar_idx[[i]]
        tot_scores_new <- tot_scores[-1]
        tot_scores_orig <- tot_scores[1]
        if(length(allowed_indices)>0){
          if(length(allowed_indices)<length(tot_scores_new)){
            tot_scores_new[-allowed_indices] <- -1e12
          }
        }
        else{
          tot_scores_new <- rep(-1e12, length(tot_scores_new))
        }
      }
      else{
        tot_scores <- banned_plus_list[[i]][1,]
        tot_scores_new <- tot_scores[-1]
        tot_scores_orig <- tot_scores[1]
      }
      tot_scores_matr <- matrix(tot_scores_orig, nrow=length(tot_scores_new), ncol=2)
      tot_scores_matr[,1] <- tot_scores_new
      tot_scores_all <- rowLogSumExps(tot_scores_matr)
      
      B_e <- as.numeric(-log(2) + tot_scores_all - tot_scores_orig)
      if(rounded){
        B_e <- ifelse(B_e>0, 0, B_e)
      }
      update_matrix[i,which(outside_matrix[i,]==1)] <- B_e
    }
    else{
      update_matrix[i, which(outside_matrix[i,]==1)] <- -1e12
    }
  }
  return(update_matrix)
}

calculate_log_death_rate <- function(H, score_list, map_pars, banned_score_list,
                                     has_order=FALSE, prec=NULL, banned_map_pars=NULL,
                                     c_star = 1, rounded=FALSE){
  if(has_order & length(banned_map_pars$banned_minus_rows)==0){
    banned_map_pars <- banned_parents_mapping(map_pars, prec, create_minus_sets=TRUE)
  }
  
  N <- nrow(H)
  update_matrix <- H
  update_matrix[H==0] <- -1e12
  for(i in 1:N){
    if(has_order){
      n_parents <- length(map_pars$par_names[[i]])
      tot_scores <- banned_score_list[[i]][banned_map_pars$banned_minus_rows[[i]],1]
      orig_score <- banned_score_list[[i]][banned_map_pars$banned_row[i],1]
    }
    
    else{
      n_parents <- length(map_pars$par_names[[i]])
      orig_score <- banned_score_list[[i]][1,1]
      if(n_parents > 0){
        tot_scores <- banned_score_list[[i]][2:(n_parents+1),1]
      }
      else{
        tot_scores <- NA
      }
    }
    
    # for numerical stability, we divide by the original score, or equivalent, 
    # subtract log original. This allows us to cancel out all terms for other nodes
    c_adjust <- log(2*c_star)
    D_e <- as.numeric(c_adjust+tot_scores - orig_score)
    if(rounded){
      D_e <- ifelse(D_e>0, 0, D_e) 
    }
    update_matrix[i, as.numeric(which(H[i,]==1))] <- D_e
  }
  return(update_matrix)
}

sample_DL_parameters <- function(G, prec, param){
  N <- nrow(G)
  n_pars <- rowsums(G)
  if(!is.null(param$UN)){
    U_N <- param$UN
  }
  else if(!is.null(param$TN)){
    U_N <- param$TN
  }
  else{
    stop("Please use a score with a valid precision matrix")
  }
  if(!is.null(param$alpha_vec)){
    alpha_vec <- param$alpha_vec
  }
  else if(!is.null(param$aw)){
    alpha_vec <- param$aw - N + 2*n_pars + 3
  }
  else{
    stop("Please use a score with a valid precision matrix")
  }
  
  D_vec <- numeric(N)
  L_matr <- matrix(0, nrow=N, ncol=N)
  
  for(i in 1:N){
    pars_i <- which(G[i,]==1)
    shape_i <- (alpha_vec[i] - n_pars[i])/2-1
    if(n_pars[i] > 0){
      U_N_pars_i_inv <- solve(U_N[pars_i,pars_i])
      scale_i <- 0.5 * (U_N[i,i] - U_N[i, pars_i] %*% U_N_pars_i_inv %*% U_N[pars_i, i])
    }
    else{
      scale_i <- 0.5 * U_N[i,i]
    }
    
    v_i <- invgamma::rinvgamma(1, shape=shape_i, scale=scale_i)
    
    if(n_pars[i] > 0){
      mu_i <- - U_N_pars_i_inv %*% U_N[pars_i, i]
      Sigma_i <- v_i * U_N_pars_i_inv
      
      b_vec_i <- Rfast::rmvnorm(1, mu_i, Sigma_i)
      L_matr[i, pars_i] <- b_vec_i
    }
    
    D_vec[i] <- v_i
    
  }
  return(list("D" = D_vec, "L" = L_matr))
}