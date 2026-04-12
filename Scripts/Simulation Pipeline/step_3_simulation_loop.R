parallelized_job <- function(seed, method, n, sim_root, pc_alg_thresh,
                             dataset_multiplier, data_model,
                             replicates, thinning_rate_full, burn_in_rate, 
                             max_spars, max_ch, space_change_prob, save_waiting_times,
                             Bs) {
  
  set.seed(seed)
  
  # True Graph Generation
  if (method == "ERs") {
    trueDAGedges <- as(pcalg::randDAG(n = n, d = 4, 
                                      wFUN = list(runif, min = 0.4, max = 2)), "matrix")
  } else if (method == "SBM_2") {
    v1 <- min(0.15, 6/n); v2 <- min(0.08, 4/n); cv <- min(0.01, 2/n)
    dagedges <- sample_sbm(n = n, pref.matrix = matrix(c(v1, cv, cv, v2), nrow = 2, ncol = 2),
                           block.sizes = c(floor(0.8 * n), ceiling(0.2 * n)), directed = TRUE)
    trueDAGedges <- g2Q(dagedges) * runif(n * n, min = 0.4, max = 2)
  } else if (method == "hSBM_3") {
    b1 <- min(4/n, 0.1); b2_v1 <- min(0.15, 6/n); b2_v2 <- min(0.08, 4/n); b2_cv <- min(0.01, 2/n)
    b3_v1 <- min(0.2, 8/n); b3_v2 <- 0; b3_v3 <- min(0.08, 4/n)
    b3_cv1 <- min(0.01, 1/n); b3_cv2 <- min(0.2, 8/n); b3_cv3 <- min(0.08, 4/n)
    dagedges <- sample_hierarchical_sbm(n = n, 
                                        m = c(round(n * 1/10), round(n * 3/10), round(n * 6/10)),
                                        rho = list(1, c(1/3, 2/3), c(1/6, 2/6, 1/2)),
                                        C = list(matrix(b1, nrow = 1),
                                                 matrix(c(b2_v1, b2_cv, b2_cv, b2_v2), nrow = 2, ncol = 2),
                                                 matrix(c(b3_v1, b3_cv1, b3_cv3, b3_cv1, b3_v2, 
                                                          b3_cv2, b3_cv3, b3_cv2, b3_v3), nrow = 3, ncol = 3)),
                                        p = 0.05)
    trueDAGedges <- g2Q(dagedges) * runif(n * n, min = 0.4, max = 2)
  }
  
  # Parallel Loops
  output <- foreach(thresh = pc_alg_thresh, .combine = 'rbind') %:%
    foreach(m = dataset_multiplier, .combine = 'rbind') %:%
    foreach(datatype = data_model, .combine = 'rbind') %:%
    foreach(sparsity_type = c("fixed", "plus1"), .combine = 'rbind') %:%
    foreach(k = replicates, .combine = 'rbind') %dopar% {
      
      trueDAGedges_curr <- spectral_rescale(trueDAGedges)
      bdmcmc_run <- (n <= 100 & k == 1) & sparsity_type=="fixed"
      Bs_2 <- max(Bs, round(n * n * log(n)))
      thinning_constant <- floor(Bs_2 / thinning_rate_full)
      N <- m * n
      
      # Path Management
      folder_name1 <- sprintf("method_%s_n_%d_N_%d", method, n, N)
      folder_name2 <- sprintf("data_%s_thresh_%.3f", datatype, thresh)
      full_path    <- file.path(sim_root, folder_name1, folder_name2)
      if (!dir.exists(full_path)) dir.create(full_path, recursive = TRUE)
      
      file_summary <- file.path(full_path, sprintf("seed_%d_rep_%d_%s_summary.Rdata", seed, k, sparsity_type))
      
      # Skip if already processed
      if (!file.exists(file_summary)) {
        
        curr_files <- dir(full_path)
        
        # Data Generation
        if (datatype == "Gaussian") {
          data <- rmvDAG(trueDAGedges_curr, N, standardise = FALSE)
        } else if (datatype == "FCM") {
          data <- rmvDAG3(trueDAGedges_curr, N, standardise = FALSE)
        } else {
          data <- rmvDAG2(trueDAGedges_curr, N, standardise = FALSE)
        }
        score_par_sim <- scoreparameters("bge", data)
        
        # Sparsity set-up:
        if(sparsity_type=="fixed"){
          max_spars_curr <- max_spars[1]
        }
        else{
          max_spars_curr <- max_spars[2]
        }
        
        # 0. Baseline skeletons (PC-based and GES-based)
        pc_fit <- skeleton(suffStat = list(C = cor(data), n = nrow(data)), indepTest = gaussCItest,
                           alpha = thresh, labels = paste0("V", 1:ncol(data)), method = "stable", m.max = max_spars_curr + 1)
        space_PC <- 1 * as(pc_fit@graph, "matrix")
        save(space_PC, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_PC.Rdata", seed, k, sparsity_type)))
        
        set.seed(seed * 10 + k)
        save(trueDAGedges_curr, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_graph.Rdata", seed, k, sparsity_type)))
        
        ges_score <- new("GaussL0penObsScore", data, lambda = log(nrow(data)))
        ges_fit   <- ges(ges_score, maxDegree = max_spars_curr + 1)
        gesCPDAG  <- 1 * as(ges_fit$essgraph, "matrix")
        
        if(sparsity_type=="fixed"){
          max_spars_pc <- max(max_spars_curr + 1, max(colSums(space_PC)) + 1)
          max_spars_ges <- max(max_spars_curr + 1, max(colSums(gesCPDAG)) + 1)
        }
        else{
          max_spars_pc <- max(colSums(space_PC)) + 1
          max_spars_ges <- max(colSums(gesCPDAG)) + 1
        }
        
        # Algorithm Executions
        
        # 1. BiDAG on PC Space (both plus and nonplus version)
        start_t <- Sys.time()
        results_BIDAG <- orderMCMC(score_par_sim, startspace = space_PC, MAP = FALSE, plus1 = TRUE, 
                                   chainout = TRUE, iterations = Bs_2, stepsave = thinning_constant, hardlimit = max_spars_pc)
        time_BIDAG <- as.numeric(Sys.time() - start_t, units = "secs")
        save(results_BIDAG, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_BiDAG_PC.Rdata", seed, k, sparsity_type)))
        
        start_t <- Sys.time()
        results_BIDAG_noplus <- orderMCMC(score_par_sim, startspace = space_PC, MAP = FALSE, plus1 = FALSE, 
                                          chainout = TRUE, iterations = Bs_2, stepsave = thinning_constant, 
                                          hardlimit = max_spars_pc-1)
        time_BIDAG_noplus <- as.numeric(Sys.time() - start_t, units = "secs")
        save(results_BIDAG_noplus, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_BiDAG_PC_noplus.Rdata", seed, 
                                                                       k, sparsity_type)))
        
        # 2. Iterative Search Space (based on Kuipers et al.)
        start_t <- Sys.time()
        bestDAGs <- iterativeMCMC(score_par_sim, scoreout = TRUE, hardlimit = max_spars_curr, 
                                  softlimit = 10, alpha = (min(0.4, 20/n)/2))
        time_MAP_p1 <- Sys.time() - start_t
        save(bestDAGs, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_MAP.Rdata", seed, k, sparsity_type)))
        
        
        if(sparsity_type=="fixed"){
          max_spars_map <- max(max_spars_curr + 1, max(colSums(bestDAGs$endspace)) + 1)
        }
        else{
          max_spars_map <- max(colSums(bestDAGs$endspace)) + 1
        }
        
        
        # 3. BROOD MCMC Runs
        # BROOD on PC Space
        start_t <- Sys.time()
        results <- graph_mcmc(space_PC, score_par_sim, iter = Bs_2, thinning = thinning_constant,
                              space_move_prob = space_change_prob, max_sparsity = max_spars_pc, 
                              save_all_weights = save_waiting_times, warm_up = floor(burn_in_rate * Bs_2), 
                              max_change = max_ch, sparse = TRUE, verbose = FALSE)
        time_BROOD <- as.numeric(Sys.time() - start_t, units = "secs")
        save(results, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_brood_PC.Rdata", seed, k, sparsity_type)))
        
        # BROOD on Iterative MAP Space
        start_t <- Sys.time()
        results_2 <- graph_mcmc(t(bestDAGs$endspace), score_par_sim, iter = Bs_2, thinning = thinning_constant,
                                space_move_prob = space_change_prob, max_sparsity = max_spars_map, 
                                save_all_weights = save_waiting_times, warm_up = floor(burn_in_rate * Bs_2), 
                                max_change = max_ch, sparse = TRUE, verbose = FALSE)
        time_BROOD2 <- as.numeric(Sys.time() - start_t, units = "secs")
        save(results_2, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_brood_MAP.Rdata", seed, k, sparsity_type)))
        
        # BROOD on GES Space
        start_t <- Sys.time()
        results_3 <- graph_mcmc(t(gesCPDAG), score_par_sim, iter = Bs_2, thinning = thinning_constant,
                                space_move_prob = space_change_prob, max_sparsity = max_spars_ges, 
                                save_all_weights = save_waiting_times, warm_up = floor(burn_in_rate * Bs_2), 
                                max_change = max_ch, sparse = TRUE, verbose = FALSE)
        time_BROOD3 <- as.numeric(Sys.time() - start_t, units = "secs")
        save(results_3, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_brood_GES.Rdata", seed, k, sparsity_type)))
        
        # 4. Further BiDAG Runs
        # BiDAG on Iterative MAP Space (both plus and nonplus version)
        start_t <- Sys.time()
        results_BIDAG_2 <- orderMCMC(score_par_sim, startspace = bestDAGs$endspace, MAP = FALSE, plus1 = TRUE, 
                                     chainout = TRUE, iterations = Bs_2, stepsave = thinning_constant, hardlimit = max_spars_map)
        time_BIDAG2 <- as.numeric((Sys.time() - start_t) + time_MAP_p1, units = "secs")
        save(results_BIDAG_2, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_BiDAG_MAP.Rdata", seed, k, sparsity_type)))
        
        start_t <- Sys.time()
        results_BIDAG_2_noplus <- orderMCMC(score_par_sim, startspace = bestDAGs$endspace, MAP = FALSE, plus1 = FALSE, 
                                            chainout = TRUE, iterations = Bs_2, stepsave = thinning_constant, 
                                            hardlimit = max_spars_map-1)
        time_BIDAG2_noplus <- as.numeric((Sys.time() - start_t) + time_MAP_p1, units = "secs")
        save(results_BIDAG_2_noplus, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_BiDAG_MAP_noplus.Rdata", 
                                                                         seed, k, sparsity_type)))
        
        # BiDAG on GES (both plus and nonplus version)
        start_t <- Sys.time()
        results_BIDAG_3 <- orderMCMC(score_par_sim, startspace = gesCPDAG, MAP = FALSE, plus1 = TRUE, 
                                     chainout = TRUE, iterations = Bs_2, stepsave = thinning_constant, hardlimit = max_spars_ges)
        time_BIDAG3 <- as.numeric(Sys.time() - start_t, units = "secs")
        save(results_BIDAG_3, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_BiDAG_GES.Rdata", seed, k, sparsity_type)))
        
        start_t <- Sys.time()
        results_BIDAG_3_noplus <- orderMCMC(score_par_sim, startspace = gesCPDAG, MAP = FALSE, plus1 = FALSE, 
                                            chainout = TRUE, iterations = Bs_2, stepsave = thinning_constant, 
                                            hardlimit = max_spars_ges-1)
        time_BIDAG3_noplus <- as.numeric(Sys.time() - start_t, units = "secs")
        save(results_BIDAG_3_noplus, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_BiDAG_GES_noplus.Rdata", 
                                                                         seed, k, sparsity_type)))
        
        # 5. BDGraph
        time_bdmcmc <- NA
        if (bdmcmc_run) {
          start_t <- Sys.time()
          results_bdmcmc <- bdgraph(scale(data, scale=F), method="ggm", algorithm="bdmcmc", 
                                    iter = Bs_2 + floor(burn_in_rate*Bs_2), burnin = floor(burn_in_rate*Bs_2))
          time_bdmcmc <- as.numeric(Sys.time() - start_t, units = "secs")
          save(results_bdmcmc, file = file.path(full_path, sprintf("seed_%d_rep_%d_%s_bdmcmc.Rdata", seed, k, sparsity_type)))
        }
        
        # Metric Calculations
        results_list <- list()
        
        # Iterations to Use in Aggregations
        total_saved <- floor(Bs_2 / thinning_constant)
        total_saved_BIDAG <- length(results_BIDAG$traceadd$incidence)
        idx_bidag <- ceiling(burn_in_rate * total_saved_BIDAG):total_saved_BIDAG
        idx_brood <- 1:total_saved
        
        # Helper Function:
        # Takes a predicted adjacency matrix and compares it to the truth across a few metrics
        calc_metrics <- function(pred_adj, true_adj, type_str, time_val, is_skel = FALSE) {
          # Vectorize
          if (is_skel) {
            # Use lower triangle for skeletons to avoid redundant symmetric pairs
            p_vec <- as.numeric(pred_adj[lower.tri(pred_adj)])
            t_vec <- as.numeric(true_adj[lower.tri(true_adj)])
          } else {
            p_vec <- as.numeric(pred_adj)
            t_vec <- as.numeric(true_adj)
          }
          
          t_f <- factor(t_vec, levels = c(1, 0))
          p_f <- factor(1 * (p_vec > 0.5), levels = c(1, 0))
          
          tibble(
            type       = type_str,
            AUC_ROC    = Rfast::auc(t_vec, p_vec),
            AUC_PR     = yardstick::pr_auc_vec(t_f, p_vec),
            F1         = yardstick::f_meas_vec(t_f, p_f),
            Pr_plus    = sum(p_vec[t_vec == 1]) / max(1, sum(t_vec == 1)),
            Pr_minus   = sum(p_vec[t_vec == 0]) / max(1, sum(t_vec == 0)),
            Time       = as.numeric(time_val)
          )
        }
        
        # Setup of true values
        trueDAG <- 1 * (trueDAGedges_curr != 0)
        trueDAG_skel <- 1 * ((trueDAG + t(trueDAG)) > 0)
        
        # Helper function:
        # Calculates estimated edge probabilities for BiDAG runs
        edge_probs_bidag <- function(trace, steps, trans = FALSE) {
          m <- Reduce('+', trace[steps]) / length(steps)
          if(trans) m <- t(m)
          return(as.matrix(m))
        }
        
        edge_probs_brood <- function(graph_list, steps) {
          # BROOD needs transpose to align with truth
          t(as.matrix(Reduce('+', graph_list[steps]))) / length(steps)
        }
        
        # BiDAG Results (PC, MAP, GES)
        bidag_configs <- list(
          list(res = results_BIDAG, lab = "PC", t = time_BIDAG),
          list(res = results_BIDAG_2, lab = "MAP", t = time_BIDAG2),
          list(res = results_BIDAG_3, lab = "GES", t = time_BIDAG3),
          list(res = results_BIDAG_noplus, lab = "PC_noplus", t = time_BIDAG_noplus),
          list(res = results_BIDAG_2_noplus, lab = "MAP_noplus", t = time_BIDAG2_noplus),
          list(res = results_BIDAG_3_noplus, lab = "GES_noplus", t = time_BIDAG3_noplus)
        )
        
        for (conf in bidag_configs) {
          m_adj <- edge_probs_bidag(conf$res$traceadd$incidence, idx_bidag)
          # Directed
          results_list[[paste0(conf$lab)]] <- calc_metrics(m_adj, trueDAG, conf$lab, conf$t)
          # Skeleton
          m_skel <- m_adj + t(m_adj)
          results_list[[paste0(conf$lab, "_skel")]] <- calc_metrics(m_skel, trueDAG_skel, paste0(conf$lab, "_skel"), 
                                                                    conf$t, is_skel = TRUE)
        }
        
        # BROOD Results (sim_PC, sim_MAP, sim_GES)
        brood_configs <- list(
          list(res = results,   lab = "sim_PC",  t = time_BROOD),
          list(res = results_2, lab = "sim_MAP", t = time_BROOD2),
          list(res = results_3, lab = "sim_GES", t = time_BROOD3)
        )
        
        for (conf in brood_configs) {
          # Main Graph
          m_adj <- edge_probs_brood(conf$res$graphs, idx_brood)
          results_list[[paste0(conf$lab)]] <- calc_metrics(m_adj, trueDAG, conf$lab, conf$t)
          
          # Main Graph Skeleton
          m_skel <- m_adj + t(m_adj)
          results_list[[paste0(conf$lab, "_skel")]] <- calc_metrics(m_skel, trueDAG_skel, paste0(conf$lab, "_skel"), 
                                                                    conf$t, is_skel = TRUE)
          
          # Plus Graph
          m_plus <- edge_probs_brood(conf$res$plusgraphs, idx_brood)
          results_list[[paste0(conf$lab, "_plusgraph")]] <- calc_metrics(m_plus, trueDAG, paste0(conf$lab, "_plusgraph"), 
                                                                         conf$t)
          
          # Plus Graph Skeleton
          m_plus_skel <- m_plus + t(m_plus)
          results_list[[paste0(conf$lab, "_plusgraph_skel")]] <- calc_metrics(m_plus_skel, trueDAG_skel, 
                                                                              paste0(conf$lab, "_plusgraph_skel"), 
                                                                              conf$t, is_skel = TRUE)
        }
        
        if (bdmcmc_run) {
          m_bdmcmc_skel <- results_bdmcmc$p_links + t(results_bdmcmc$p_links)
          results_list[["bdmcmc_skel"]] <- calc_metrics(m_bdmcmc_skel, trueDAG_skel, "bdmcmc_skel", time_bdmcmc, is_skel = TRUE)
        }
        
        output_full <- bind_rows(results_list) |>
          add_column(n_nodes = n, n_obs = N, graph_model = method, 
                     data_model = datatype, pc_thresh = thresh, 
                     sparsity = sparsity_type, seed = seed, replicate = k, .before = 1)
        save(output_full, file=file_summary)
        print("saved summary")
      }
      else{
        load(file_summary)
      }
      output_full
    }
  return(output)
}
        