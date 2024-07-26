library(tidyverse)
library(BiDAG)
library(matrixStats)
library(gtools)
library(Rfast)
library(pcalg)
library(RBGL)
library(foreach)
library(doParallel)
library(igraph)

source("./Functions.R")

sim_root <- "./Output/Simulation_5"

theme_tuned <- ggplot2::theme(legend.position="bottom",
                              legend.direction="horizontal",
                              legend.key = element_rect(fill="white"),
                              legend.box="vertical", 
                              legend.margin=margin(),
                              text = element_text(family="serif", size = 11, color = "black"),
                              panel.background = element_blank(),
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(),
                              plot.title = element_text(family="serif", size = 12, face = "bold", 
                                                        hjust=0.5, margin = margin(b = 20, r=0, l=0, t=20)),
                              plot.caption = element_text(family="serif", hjust = 0),
                              strip.background = element_blank(), 
                              axis.line = element_line(color="black", size = .25),
                              axis.ticks.x  = element_blank(),
                              axis.text.x=element_text (color = "black"),
                              axis.text.y=element_text (color = "black"))


rmvDAG <- function(trueDAGedges, N, standardise = TRUE) {
  trueDAG <- 1*(trueDAGedges != 0) # the edge presence in the DAG
  n <- ncol(trueDAG) # number of variables
  data <- matrix(0, nrow = N, ncol = n) # to store the simulated data
  top_order <- rev(BiDAG:::DAGtopartition(n, trueDAG)$permy) # go down order
  for (jj in top_order) {
    parents <- which(trueDAG[, jj] == 1) # find parents
    lp <- length(parents) # number of parents
    if (lp == 0) { # no parents
      data[, jj] <- 0
    } else if (lp == 1) { # one parent
      data[, jj] <- data[, parents]*trueDAGedges[parents, jj]
    } else { # more than one parent
      data[, jj] <- colSums(t(data[, parents])*trueDAGedges[parents, jj])
    }
    # add random noise
    data[, jj] <- data[, jj] + rnorm(N)
  }
  if(standardise) { # whether to standardise
    scale(data)
  } else {
    data
  }
}


rmvlogDAG <- function(trueDAGedges, N, standardise = TRUE) {
  trueDAG <- 1*(trueDAGedges != 0) # the edge presence in the DAG
  n <- ncol(trueDAG) # number of variables
  data <- matrix(0, nrow = N, ncol = n) # to store the simulated data
  top_order <- rev(BiDAG:::DAGtopartition(n, trueDAG)$permy) # go down order
  for (jj in top_order) {
    parents <- which(trueDAG[, jj] == 1) # find parents
    lp <- length(parents) # number of parents
    if (lp == 0) { # no parents
      data[, jj] <- 0
    } else if (lp == 1) { # one parent
      data[, jj] <- data[, parents]*trueDAGedges[parents, jj]
    } else { # more than one parent
      data[, jj] <- colSums(t(data[, parents])*trueDAGedges[parents, jj])
    }
    # add random noise
    data[, jj] <- data[, jj] + rnorm(N)
  }
  data <- exp(data)
  if(standardise) { # whether to standardise
    scale(data)
  } else {
    data
  }
}


rmvlogexpDAG <- function(trueDAGedges, N, standardise = TRUE) {
  trueDAG <- 1*(trueDAGedges != 0) # the edge presence in the DAG
  n <- ncol(trueDAG) # number of variables
  data <- matrix(0, nrow = N, ncol = n) # to store the simulated data
  top_order <- rev(BiDAG:::DAGtopartition(n, trueDAG)$permy) # go down order
  for (jj in top_order) {
    parents <- which(trueDAG[, jj] == 1) # find parents
    lp <- length(parents) # number of parents
    if (lp == 0) { # no parents
      data[, jj] <- 0
    } else if (lp == 1) { # one parent
      data[, jj] <- log(data[, parents]*trueDAGedges[parents, jj])
    } else { # more than one parent
      data[, jj] <- log(colSums(t(data[, parents])*trueDAGedges[parents, jj]))
    }
    # add random noise
    data[, jj] <- data[, jj] + rnorm(N)
  }
  data <- exp(data)
  if(standardise) { # whether to standardise
    scale(data)
  } else {
    data
  }
}

g2Q <- function(g, sparse = FALSE) {
  Q <- get.adjacency(g, sparse = sparse)
  perm <- sample.int(nrow(Q))
  Q2 <- Q 
  Q2[perm, perm] <- Q * upper.tri(Q)
  Q2
}

# g2Q <- function(g, sparse = FALSE) {
#   Q <- get.adjacency(g, sparse = sparse)
#   perm <- sample.int(nrow(Q))
#   Q <- Q[perm, perm]
#   Q * upper.tri(Q)
# }


seeds <- 101:150
replicates <- 1:5
# pc_alg_thresh <- c(0.05, 0.5, 0.95)
pc_alg_thresh <- c(0.95, 0.05)
graph_model <- c("ERs", "iER", "SBM_2pos", "SBM_2neg", "SBM_2both")
#nodes <- c(20, 200)
nodes <- 20
dataset_multiplier <- c(2, 10)
# data_model <- c("Gaussian", "Log_Gaussian", "LogExp_Gaussian")
data_model <- c("Gaussian")
Bs <- 20000


cl <- makeCluster(detectCores()-1, 'PSOCK')
clusterExport(cl, varlist=c("Bs", "seeds", "pc_alg_thresh", "graph_model", "nodes",
                            "dataset_multiplier", "data_model", "replicates", "sim_root"))
clusterExport(cl, varlist=c("graph_mcmc", "rmvDAG", "g2Q", "create_weights", 
                            "shrink_search_space_v2", "calculate_set_size", "expand_search_space", 
                            "sample_graph", "sample_plus_graph", "sample_from_2_orders",
                            "sample_from_multiple_orders", "implement_order_v1",
                            "implement_order_v2", "dettwobytwo", "bge_score_node",
                            "logMinusExp", "bge_score_plus_parent", "score_plus_space_new",
                            "create_banned_plus_parent_table_new", "parents_mapping",
                            "create_parent_table_idx", "create_parent_table",
                            "banned_parents_mapping", "plus_parents_mapping",
                            "mcmc_sampler_step", "powerset", "index_finder_plus",
                            "implement_order_random", "sample_graph_random",
                            "sample_minus_graph", "sample_from_2_graphs",
                            "rmvlogDAG", "shrink_search_space", "calculate_birth_rate",
                            "calculate_death_rate", "rmvlogexpDAG", "are_equivalent"),
              envir = .GlobalEnv)


clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(BiDAG))
clusterEvalQ(cl, library(matrixStats))
clusterEvalQ(cl, library(gtools))
clusterEvalQ(cl, library(Rfast))
clusterEvalQ(cl, library(pcalg))
clusterEvalQ(cl, library(RBGL))
clusterEvalQ(cl, library(pROC))
clusterEvalQ(cl, library(igraph))

registerDoParallel(cl)

op <- foreach(method=graph_model, .errorhandling='pass', .combine='rbind') %:%
  foreach(thresh=pc_alg_thresh, .combine='rbind') %:%
    foreach(n=nodes, .combine='rbind') %:%
      foreach(m=dataset_multiplier, .combine='rbind') %:%
        foreach(datatype=data_model, .combine='rbind') %:%
          foreach(i=seeds, .combine='rbind') %:%
            foreach(k=replicates, .combine='rbind') %dopar% {
              Bs_2 <- max(Bs, round(n*n*log(n)/2))
              N <- m*n
              folder_name1 <- paste0("_method_", method, "_n_", n, "_N_", N)
              folder_name2 <- paste0("_data_model_", datatype, "_pcthresh_", thresh)
              
              if(!dir.exists(paste(sim_root, folder_name1, sep="/"))){
                dir.create(paste(sim_root, folder_name1, sep="/"))
              }
              if(!dir.exists(paste(sim_root, folder_name1, folder_name2, sep="/"))){
                dir.create(paste(sim_root, folder_name1, folder_name2, sep="/"))
              }
              file_name1 <- paste0("_seed_", i, "_rep_", k, "_sim_PC.Rdata")
              file_name2 <- paste0("_seed_", i, "_rep_", k, "_summary.Rdata")
              
              curr_files <- dir(paste(sim_root, folder_name1, folder_name2, sep="/"))
              
              if(!(file_name2 %in% curr_files)){
                file_name3 <- paste0("seed_", i, "_rep_", k, "_graph.Rdata")
                file_name4 <- paste0("seed_", i, "_rep_", k, "_MAP.Rdata")
                file_name5 <- paste0("seed_", i, "_rep_", k, "_PC.Rdata")
                file_name6 <- paste0("seed_", i, "_rep_", k, "_sim_MAP.Rdata")
                file_name7 <- paste0("seed_", i, "_rep_", k, "_BiDAG_PC.Rdata")
                file_name8 <- paste0("seed_", i, "_rep_", k, "_BiDAG_MAP.Rdata")
                set.seed(i)
                if(method=="ERs"){
                  trueDAGedges <- as(pcalg::randDAG(n = n, d = 4, 
                                                    wFUN = list(runif, min=0.4, max=2)), "matrix")
                }
                else if(method=="iER"){
                  trueDAGedges <- as(pcalg::randDAG(n = n, d = 4, 
                                                    method = "interEr", par1 = 2, par2 = 0.1), "matrix")
                }
                else if(method=="SBM_2pos"){
                  dagedges <- sample_sbm(n=n, pref.matrix=matrix(c(0.15, 0.05, 0.05, 0.08), 
                                                                 nrow=2, ncol=2),
                                         block.sizes=c(floor(0.8*n), ceiling(0.2*n)),
                                         directed=TRUE)
                  dagweights <- runif(n*n, min=0.4, max=2)
                  trueDAGedges <- g2Q(dagedges) * dagweights
                }
                else if(method=="SBM_2neg"){
                  dagedges <- sample_sbm(n=n, pref.matrix=matrix(c(0.15, 0.05, 0.05, 0.08), 
                                                                 nrow=2, ncol=2),
                                         block.sizes=c(floor(0.8*n), ceiling(0.2*n)),
                                         directed=TRUE)
                  dagweights <- runif(n*n, min=-2, max=0.4)
                  trueDAGedges <- g2Q(dagedges) * dagweights
                }
                else{
                  dagedges <- sample_sbm(n=n, pref.matrix=matrix(c(0.15, 0.05, 0.05, 0.08), 
                                                                 nrow=2, ncol=2),
                                         block.sizes=c(floor(0.8*n), ceiling(0.2*n)),
                                         directed=TRUE)
                  dagweights <- runif(n*n, min=-2, max=2)
                  trueDAGedges <- g2Q(dagedges) * dagweights
                }
                trueDAG <- 1*(trueDAGedges != 0)
                trueCPDAG <- BiDAG:::dagadj2cpadj(trueDAG)
                
                if(datatype=="Gaussian"){
                  data <- rmvDAG(trueDAGedges, N, standardise = FALSE)
                }
                else if(datatype=="Log_Gaussian"){
                  data <- rmvlogDAG(trueDAGedges, N, standardise = FALSE)
                }
                else{
                  data <- rmvlogexpDAG(trueDAGedges, N, standardise = FALSE)
                }
                
                score_par_sim <- scoreparameters("bge", data)
                
                cor_mat <- cor(data)
                pc_fit <- skeleton(suffStat = list(C = cor_mat, n = nrow(data)),
                                   indepTest = gaussCItest, ## indep.test: partial correlations
                                   alpha = thresh, labels = paste0("V", 1:ncol(data)),
                                   method="stable", verbose = FALSE)
                
                space_PC <- 1*as(pc_fit@graph, "matrix")
                
                save(space_PC, file=paste(sim_root, folder_name1, folder_name2,
                                          file_name5, sep="/"))
                
                set.seed(i*10+k)
                save(trueDAGedges, file=paste(sim_root, folder_name1, folder_name2,
                                              file_name3, sep="/"))
                
                results <- graph_mcmc(space_PC, score_par_sim, B=Bs_2, 
                                      verbose=FALSE)
                
                save(results, file=paste(sim_root, folder_name1, folder_name2, 
                                         file_name1, sep="/"))
                
                results_BIDAG <- orderMCMC(score_par_sim, startspace = space_PC, 
                                           MAP = FALSE, plus1 = TRUE, chainout = TRUE, 
                                           startorder = 1:n,
                                           iterations = Bs_2,
                                           hardlimit = max(colSums(space_PC)))
                
                save(results_BIDAG, file=paste(sim_root, folder_name1, folder_name2, 
                                               file_name7, sep="/"))
                
                default_alpha <- min(0.4, 20/n)/2
                
                bestDAGs <- iterativeMCMC(score_par_sim, scoreout = TRUE, hardlimit = 16, 
                                          softlimit = 10, alpha = default_alpha)
                
                
                save(bestDAGs, file=paste(sim_root, folder_name1, folder_name2,
                                          file_name4, sep="/"))
                
                results_2 <- graph_mcmc(t(bestDAGs$endspace), score_par_sim, B=Bs_2, 
                                        verbose=FALSE)
                
                
                results_BIDAG_2 <- orderMCMC(score_par_sim, startspace = bestDAGs$endspace, 
                                             MAP = FALSE, plus1 = TRUE, chainout = TRUE, 
                                             startorder = 1:n,
                                             iterations = Bs_2,
                                             hardlimit = max(colSums(bestDAGs$endspace)))
                
                
                save(results_2, file=paste(sim_root, folder_name1, folder_name2, 
                                           file_name6, sep="/"))
                save(results_BIDAG_2, file=paste(sim_root, folder_name1, folder_name2, 
                                                 file_name8, sep="/"))
                
                
                # levels <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
                # output_full <- vector(mode="numeric", length=4*length(levels))
                burn_in <- ceiling(Bs_2/3)
                trueDAG_skel <- (trueDAG + t(trueDAG) > 0)*1
                trueDAG_skel_lower <- trueDAG_skel[lower.tri(trueDAG_skel)]
                truegraph_obj <- m2graph(trueDAG)
                mean_edge_est <- t(apply(results$skeletons[,,burn_in:Bs_2], c(1,2), mean))
                w_mean_edge_est <- t(apply(results$skeletons[,,burn_in:Bs_2], c(1,2), 
                                           weighted.mean, w=results$weight[burn_in:Bs]))
                auc_orig <- pROC::auc(pROC::roc(trueDAG_skel_lower,
                                                mean_edge_est[lower.tri(mean_edge_est)]))
                auc_weighted <- pROC::auc(pROC::roc(trueDAG_skel_lower,
                                                    w_mean_edge_est[lower.tri(w_mean_edge_est)]))
                # pct_equiv <- sum(sapply(burn_in:Bs, function(i){
                #   are_equivalent(trueDAG, t(results$graphs[,,i]))}))/(Bs-burn_in+1)
                # shd_totals <- sapply(seq(burn_in, Bs_2, n), function(i){
                #   shd(truegraph_obj, m2graph(t(results$graphs[,,i])))
                # })
                # mean_edge_est2 <- t(apply(results_2$skeletons[,,burn_in:Bs_2], c(1,2), mean))
                # w_mean_edge_est2 <- t(apply(results_2$skeletons[,,burn_in:Bs_2], c(1,2), 
                #                             weighted.mean, w=results_2$weight[burn_in:Bs_2]))
                # auc_orig2 <- pROC::auc(pROC::roc(trueDAG_skel_lower,
                #                           mean_edge_est2[lower.tri(mean_edge_est2)]))
                # auc_weighted2 <- pROC::auc(pROC::roc(trueDAG_skel_lower,
                #                               w_mean_edge_est2[lower.tri(w_mean_edge_est2)]))
                # shd_totals2 <- sapply(seq(burn_in, Bs_2, n), function(i){
                #   shd(truegraph_obj, m2graph(t(results_2$graphs[,,i])))
                # })
                # output_full <- c(auc_orig, auc_weighted, shd_totals)
                output_full <- c(auc_orig, auc_weighted)
                save(output_full, file=paste(sim_root, folder_name1, folder_name2,
                                             file_name2, sep="/"))
              }
              else{
                output_full <- c(NULL, NULL)
              }
              output_full
              
            }