library(tidyverse)
library(BiDAG)
library(matrixStats)
library(gtools)
library(Rfast)
library(pcalg)
library(RBGL)
library(foreach)
library(doParallel)

source("./Functions.R")

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



seeds <- 101:200
replicates <- 10
pc_alg_thresh <- c(0.05, 0.5, 0.95)
graph_model <- c("ERs", "iER")
nodes <- c(20, 200)
dataset_multiplier <- c(2, 10)
data_model <- c("Gaussian", "Log_Gaussian")
Bs <- 10000


cl <- makeCluster(detectCores()-1, 'PSOCK')
clusterExport(cl, varlist=c("Bs", "seeds", "pc_alg_thresh", "graph_model", "nodes",
                            "dataset_multiplier", "data_model", "replicates"))
clusterExport(cl, varlist=c("graph_mcmc", "rmvDAG","create_weights", 
                            "shrink_search_space_v2", "calculate_set_size", "expand_search_space", 
                            "sample_graph", "sample_plus_graph", "sample_from_2_orders",
                            "sample_from_multiple_orders", "implement_order_v1",
                            "implement_order_v2", "dettwobytwo", "bge_score_node",
                            "score_full_space", "score_plus_space", "create_banned_parent_table",
                            "create_banned_plus_parent_table", "parents_mapping",
                            "create_parent_table_idx", "create_parent_table",
                            "banned_parents_mapping", "plus_parents_mapping",
                            "mcmc_sampler_step", "powerset", "index_finder_plus",
                            "implement_order_random", "sample_graph_random",
                            "sample_minus_graph", "sample_from_2_graphs",
                            "rmvlogDAG", "shrink_search_space", "calculate_birth_rate",
                            "calculate_death_rate"),
              envir = .GlobalEnv)


clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(BiDAG))
clusterEvalQ(cl, library(matrixStats))
clusterEvalQ(cl, library(gtools))
clusterEvalQ(cl, library(Rfast))
clusterEvalQ(cl, library(pcalg))
clusterEvalQ(cl, library(RBGL))
clusterEvalQ(cl, library(pROC))

registerDoParallel(cl)

op <- foreach(method=graph_model, .errorhandling='pass', .combine='rbind') %:%
  foreach(thresh=pc_alg_thresh, .combine='rbind') %:%
    foreach(n=nodes, .combine='rbind') %:%
      foreach(m=dataset_multiplier, .combine='rbind') %:%
        foreach(datatype=data_model, .combine='rbind') %:%
          foreach(i=seeds, .combine='rbind') %:%
            foreach(k=replicates, .combine='rbind') %dopar% {
              N <- m*n
              file_name1 <- paste0("_method_", method, "_n_", n, "_N_", N, "_seed_", i, 
                                   "_pcthresh_", thresh, "_data_model_", datatype, "_sim.Rdata")
              set.seed(i)
              if(method=="ERs"){
                trueDAGedges <- as(pcalg::randDAG(n = n, d = 4, 
                                                  wFUN = list(runif, min=0.4, max=2)), "matrix")
              }
              else{
                trueDAGedges <- as(pcalg::randDAG(n = n, d = 4, 
                                                  method = "interEr", par1 = 2, par2 = 0.1), "matrix")
              }
              trueDAG <- 1*(trueDAGedges != 0)
              trueCPDAG <- BiDAG:::dagadj2cpadj(trueDAG)
              
              if(datatype=="Gaussian"){
                data <- rmvDAG(trueDAGedges, N, standardise = FALSE)
              }
              else{
                data <- rmvlogDAG(trueDAGedges, N, standardise = FALSE)
              }
              
              score_par_sim <- scoreparameters("bge", data)
              
              cor_mat <- cor(data)
              pc_fit <- skeleton(suffStat = list(C = cor_mat, n = nrow(data)),
                                 indepTest = gaussCItest, ## indep.test: partial correlations
                                 alpha = 0.4, labels = paste0("V", 1:ncol(data)),
                                 method="stable", verbose = FALSE)
              
              space_PC <- 1*as(pc_fit@graph, "matrix")
              
              results <- graph_mcmc(space_PC, score_par_sim, B=Bs, 
                                    verbose=FALSE, bounce = 0.000000001)
              save(results, file=paste("./Output/Simulation_3", file_name1, sep="/"))
              
              
              
              levels <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
              output_full <- vector(mode="numeric", length=4*length(levels))
              for(ls in 1:length(levels)){
                mean_edge_est <- t(apply(results$graphs[,,1000:5000], c(1,2), mean))
                graph_est <- 1*(mean_edge_est>levels[ls])
                output <- as.numeric(table(graph_est, trueDAG))
                output_full[(ls-1)*4+1:ls*4] <- output
              }
              output_full
            }