library(tidyverse)
library(BiDAG)
library(matrixStats)
library(gtools)
library(Rfast)
library(pcalg)
library(RBGL)

source("./Functions.R")

set.seed(101)
trueDAGedges <- as(pcalg::randDAG(n = 20, d = 4, 
                                  wFUN = list(runif, min=0.4, max=2)), "matrix")
trueDAG <- 1*(trueDAGedges != 0)
trueCPDAG <- BiDAG:::dagadj2cpadj(trueDAG)

### This function generates Gaussian data from a DAG
# following the topological order

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



data <- rmvDAG(trueDAGedges, 200, standardise = FALSE)

score_par_test <- scoreparameters("bge", data)

cor_mat <- cor(data)

pc_fit <- pc(suffStat = list(C = cor_mat, n = nrow(data)),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha = 0.4, labels = paste0("V", 1:ncol(data)),
             skel.method="stable", verbose = FALSE)


space_PC <- 1*as(pc_fit@graph, "matrix")

probDAGs <- iterativeMCMC(score_par_test, scoreout = TRUE, 
                          hardlimit = 16, softlimit = 10,
                          MAP = FALSE, posterior=0.2, alpha=0.4)

mcmc_run <- graph_mcmc(probDAGs$startspace, score_par_test, B=2000)
