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

rmvDAG2 <- function(trueDAGedges, N, standardise = TRUE, dft = 4) {
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
    data[, jj] <- data[, jj] + rt(N, dft)
  }
  if(standardise) { # whether to standardise
    scale(data)
  } else {
    data
  }
}

rmvDAG3 <- function(trueDAGedges, N, standardise = TRUE, sd2=2) {
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
    mixture_component <- sample(c(TRUE, FALSE), N, replace=TRUE)
    noise1 <- rnorm(N)
    noise2 <- rnorm(N, sd=sd2)
    data[,jj] <- ifelse(mixture_component, data[,jj]+noise2, data[,jj]+noise1)
    if(standardise) { # whether to standardise
      scale(data)
    } else {
      data
    }
  }
  return(data)
}

g2Q <- function(g, sparse = FALSE) {
  Q <- as_adjacency_matrix(g, sparse = sparse)
  perm <- sample.int(nrow(Q))
  Q2 <- Q 
  Q2[perm, perm] <- Q * upper.tri(Q)
  Q2
}

spectral_rescale <- function(W, tolerance=0.001){
  eigvals <- eigen(W, only.values=TRUE)$values
  s <- max(Mod(eigvals))  # spectral radius
  if(s >= 1) {
    W <- W / (s + tolerance)  # shrink so rho(W) < 1
  }
  return(W)
}


graphnode_edge_prior_unnormalized <- function(k, p, q, log=TRUE){
  coef_degree <- pnbinom(p-1-k, size=k+1, prob=q, log.p=log)
  if(log){
    return(-log(q) -log(q) + coef_degree)
  }
  return(1/(p*q)*coef_degree)
}

graphnode_edge_prior <- function(k, p, q, log=TRUE){
  node_normalizing_c <- sapply(0:(p-1), function(i){graphnode_edge_prior_unnormalized(i, p, q, log)})
  if(log){
    return(graphnode_edge_prior_unnormalized(k, p, q, log)-logSumExp(node_normalizing_c))
  }
  return(graphnode_edge_prior_unnormalized(k, p, q, log)/sum(node_normalizing_c))
}

graph_edge_prior <- function(G, q, log=TRUE){
  p <- nrow(G)
  k_vec <- colSums(G)
  
  if(log){
    return(sum(sapply(k_vec, function(i){graphnode_edge_prior(i, p, q, log)})))
  }
  return(prod(sapply(k_vec, function(i){graphnode_edge_prior(i, p, q, log)})))
}

usrscoreparameters <- function(initparam, usrpar=list(pctesttype="usrCItest", edgepf=NULL, 
                                                      chi=0.5, delta=NULL, eta=NULL)){
  
  initparam$UN <- initparam$U0 + crossprod(initparam$data)
  initparam$alpha_post <- initparam$alpha_vec + initparam$N
  
  p <- initparam$n
  N <- initparam$N
  scoreconstlist <- vector(mode="list", length=p)
  for(i in 0:(p-1)){
    a_i_vec <- initparam$alpha_vec-i
    scoreconstlist[[i+1]] <- -N/2*log(pi)-lgamma(a_i_vec/2-1)+lgamma((a_i_vec+N)/2-1)
  }
  
  initparam$scoreconstlist <- scoreconstlist
  if(!is.null(usrpar$edgepf)){
    edgepf <- usrpar$edgepf
    if(edgepf < 1 & edgepf >= 0){
      logedgepf <- numeric(p)
      for(j in 0:(p-1)){
        logedgepf[j+1] <- graphnode_edge_prior(j, p, edgepf)
      }
    }
    else{
      stop("Need for the edge penalty factor to be in range [0,1)")
    }
    initparam$logedgepvec <- logedgepf
  }
  
  initparam
}



usrDAGcorescore <- function(j, parentnodes, n, param) {
  ll <- dagwishart_score_node(j = j, parentnodes=parentnodes,
                              N=n, param=param)
  return(ll)
}

