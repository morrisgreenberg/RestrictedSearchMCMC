library(tidyverse)
library(BiDAG)
library(matrixStats)
library(gtools)
library(Rfast)
library(pcalg)
library(RBGL)

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

mcmc_run <- graph_mcmc(probDAGs$startspace, score_par_test, B=5000)


Bs <- 5000
dater <- mcmc_run

total_size <- sapply(1:Bs, function(i){sum(t(dater$spaces[,,i]))})
true_edges_space <- sapply(1:Bs, function(i){sum(which(trueDAGedges>0) %in% 
                                                   which(t(dater$spaces[,,i])==1))})
false_edges_space <- sapply(1:Bs, function(i){sum(!(which(t(dater$spaces[,,i])==1) %in% 
                                                      which(trueDAGedges>0)))})

changed_space_idx <- which(sapply(1:(Bs-1), function(i){
  return(max(t(dater$spaces[,,i])!=t(dater$spaces[,,i+1])))
})==1)

changed_space <- rep(0, Bs)
changed_space[changed_space_idx+1] <- 1

num_shrunk <- sapply(1:(Bs-1), function(i){
  return(sum(!(which(t(dater$spaces[,,i])==1) %in% which(t(dater$spaces[,,i+1])==1))))
})
num_shrunk <- c(0, num_shrunk)

num_expand <- sapply(1:(Bs-1), function(i){
  return(sum(!(which(t(dater$spaces[,,i+1])==1) %in% which(t(dater$spaces[,,i])==1))))
})
num_expand <- c(0, num_expand)


total_true_edges <- sum(trueDAGedges>0)

temp_df <- tibble(total_size, true_edges_space, false_edges_space,
                  changed_space, num_expand, num_shrunk) %>% 
  mutate(iter=row_number(),
         lag_ts = lag(total_size),
         lag_tes = lag(true_edges_space),
         lag_fes = lag(false_edges_space),
         ext_trues = total_true_edges-lag_tes,
         ext_falses = 380-total_true_edges-lag_fes,
         ext_edges = 380-lag_ts,
         expand_max = ifelse(num_expand < ext_trues, num_expand, ext_trues),
         shrink_max = ifelse(num_shrunk < lag_tes, num_shrunk, lag_tes))

expected_add_vec <- sapply(1:Bs, 
                           function(j){
                             ifelse(temp_df$num_expand[j]==0, 0, 
                                    sum(sapply(0:temp_df$expand_max[j], 
                                               function(i){
                                                 choose(temp_df$ext_falses[j],temp_df$num_expand[j]-i)*choose(temp_df$ext_trues[j],i)/choose(temp_df$ext_edges[j], temp_df$num_expand[j])*i
                                               })
                                    )
                             )
                           })

expected_shr_vec <- sapply(1:Bs, 
                           function(j){
                             ifelse(temp_df$num_shrunk[j]==0, 0, 
                                    sum(sapply(0:temp_df$shrink_max[j], 
                                               function(i){
                                                 choose(temp_df$lag_fes[j],temp_df$num_shrunk[j]-i)*choose(temp_df$lag_tes[j],i)/choose(temp_df$lag_ts[j], temp_df$num_shrunk[j])*i
                                               })
                                    )
                             )
                           })
full_sim_data <- temp_df %>%
  mutate(expected_add=expected_add_vec,
         expected_shr=expected_shr_vec,
         exp_true_edges=ifelse(changed_space==1, 
                               lag_tes + expected_add - expected_shr,
                               NA)) %>% 
  dplyr::select(iter, 
                "Total Edges in Space"=total_size, 
                "True Positive Edges"=true_edges_space, 
                "False Positive Edges"=false_edges_space,
                "Expected Positive Edges"=exp_true_edges)

full_sim_data_long <- full_sim_data %>%
  pivot_longer(-iter, names_to="edge_count_type", values_to="count")

p1 <- ggplot(full_sim_data_long, aes(x=iter, y=count, color=edge_count_type, 
                                     shape=edge_count_type))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept=sum(trueDAGedges>0), linetype="dashed", color="#33CCFF")+
  geom_text(x=600, y=56, color="#33CCFF", label="Total True Positive Edges")+
  labs(color="", shape="", x="Step", y="Count")+
  scale_color_manual(values=c("#333333", "#FF9933", "#FF3399", "#6633FF"))+
  scale_shape_manual(values=c("x", ".", ".", "."))+
  ggtitle(bquote("Space Size after"~.(Bs)~"Steps in the Markov Chain"))+
  theme_tuned


index_post_burn <- which(full_sim_data$`True Positive Edges`==sum(trueDAGedges>0))

index_post_burn <- ifelse(length(index_post_burn)>0, index_post_burn[1], floor(1/3*Bs))

mean_edge_est <- t(apply(dater$graphs[,,index_post_burn:Bs], c(1,2), mean))

library(pROC)
plot(roc(as.numeric(trueDAGedges > 0),as.numeric(mean_edge_est)))
auc(roc(as.numeric(trueDAGedges > 0),as.numeric(mean_edge_est)))
