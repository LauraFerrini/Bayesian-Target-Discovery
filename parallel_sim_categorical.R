# simulations over samples

## q : number of nodes/variables
## pr : prior probability of edge inclusion
## nk : number of observations for experimental context

## size_target : number of intervened nodes for each experimental context


nc   = detectCores()   # detect cores
clus = makeCluster(min(nc,maxcore)) # create cluster

clusterExport(clus, ls(), envir = .GlobalEnv) # export vars to nodes

out = parLapply(clus, 1:N, function(i){
  
  library(BCDAG)
  library(gtools)
  library(pcalg)
  
  n_all = rep(nk, K)
  
  ## Generate intervention targets I.cal, i.e. fix nodes that are intervened in the experimental context and
  
  ## Generate the data (guarantee that each variable has at least two different levels)
  
  file_to_load = paste0("data_sim_K4_size",size_target,".RData")
  
  load(file_to_load)
  
  X = datasets$X[[which(n_set == nk)]][[i]][1:sum(n_all),]
  
  dag = datasets$dag[[i]]
  
  I.cal = datasets$I[[i]][,1:K]
  
  X.lev = sapply(1:ncol(X), function(j) length(unique(X[,j])))
  
  #############################
  ## [1] Our Bayesian method ##
  #############################
  
  source("mcmc_dags_targets_constrained.R")
  
  set.seed(123)
  out = mcmc_dag_targets(X, S, burn, ne = NULL, a = a, a_k = a_k, b_k = b_k, a_dag, b_dag, n_all = n_all, X.lev = X.lev)
  
  ## Posterior probabilities of intervention
  
  I.prob = apply(X = out$Targets_post, FUN = mean, MARGIN = c(1,2))
  
  ## Posterior probabilities of edge inclusion
  
  PPI = apply(X = out$DAG_post, FUN = mean, MARGIN = c(1,2))
  
  
  #####################################
  ## [2] He & Geng Algorithm 1 and 2 ##
  #####################################
  
  source("He_Geng_method_categorical.R")
  
  index = rep(1:K, n_all)
  alpha = 0.01
  
  set.seed(123)
  out_HG_1 = HG_alg_1(X, index, alpha)
  
  set.seed(123)
  out_HG_2 = HG_alg_2(X, index, alpha, S = 100, s_dim = round(K/2))
  
  
  #########################################################################
  ## [3] Oracle version of our Bayesian method for target identification ##
  #########################################################################
  
  source("mcmc_dags_targets_constrained_oracle.R")
  
  set.seed(123)
  out_oracle = mcmc_dag_targets_oracle(X, S, true_dag = dag, burn, ne = NULL, a = NULL, a_k = a_k, b_k = b_k, a_dag, b_dag, n_all = n_all, X.lev = X.lev, I_bar = NULL)
  
  I.prob.oracle = apply(X = out_oracle$Targets_post, FUN = mean, MARGIN = c(1,2))
  
  
  ##################################################################
  ## [4] Oracle version of our Bayesian method for DAG estimation ##
  ##################################################################
  
  source("mcmc_dags_targets_constrained_oracle_targets.R")
  
  set.seed(123)
  out_oracle = mcmc_dag_targets_oracle_targets(X, S, true_targets = I.cal, burn, ne = NULL, a = NULL, a_k = a_k, b_k = b_k, a_dag, b_dag, n_all = n_all, X.lev = X.lev, I_bar = NULL)
  
  PPI.oracle = apply(X = out_oracle$DAG_post, FUN = mean, MARGIN = c(1,2))
  
  
  ###################
  ## Return output ##
  ###################
  
  return(out = list(dag   = dag,
                    I.cal = I.cal,
                    X   = X,
                    index = index,
                    I.prob = I.prob,
                    PPI    = PPI,
                    I.prob.oracle = I.prob.oracle,
                    PPI.oracle = PPI.oracle,
                    HG_1 = out_HG_1,
                    HG_2 = out_HG_2))
  
})
