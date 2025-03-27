mcmc_dag_targets = function(X, S, burn, ne = NULL, a = NULL, a_k = a_k, b_k = b_k, a_dag, b_dag, n_all = n_all, X.lev = X.lev, int_cont = NULL){
  
  ###########
  ## INPUT ##
  ###########
  
  # X : (n,q) data matrix
  
  # S    : number of MCMC iterations
  # burn : burn in period
  
  # w : prior probability of edge inclusion in p(D)
  # a : common hyper-parameter of the Dirichlet prior
  
  # a_k, b_k     : hyper-parameters of the Beta(a_k, b_k) prior on the probability of intervention for any node
  # a_dag, b_dag : hyper-parameters of the Beta(a_dag, b_dag) prior on the probability of edge inclusion
  
  # n_all : (K,1) vector with group sample sizes (n_1, ..., n_K)
  # X.lev : (q,1) vector with the number of levels of each categorical variable
  
  # int_cont : indexes (between 1 and K) of interventional datasets
  
  ############
  ## OUTPUT ##
  ############
  
  # X : (n,q) input data matrix
  
  # DAG_post    : (q,q,T) array collecting the T = (S - burn) adjacency matrices of the DAGs visited by the MCMC chain
  # Target_post : (q,K,T) array collecting the T = (S - burn) binary matrices representing the K intervention targets
  
  
  #########################
  ## Auxiliary functions ##
  #########################
  
  source("move_dag.R")
  source("marg_like_node.R")
  
  K = length(n_all)
  
  if(is.null(a)){a = 1}
  
  n = dim(X)[1]
  q = dim(X)[2]
  
  X[,1:ncol(X)] = lapply(X[,1:ncol(X)], as.factor)
  
  if(is.null(int_cont)){
    int_cont = 1:K
  }
  
  DAG_post     = array(NA, c(q, q, S))
  Targets_post = array(NA, c(q, K, S))
  
  num_edges = c()
  
  # Initialize the chain
  
  dag   = matrix(0, q, q)
  I.cal = matrix(0, q, K)
  
  D = lapply(1:K, function(k) matrix(0, n_all[k], q))
  
  for(k in 1:K){
    D[[k]][,which(I.cal[,k]==1)] = k
  }
  
  D = do.call(rbind, D)
  
  DAG_post[,,1] = dag
  
  out_D = lapply(1:K, function(k) matrix(0, n_all[k], q))
  
  for(k in 1:K){
    out_D[[k]][,which(I.cal[,k] == 1)] = k
  }
  
  cat("MCMC sampling")
  pb = utils::txtProgressBar(min = 2, max = S, style = 3)
  
  for(s in 1:S){
    
    
    ###################
    ## 1. Update DAG ##
    ###################
    
    ## 1.1 Propose new DAG
    
    dag_move = move(A = dag, ne = ne)
    
    dag_prop = dag_move$A_new
    nodes_prop = dag_move$nodes
    
    type.operator = dag_move$type.operator
    
    ## 1.2 Compute (log)priors p(dag), p(dag_prop) and (log)prior ratio
    
    logprior.star = lgamma(sum(dag_prop) + a_dag) +
      lgamma(q*(q-1)/2 - sum(dag_prop) + b_dag - 1)

    logprior = lgamma(sum(dag) + a_dag) +
      lgamma(q*(q-1)/2 - sum(dag) + b_dag - 1)

    logprior.ratio = logprior.star - logprior
    
    
    ## 1.3 Compute the marginal likelihood of Dag and dag_prop
    ## We can distinguish among the three types of operators (1,2,3) applied in the local move
    
    if(type.operator == 1 | type.operator == 2){
      
      ## (1,2) Insert or Delete a directed edge
      
      j = nodes_prop[2]
      
      marg_star = sum(sapply(unique(D[,j]), function(k) marg_j(j = j, dag = dag_prop, X.lev = X.lev, I = which(I.cal[,k] == 1), X = X[D[,j] == k,], a)))
      marg      = sum(sapply(unique(D[,j]), function(k) marg_j(j = j, dag = dag, X.lev = X.lev, I = which(I.cal[,k] == 1), X = X[D[,j] == k,], a)))
      
    }else{
      
      ## (3) Reverse a directed edge
      
      i = nodes_prop[1]
      j = nodes_prop[2]
      
      marg_star = sum(sapply(unique(D[,j]), function(k) marg_j(j = i, dag = dag_prop, X.lev = X.lev, I = which(I.cal[,k] == 1), X = X[D[,j] == k,], a))) +
                  sum(sapply(unique(D[,j]), function(k) marg_j(j = j, dag = dag_prop, X.lev = X.lev, I = which(I.cal[,k] == 1), X = X[D[,j] == k,], a)))
      
      marg      = sum(sapply(unique(D[,j]), function(k) marg_j(j = i, dag = dag, X.lev = X.lev, I = which(I.cal[,k] == 1), X = X[D[,j] == k,], a))) +
                  sum(sapply(unique(D[,j]), function(k) marg_j(j = j, dag = dag, X.lev = X.lev, I = which(I.cal[,k] == 1), X = X[D[,j] == k,], a)))
      
    }
    
    ## 1.4 Compute MH ratio and accept/reject proposed DAG
    
    ratio_dag = min(0, marg_star - marg + logprior.ratio)
    
    if(log(runif(1)) < ratio_dag){
      
      dag = dag_prop
      
    }
    
    DAG_post[,,s] = dag
    
    
    #######################
    ## 2. Update targets ##
    #######################
    
    ## 2.1 Propose and accept/reject new target
    
    I.cal.prop = I.cal
    
    target_prop = rep(NA, K)


    for(k in int_cont){

      # Propose a target I_k (add or remove a target node from the set I_k)

      target_prop[k] = sample(1:q, 1)

      if(I.cal[target_prop[k],k] == 0){

        I.cal.prop[target_prop[k],k] = 1

      }else{

        I.cal.prop[target_prop[k],k] = 0

      }

      D.k.prop = matrix(0, n_all[k], q)

      D.k.prop[,which(I.cal.prop[,k] == 1)] = k

      out_D_tmp = out_D

      out_D_tmp[[k]] = D.k.prop

      D.prop = do.call(rbind, out_D_tmp)


      j = target_prop[k]
      
      marg_I_k_prop = sum(sapply(unique(D.prop[,j]), function(h) marg_j(j = j, dag = dag, X.lev = X.lev, I = which(I.cal.prop[,h] == 1), X = X[D.prop[,j] == h,], a)))
      marg_I_k      = sum(sapply(unique(D[,j]), function(h) marg_j(j = j, dag = dag, X.lev = X.lev, I = which(I.cal[,h] == 1), X = X[D[,j] == h,], a)))

      logprior_k = lgamma(a_k + sum(I.cal.prop[,k])) + lgamma(q - sum(I.cal.prop[,k]) + b_k) -
                   lgamma(a_k + sum(I.cal[,k])) - lgamma(q - sum(I.cal[,k]) + b_k)


      ratio_k = min(0, marg_I_k_prop - marg_I_k + logprior_k)

      # accept move

      if(log(runif(1)) < ratio_k){

        I.cal[,k] = I.cal.prop[,k]

      }

      out_D[[k]] = matrix(0,n_all[k],q)
      out_D[[k]][,which(I.cal[,k] == 1)] = k

      D = do.call(rbind, out_D)

    }
    
    Targets_post[,,s] = I.cal
      
    
    utils::setTxtProgressBar(pb, s)
    close(pb)
    
  }
  
  return(out = list(X = X,
                    DAG_post = DAG_post[,,(burn + 1):S],
                    Targets_post = Targets_post[,,(burn + 1):S]))
  
}

