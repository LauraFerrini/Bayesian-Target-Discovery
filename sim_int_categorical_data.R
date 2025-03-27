library(BCDAG)
source("sim_categorical_data.R")

gen.int.data = function(theta = theta, A = A, a = a, n_k = n_k, I_k, q){
    
    A.int = A
    A.int[, I_k] = 0   # I remove all edges "pointing" to the intervened node (all j --> I)
    
    colnames(A.int) = 1:q
    
    theta_k = theta
    
    if(length(I_k) != 0){
      
      set.seed(123)
      theta_tmp = gen.theta(q = q, A = A.int, a = a)
      theta_k[I_k] = theta_tmp[I_k]
      
    }
    
    X.k = gen.X(n = n_k, q = q, theta = theta_k, A = A.int)
    
    return(X.k = X.k)
    
}


gen.dataset = function(theta, A, a, n_all, I.cal){

  # INPUT:
  
  # theta : true parameters obtained using the function gen.theta()
  # A     : true DAG
  # a     : parameter pf the Dirichlet priors
  # n_all : (K, 1) vector collecting the number of observations per scenario
  # I.cal : (q, K) matrix of indicators for interventional nodes
  
  K = length(n_all)
  q = ncol(A)
  out_X = lapply(1:K, function(k) gen.int.data(theta = theta, A = A, a = a, 
                                               n_k = n_all[k], I_k = which(I.cal[,k] == 1),
                                               q = q))
  out_D = lapply(1:K, function(k) matrix(0, n_all[k], ncol(A)))
  
  for(k in 1:K){
    
    out_D[[k]][,which(I.cal[,k] == 1)] = k
    
  }
  
  X = do.call(rbind, out_X)
  X = data.frame(X)
  
  X_labs = X
  
  X[,1:ncol(X)] = lapply(X[,1:ncol(X)], as.factor)
  
  
  colnames(X) = 1:ncol(A)
  
  return(list(X = X, D = do.call(rbind, out_D)))
  
  
}


## Example

# set.seed(1234)
# q = 10
# 
# A = rDAG(q = q, w = 0.2)
# rownames(A) = colnames(A) = 1:q
# a = 1
# 
# set.seed(123)
# theta = gen.theta(q, A, a)
# X = gen.X(n = 200, q, theta, A)
# X 
# 
# # I_k is the collection of intervention nodes in scenario k
# theta
# I_k = 1
# 
# theta
# n_all = c(200, 100)
# I.cal = cbind(c(0,0,1,0,0,0,0,0,0,0),
#               c(0,0,0,1,0,0,0,0,0,0))
# 
# data_all = gen.dataset(theta, A, a, n_all, I.cal)