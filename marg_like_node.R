#########################
## Auxiliary functions ##
#########################

library(gRbase)

## find parents of node j in dag

pa = function(j, dag){
  ifelse(all(dag[,j] == 0),
         return(NULL),
         return(as.numeric(which(dag[,j] != 0))))
}

## find the family of node j in dag

fa = function(j, dag){
  return(as.numeric(c(j, which(dag[,j] != 0))))
}


#########################################
## Compute the DAG marginal likelihood ##
#########################################

## In our MCMC scheme we compare two DAGs (Dag and Dag_star) which differ locally (by the insertion, removal, reversal of a single edge u -> j)
## Because the marginal likelihood admits a node-by-node factorization this will simplify in the ratio between Dag and Dag_star except for component j
## We can implement a formula which computes the marginal likelihood relative to component j only

library(plyr)

marg_j = function(j, a, X, D, I, X.lev, dag){

  ## Notice that there are configurations that are never observed in the data and therefore not contained in table N
  ## We need to include for these configurations the corresponding zero frequencies in the formulas below
  ## Therefore I distinguish between observed and unobserved configurations/frequencies
  
  dag.int = dag
  dag.int[,I] = 0
  
  paj = pa(j = j, dag = dag.int)
  faj = fa(j = j, dag = dag.int)
  
  X.faj  = count(X[,faj])
  
  N.faj.obs   = X.faj$freq
  N.faj.zeros = rep(0, prod(X.lev[faj]) - length(N.faj.obs))
  
  N.faj = c(N.faj.obs, N.faj.zeros)
  
  if(length(paj) == 0){
    
    return((lgamma(a) - lgamma(a + sum(N.faj))) + 
             sum(lgamma(a/length(N.faj) + N.faj)) - length(N.faj)*lgamma(a/length(N.faj)))
    
  }else{
    
    N.paj.obs   = aggregate(X.faj$freq, by = as.list(X.faj[1+1:length(paj)]), FUN = sum)$x
    N.paj.zeros = rep(0, prod(X.lev[paj]) - length(N.paj.obs))
    
    N.paj = c(N.paj.obs, N.paj.zeros)
    
    return((length(N.paj)*lgamma(a/length(N.paj)) - sum(lgamma(a/length(N.paj) + N.paj))) + 
             sum(lgamma(a/length(N.faj) + N.faj)) - length(N.faj)*lgamma(a/length(N.faj)))
    
  }
  
}

score_int = function(dag_mat, j, I, data, a){
  
  if(length(I) == 0){
    
    dag_bn = bn_tmp
    amat(dag_bn) = dag_mat
    
    score(dag_bn, data, targets = names(dag_bn$nodes)[j], type = "bde", by.node = FALSE, iss = a)
    
  }else{
    
    N.j = count(data[,j])$freq
    
      return((lgamma(a) - lgamma(a + sum(N.j))) + 
               sum(lgamma(a/length(N.j) + N.j)) - length(N.j)*lgamma(a/length(N.j)))
    
  }
  
}
