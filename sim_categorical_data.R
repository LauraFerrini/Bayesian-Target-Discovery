library(gtools)
library(prodlim)

## find the parent set of node j in dag

pa = function(j, dag){
  ifelse(all(dag[,j] == 0),
         return(NULL),
         return(as.numeric(which(dag[,j] != 0))))
}

## find the family of node j in dag

fa = function(j, dag){
  return(as.numeric(c(j, which(dag[,j] != 0))))
}


# The following function generates the data when all the variables are dichotomous 
# it should be extended for cases in which |Xj| > 2. 

# lev_set = rbind(rep(0,q), rep(1,q)); colnames(lev_set) = 1:q # input 

gen.theta <- function(q, A, a) {
  
  lev_size = rep(2, q) # input for |Xj| > 2
  
  fa_set   = sapply(1:q, function(j) fa(j, dag = A))
  fa_size  = sapply(1:q, function(j) length(fa_set[[j]]))
  Xfa_size = sapply(1:q, function(j) prod(lev_size[fa(j, dag = A)]))
  
  pa_set   = sapply(1:q, function(j) pa(j, dag = A))
  pa_size  = sapply(1:q, function(j) length(pa_set[[j]]))
  Xpa_size = sapply(1:q, function(j) prod(lev_size[pa(j, dag = A)]))
  
  theta = vector(mode = "list", length = q)
  names(theta) = as.character(1:q)
  

  # this creates theta
  
  for(j in 1:q){
    
    a_j = a/Xfa_size[j]
    
    theta_j = cbind(expand.grid(rep(list(0:1), fa_size[[j]])), NA)
    colnames(theta_j) = c(fa_set[[j]], "theta")
    
    theta_j[,fa_size[j]+1] = c(t(rdirichlet(n = Xpa_size[j], rep(a_j, lev_size[j]))))
    
    theta[[j]] = theta_j
    
  }
  
  return(theta)
  
}
  
# this creates the dataset
  
gen.X <- function(n, q, theta, A) {
  
  fa_set   = sapply(1:q, function(j) fa(j, dag = A))
  fa_size  = sapply(1:q, function(j) length(fa_set[[j]]))
  
  pa_set   = sapply(1:q, function(j) pa(j, dag = A))
  pa_size  = sapply(1:q, function(j) length(pa_set[[j]]))
  
  X = matrix(NA, n, q)
    
  for(i in 1:n){
    
    for(j in rev(1:q)){
        
      if(pa_size[j] == 0){
        X[i,j] = rbinom(n = 1, size = 1, prob = theta[[j]]$theta[2]) # this does not work for |Xj| > 2
      }
        
      else{
          
        paj_obs = X[i, pa_set[[j]]]
        paj_conf = theta[[j]][, 2:(pa_size[j] +1)]
          
        if(length(pa_set[[j]]) == 1 ){
          
          idx_match = which(paj_obs == paj_conf)[2] # this does not work for |Xj| > 2
        }
          
        else{
          
          idx_match = row.match(paj_obs, paj_conf) + 1 # this does not work for |Xj| > 2
        }
          
        X[i,j] = rbinom(n = 1, size = 1, prob = theta[[j]]$theta[idx_match])
      }  
    }
  }
  return(X)
}

