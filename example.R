###############################################
## Simulated categorical interventional data ##
###############################################

library(BCDAG)

## Example

set.seed(1234)
q = 10

A = rDAG(q = q, w = 0.2)
rownames(A) = colnames(A) = 1:q
a = 1

source("sim_int_categorical_data.R")

set.seed(123)
theta = gen.theta(q, A, a)

# theta is the parameter indexing the observational distributions

# Fix sample sizes and targets

n_all = c(500, 500)
I.cal = cbind(c(0,0,0,0,0,0,0,0,0,0),
              c(0,0,0,1,0,0,0,1,0,0))

# Generate interventional datasets

data_int = gen.dataset(theta, A, a, n_all, I.cal)

X = data_int$X
D = data_int$D

# D is an (n, q) matrix associating observations and nodes to intervention targets:
# D[i,j] = k means that i-th observations has been produced under context k with intervention on node j

# Finds the level (unique values) of each cateogorical variable:

X.lev = sapply(1:ncol(X), function(j) length(unique(X[,j])))


###########################################################################
## Run the MCMC for posterior inference on DAGs and intervention targets ##
###########################################################################

# Set prior hyperparameters

a_dag = b_dag = 1
a_k = b_k = 1

ne = 5
a  = 1

## Fix number of MCMC iterations and burn-in

S    = 6000
burn = 1000

source("mcmc_dags_targets.R")

out = mcmc_dag_targets(X, S, burn, ne = ne, a = NULL, a_k = a_k, b_k = b_k, a_dag, b_dag, n_all = n_all, X.lev = X.lev, int_cont = 2)

dim(out$DAG_post)
dim(out$Targets_post)


#########################
## Posterior summaries ##
#########################

## Posterior probabilities of intervention

I.cal.hat = apply(X = out$Targets_post, FUN = mean, MARGIN = c(1,2))
round(I.cal.hat, 2)

I.cal

## Posterior probabilities of edge inclusion

P_hat = apply(X = out$DAG_post, FUN = mean, MARGIN = c(1,2))
round(P_hat, 2)

## Median Probability DAG Model (MPM)

A_hat = round(P_hat)
A_hat

## Compare true DAG and posterior probabilities of edge inclusion

par(mfrow = c(1,2))

image(A)
image(P_hat)

## Compare true target and posterior probabilities of intervention

par(mfrow = c(1,2))

image(I.cal)
image(I.cal.hat)
