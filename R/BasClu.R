library(Rcpp)

source("R/init_params.R")
sourceCpp("src/sampling.cpp")

rep_row = function(vec, nrow) {
  matrix(vec,
         nrow = nrow,
         ncol = length(vec),
         byrow = T)
}

rep_col = function(vec, ncol) {
  matrix(vec, nrow = length(vec), ncol = ncol)
}

#############################################
# MCMC Parameters
#############################################
# NT: number of MCMC iterations
NT <- 10000

# nthin: for chain thinning; one sample is kept for every nthin iterations
nthin <- 5

# nupd: adaptive adjustment of random walk step in the burn-in stage;
#  random walk steps are updated once every nupd iterations in the burn-in stage
nupd <- 100

# Nburn: number of iterations for burn-in
Nburn <- 100


## DEBUG
df_mat <- df[, -(1:3)]
x <- log(df_mat + 1)
s <- log(colSums(df_mat))

ng <- dim(df_mat)[0]
ns <- dim(df_mat)[1]

nrcd = (NT - Nburn) / nthin

# g: initial values for vector g (length equals number of genes; detailed definition of
#   this variable and following ones are given in the method session of the main
#   manuscript)

g <- rnorm(ng, nu.g, tausq.g)


# clu: initial value for cell clustering (length equals number of cells)
clu <-


## Posterior samples
g.seq = matrix(0.0, ng, nrcd)
delta.seq = matrix(0.0, ng, nrcd)
sigma0sq.seq = matrix(0.0, ng, nrcd)
sigma1sq.seq = matrix(0.0, ng, nrcd)
c.seq = matrix(as.integer(0), ns, nrcd)
z.seq = matrix(0.0, ng * ns, nrcd)
gamma.seq = matrix(0.0, 2, nrcd)
mu.seq = rep(0.0, nrcd)
alpha.lambda.seq = rep(0.0, nrcd)
beta.pi.seq = rep(0.0, nrcd)
nu.g.seq = rep(0.0, nrcd)
tausq.g.seq = rep(0.0, nrcd)
beta.delta.seq = rep(0.0, nrcd)
beta0.seq = rep(0.0, nrcd)
beta1.seq = rep(0.0, nrcd)
llambda = rep(0.0, nrcd)
lpi = rep(0.0, nrcd)
lobs = rep(0.0, nrcd)
lother = rep(0.0, nrcd)
lpost = rep(0.0, nrcd)
nupd = rep(0.0, 4 * ng + ns + 4)


# hand over to cpp
sampling(
  x,
  s,
  NT,
  nthin,
  nupd,
  Nburn,
  g,
  nu.g,
  tausq.g,
  delta,
  beta.delta,
  sigma0sq,
  beta0,
  sigma1sq,
  beta1,
  z,
  gammam,
  clu,
  beta.pi,
  alpha.lambda,
  g.seq,
  delta.seq,
  sigma0sq.seq,
  sigma1sq.seq,
  z.seq,
  c.seq,
  gamma.seq,
  beta.pi.seq,
  alpha.lambda.seq,
  nu.g.seq,
  tausq.g.seq,
  beta.delta.seq,
  beta0.seq,
  beta1.seq,
  llambda,
  lpi,
  lobs,
  lother,
  lpost,
  m.g,
  dsq.g,
  a.g,
  b.g,
  alpha.delta,
  alpha0,
  alpha1,
  a.l,
  b.l,
  nu1,
  tau1sq,
  nu2,
  tau2sq,
  alpha.pi
)
