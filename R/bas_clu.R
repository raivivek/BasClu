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

NT <- 10000   # NT: number of MCMC iterations
nthin <- 5    # nthin: for chain thinning; one sample is kept for every nthin iterations
nupd <- 10   # nupd: random walk steps are updated once every nupd iterations in the burn-in stage
Nburn <- 100  # Nburn: number of iterations for burn-in
nrcd = (NT - Nburn) / nthin


## DEBUG
df_mat <- as.matrix(df[, -(1:3)])
x <- log(df_mat + 1)
x <- x[1:100, ]

s <- log(matrix(colSums(df_mat), nrow = 1))
#s <- log(colSums(df_mat))

ng <- dim(x)[1]
ns <- dim(x)[2]

#############################################
# Input dependent parameters
#############################################

g <- rnorm(ng, nu.g, tausq.g) # g: initial value of g (length = no. of genes)
sigma0sq <- invgamma::rinvgamma(ng, shape = alpha0, rate = beta0) # sigma0sq: initial value
sigma1sq <- invgamma::rinvgamma(ng, shape = alpha1, rate = beta1) # sigma1sq: initial value
delta <- rgamma(ng, alpha.delta, beta.delta) # delta: initial value of delta (length = no. of genes)

clu <-  rbinom(n = ns, size = 4, prob = 1/4) # clu: initial value for cell clustering (length = number of cells)
z <- matrix(0, ng, ns)



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
print("Starting sampling now...")
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
