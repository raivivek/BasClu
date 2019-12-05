rep_row = function(vec, nrow) {
  matrix(vec,
         nrow = nrow,
         ncol = length(vec),
         byrow = T)
}

rep_col = function(vec, ncol) {
  matrix(vec, nrow = length(vec), ncol = ncol)
}


init_params <<- function() {
  #############################################
  # Hyperparameters
  #############################################
  # alpha.pi: value of hyper-parameter alpha.pi
  alpha.pi <<- 1.0
  a.d <<- 2.0
  b.d <<- 1.0

  # a.l: value of hyper-parameter a.l
  # b.l: value of hyper-parameter b.l
  a.l <<- 2.0
  b.l <<- 1.0

  # a.p: value of hyper-parameter a.p
  # b.p: value of hyper-parameter b.p
  a.p <<- 2.0
  b.p <<- 1.0

  # m.g: value of hyper-parameter m.g
  # dsq.g: value of hyper-parameter d.g^2
  m.g <<- 0.0
  dsq.g <<- 100.0

  # nu.g: initial value for parameter nu.g
  # USES: m.g, dsq.g
  nu.g <<- stats::rnorm(m.g, dsq.g)

  # a.g: value of hyper-parameter a.g
  # b.g: value of hyper-parameter b.g
  a.g <<- 1.0
  b.g <<- 0.01

  #alpha.delta: value of hyper-parameter alpha.delta
  #alpha0: value of hyper-parameter alpha0
  #alpha1: value of hyper-parameter alpha1
  alpha.delta <<- 2.0
  alpha0 <<- 1.0
  alpha1 <<- 1.0

  # nu1: value of hyper-parameter nu1
  # tau1sq: value of hyper-parameter tau.1^2
  nu1 <<- 0.0
  tau1sq <<- 100.0

  # nu2: value of hyper-parameter nu2
  # tau2sq: value of hyper-parameter tau.2^2
  nu2 <<- 0.0
  tau2sq <<- 100.0


  ################################################
  # Priors
  ################################################

  # beta.delta: initial value for parameter beta.delta, beta0, beta1
  beta.delta <<- stats::rgamma(1, shape = a.d, rate = b.d)
  beta0 <<- stats::rgamma(1, shape = a.d, rate = b.d)
  beta1 <<- stats::rgamma(1, shape = a.d, rate = b.d)

  # nu.g: initial value for parameter nu.g
  nu.g <<- stats::rnorm(1, m.g, dsq.g)

  # beta.pi: initial value for parameter beta.pi
  beta.pi <<- stats::rgamma(1, shape = a.p, rate = b.p)

  # tausq.g: initial value for parameter tausq.g
  tausq.g <<- invgamma::rinvgamma(n = 1, shape = a.g, rate = b.g)

  # alpha.lambda: initial value for parameter alpha.lambda
  alpha.lambda <<- stats::rgamma(1, shape = a.l, rate = b.l)

  # gammam: initial value for (gamma1, gamma2)
  gammam <<- mapply(stats::rnorm, c(1, 1), c(nu1, nu2), c(tau1sq, tau2sq))
}


initialize_values <- function() {
  #############################################
  # Input dependent parameters
  #############################################
  g <<- stats::rnorm(ng, nu.g, tausq.g) # g: initial value of g (length = no. of genes)
  sigma0sq <<- invgamma::rinvgamma(ng, shape = alpha0, rate = beta0) # sigma0sq: initial value
  sigma1sq <<- invgamma::rinvgamma(ng, shape = alpha1, rate = beta1) # sigma1sq: initial value
  delta <<- stats::rgamma(ng, alpha.delta, beta.delta) # delta: initial value of delta (length = no. of genes)

  clu <<-  stats::rbinom(n = ns, size = 4, prob = 1/4) # clu: initial value for cell clustering (length = number of cells)
  z <<- matrix(do.call(cbind, lapply(1:ns, function(x) stats::rbeta(ng, alpha.pi, beta.pi))))
}

## Posterior samples
initialize_posteriors <- function() {
  g.seq <<- matrix(0.0, ng, nrcd)
  delta.seq <<- matrix(0.0, ng, nrcd)
  sigma0sq.seq <<- matrix(0.0, ng, nrcd)
  sigma1sq.seq <<- matrix(0.0, ng, nrcd)
  c.seq <<- matrix(as.integer(0), ns, nrcd)
  z.seq <<- matrix(0.0, ng * ns, nrcd)
  gamma.seq <<- matrix(0.0, 2, nrcd)
  mu.seq <<- rep(0.0, nrcd)
  alpha.lambda.seq <<- rep(0.0, nrcd)
  beta.pi.seq <<- rep(0.0, nrcd)
  nu.g.seq <<- rep(0.0, nrcd)
  tausq.g.seq <<- rep(0.0, nrcd)
  beta.delta.seq <<- rep(0.0, nrcd)
  beta0.seq <<- rep(0.0, nrcd)
  beta1.seq <<- rep(0.0, nrcd)
  llambda <<- rep(0.0, nrcd)
  lpi <<- rep(0.0, nrcd)
  lobs <<- rep(0.0, nrcd)
  lother <<- rep(0.0, nrcd)
  lpost <<- rep(0.0, nrcd)
  nupd <<- rep(0.0, 4 * ng + ns + 4)
}


#' Run BasClu.
#'
#' @param x Log-transformed count matrix
#' @param s Log library sizes
#' @param NT Number of MCMC iterations
#' @param nthin For thinning: keep sample every nthin iteration
#' @param nupd For updating random walk during burn-in
#' @param Nburn Number of iterations for burn-in
#' @param clustering Clustering method to use
#' @param ... Additional params
#'
#' @return z Z-matrix
#' @export
basclu <- function(x, s, NT, nthin, nupd, Nburn, clustering = "BasCluZ", ...) {
  stopifnot(is.matrix(x) && is.matrix(s))

  # MCMC Parameters
  NT <- 10000   # NT: number of MCMC iterations
  nthin <- 5    # nthin: for chain thinning; one sample is kept for every nthin iterations
  nupd <- 10   # nupd: random walk steps are updated once every nupd iterations in the burn-in stage
  Nburn <- 1000  # Nburn: number of iterations for burn-in
  nrcd = (NT - Nburn) / nthin

  # Get dimensions
  ng <<- dim(x)[1]
  ns <<- dim(x)[2]

  # Setup all parameters
  init_params()
  initialize_values()
  initialize_posteriors()

  # Call CPP code
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

  return (z.seq)
}
