library(Rcpp)
sourceCpp("sampling.cpp")

rep_row = function(vec, nrow){
    matrix(vec, nrow = nrow, ncol = length(vec), byrow = T)
}

rep_col = function(vec, ncol){
    matrix(vec, nrow = length(vec), ncol = ncol)
}


nrcd = (NT - Nburn) / nthin
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
upd = rep(0.0, 4 * ng + ns + 4)

sampling(x, s, NT, nthin, nupd, Nburn, g, nu.g, tausq.g, delta, beta.delta, sigma0sq, beta0, sigma1sq, beta1, z, gammam, clu, beta.pi, alpha.lambda, g.seq, delta.seq, sigma0sq.seq, sigma1sq.seq, z.seq, c.seq, gamma.seq, beta.pi.seq, alpha.lambda.seq, nu.g.seq, tausq.g.seq, beta.delta.seq, beta0.seq, beta1.seq, llambda, lpi, lobs, lother, lpost, m.g, dsq.g, a.g, b.g, alpha.delta, alpha0, alpha1, a.l, b.l, nu1, tau1sq, nu2, tau2sq, alpha.pi)
