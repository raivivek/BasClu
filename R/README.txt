##### variables need to define in R before running BasClu.R
# x: a matrix of log(counts), with rows representing genes and columns representing cells 
# s: log library size
# NT: number of MCMC iterations
# nthin: for chain thinning; one sample is kept for every nthin iterations
# nupd: for adaptive adjustment of random walk step in the burn-in stage; random walk steps are updated once every nupd iterations in the burn-in stage 
# Nburn: number of iterations for burn-in
# g: initial values for vector g (length equals number of genes; detailed definition of this variable and following ones are given in the method session of the main manuscript)
# nu.g: initial value for parameter nu.g
# tausq.g: initial value for parameter tausq.g
# delta: initial value for vector delta (length equals number of genes)
# beta.delta: initial value for parameter beta.delta sigma0sq, beta0, sigma1sq, beta1, z, # gammam: initial value for (gamma1, gamma2)
# clu: initial value for cell clustering (length equals number of cells)
# beta.pi: initial value for parameter beta.pi
# alpha.lambda: initial value for parameter alpha.lambda
# m.g: value of hyper-parameter m.g
# dsq.g: value of hyper-parameter d.g^2
# a.g: value of hyper-parameter a.g
# b.g: value of hyper-parameter b.g
# alpha.delta: value of hyper-parameter alpha.delta
# alpha0: value of hyper-parameter alpha0
# alpha1: value of hyper-parameter alpha1
# a.l: value of hyper-parameter a.l
# b.l: value of hyper-parameter b.l
# nu1: value of hyper-parameter nu1
# tau1sq: value of hyper-parameter tau.1^2
# nu2: value of hyper-parameter nu2
# tau2sq: value of hyper-parameter tau.2^2
# alpha.pi: value of hyper-parameter alpha.pi

##### source( "BasClu.R" )
# g.seq, delta.seq, sigma0sq.seq, sigma1sq.seq, z.seq, c.seq, gamma.seq, beta.pi.seq, alpha.lambda.seq, nu.g.seq, tausq.g.seq, beta.delta.seq, beta0.seq, beta1.seq
# are the posterior samples