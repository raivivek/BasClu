#############################################
# Hyperparameters
#############################################

# alpha.pi: value of hyper-parameter alpha.pi
alpha.pi <- 1.0
a.d <- 2.0
b.d <- 1.0

# a.l: value of hyper-parameter a.l
# b.l: value of hyper-parameter b.l
a.l <- 2.0
b.l <- 1.0

# a.p: value of hyper-parameter a.p
# b.p: value of hyper-parameter b.p
a.p <- 2.0
b.p <- 1.0

# m.g: value of hyper-parameter m.g
# dsq.g: value of hyper-parameter d.g^2
m.g <- 0.0
dsq.g <- 100.0

# nu.g: initial value for parameter nu.g
# USES: m.g, dsq.g
nu.g <- rnorm(m.g, dsq.g)

# a.g: value of hyper-parameter a.g
# b.g: value of hyper-parameter b.g
a.g <- 1.0
b.g <- 0.01

#alpha.delta: value of hyper-parameter alpha.delta
#alpha0: value of hyper-parameter alpha0
#alpha1: value of hyper-parameter alpha1
alpha.delta <- 2.0
alpha0 <- 1.0
alpha1 <- 1.0

# nu1: value of hyper-parameter nu1
# tau1sq: value of hyper-parameter tau.1^2
nu1 <- 0.0
tau1sq <- 100.0

# nu2: value of hyper-parameter nu2
# tau2sq: value of hyper-parameter tau.2^2
nu2 <- 0.0
tau2sq <- 100.0


################################################
# Priors
################################################

# beta.delta: initial value for parameter beta.delta, beta0, beta1
beta.delta <- rgamma(1, shape = a.d, rate = b.d)
beta0 <- rgamma(1, shape = a.d, rate = b.d)
beta1 <- rgamma(1, shape = a.d, rate = b.d)

# nu.g: initial value for parameter nu.g
nu.g <- rnorm(1, m.g, dsq.g)

# beta.pi: initial value for parameter beta.pi
beta.pi <- rgamma(1, shape = a.p, rate = b.p)

# tausq.g: initial value for parameter tausq.g
tausq.g <- invgamma::rinvgamma(n = 1, shape = a.g, rate = b.g)

# alpha.lambda: initial value for parameter alpha.lambda
alpha.lambda <- rgamma(1, shape = a.l, rate = b.l)

# gammam: initial value for (gamma1, gamma2)
gammam <- mapply(rnorm, c(1, 1), c(nu1, nu2), c(tau1sq, tau2sq))
