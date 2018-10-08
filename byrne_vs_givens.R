samples <- emhmc_sample_eigennetwork(nsamples = 5, U, L, c, lp, grad_lp, h = c(0.005, 0.1, 0.001), nsteps = 1)

Y[Y == -1] <- 0
#fit <- stan(file = "Stan/eigennetwork.stan", data = list(N=270, P=3, Y=Y), chains = 1, iter = 1000, refresh = 1)

N <- 270
P <- 3
