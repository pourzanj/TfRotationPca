system.time(samples <- emhmc_sample_eigennetwork(nsamples = 1000, U, L, c, lp, grad_lp, h = c(0.005, 0.1, 0.001), nsteps = 1))
qplot(1:1001, samples$U1_1)

Y[Y == -1] <- 0
Y <- Y_Pro
Y[is.na(Y)] <- 0

dat <- list(n=230, p=3, Y=Y)

# set init
n <- 270
p <- 3
d <- n*p - p*(p+1)/2

x <- rep(0.0, D)
y <- rep(0.0, D)
L_rev <- c(3, 2, 1)
c <- 0.1
stan_rdump(c('x', 'y', 'L_rev', 'c'), "eigennetwork_init.Rdump")

# set metric
inv_metric <- c(rep(0.02, D), rep(0.02, D), rep(1, 3), 0.001)
stan_rdump(c('inv_metric'), "eigennetwork_mass_matrix.Rdump")

m <- stan_model(file = "Stan/eigennetwork.stan")

fit_optim <- optimizing(m, data = dat)
fit_vb <- vb(m, data = dat)

fit <- stan(file = "Stan/eigennetwork.stan", data = dat, chains = 1, iter = 1000, refresh = 1)
fit <- stan(file = "../../Stan/new_khatri/eigennetwork.stan", data = dat, chains = 1, iter = 1000, refresh = 1)


