library(rstan)

n <- 3
p <- 2
d <- n*p - p*(p+1)/2
data <- list(n = n, p = p, d = d)

fit <- stan(file = "Stan/UniformStiefel.stan", data = data, chains = 4, iter = 100000)

#fit <- stan(file = "Stan/StiefelFunctions.stan", data = data, chains = 4, iter = 100000)


