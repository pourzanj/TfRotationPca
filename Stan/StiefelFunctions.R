library(rstan)

n <- 3
p <- 2
d <- n*p - p*(p+1)/2

data <- list(n = n, p = p, d = d)

fit <- stan(file = "StiefelFunctions.stan", data = data, chains = 4, iter = 100000)
