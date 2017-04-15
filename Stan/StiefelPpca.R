library(rstan)

n <- 3
p <- 2
d <- n*p - p*(p+1)/2
sigmaSqHyperPrior <- 10
N <- 15
load(file = "experiments/uncertaintyQuantification/unstableEstimate.Rdata")
xMat <- as.matrix(x[,1:3])
SigmaHat <- (1/N)*t(xMat) %*% xMat
data <- list(n = n, p = p, d = d, sigmaSqHyperPrior = sigmaSqHyperPrior, N = N, SigmaHat = SigmaHat)

fit <- stan(file = "Stan/StiefelPpca.stan", data = data, chains = 4, iter = 10000)

