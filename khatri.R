# expose_stan_functions(stanc("Stan/givens_optimized.stan"))

# equivalent to \Gamma_n(p/2)
gamma_n <- function(n,p) pi^(n/4*(n-1)) * prod(gamma(1/2*(p - 1:n + 1)))

dX <- function(n, p, theta) gamma_n(n,p) * 2^(-n) * pi^(-p*n/2) * cos(theta[2])

dX <- function(p, n, theta) {
  cos_theta_power <- theta
  idx <- 1
  for(i in 1:n) {
    for(j in (i+1):p) {
      cos_theta_power[idx] <- cos(theta[idx])^(j-i-1)
      idx <- idx + 1
    }
  }
  sum(log(cos_theta_power))
}

# theta <- c(0, pi/4)
# dX(3, 1, theta)
# givens_lp(3, 1, theta)

theta <- c(pi/5, pi/5, pi/5)
dX(3, 2, theta)
givens_lp(3, 2, theta)

####
fit <- stan("Stan/new_khatri/given_uniform.stan", data = list(n=3,p=1), chains=1, iter=1e5)
fit <- stan("givens_uniform_patch.stan", data = list(n=3,p=1), chains=1, iter=1)

#######
# test patch when there's significant mass at the poles
#######
n <- 3
p <- 2
d <- n*p-p*(p+1)/2

source("../../Givens.R")
source("../../GenerateSynthetic.R")
W <- InverseGivensTransform(c(0,0,0), n, p)

N <- 15
sim_3_1 <- GenerateHighDimData(n = n, p = p, W = W, LambdaVec = c(100, 50, 10), sd = 1.0, N = N)

sigmaSqHyperPrior <- 10
xMat <- as.matrix(sim_3_1$x)
SigmaHat <- (1/N)*t(xMat) %*% xMat

dat <- list(n = n, p = p, N = N, SigmaHat = SigmaHat, sigmaSqHyperPrior)
fit <- stan("ppca.stan", data = dat, chains = 1, iter = 1e4)
fit_patch <- stan("ppca_patch.stan", data = dat, chains = 1, iter = 1e4)
