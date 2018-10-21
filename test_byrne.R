source("byrne.R")

#####
# Generate
print("~~~Generating data and lp")

m <- 270
p <- 3

U <- diag(m)[,1:p]
L <- c(120, 100, 80)
c <- 60.0

Y <- generate_synthetic_data(U, L, c)
print("Y: ")
print(Y[1:11,1:11])

lp <- generate_lp(Y)
print(lp(U,L,c))

#####
# Gradients
print("~~~Testing Gradients")

grad_lp <- generate_grad_lp(Y)
g <- grad_lp(U,L,c)

eps <- 1e-4

print("~~Testing dlp/dU_{11,2}...")

U_plus_eps <- U
U_plus_eps[11,2] <- U[11,2] + eps

print(paste("Finite Difference: ", (lp(U_plus_eps, L, c) - lp(U, L, c))/eps))
print(paste("Analytic: ",g$dlp_dU[11,2]))

print("~~Testing dlp/lambda_2...")

L_plus_eps <- L
L_plus_eps[2] <- L_plus_eps[2] + eps

print(paste("Finite Difference: ", (lp(U, L_plus_eps, c) - lp(U, L, c))/eps))
print(paste("Analytic: ",g$dlp_dL[2]))

print("~~Testing dlp/dc...")

eps <- 1e-7
c_plus_eps <- c + eps

print(paste("Finite Difference: ", (lp(U, L, c_plus_eps) - lp(U, L, c))/eps))
print(paste("Analytic: ",g$dlp_dc))

#####
# HMC
print("~~~Testing HMC")

print("~~Testing orthogonalization of stiefel manifold")

# draw random velocity and orthogonalize
V <- rnorm(m*p) %>% matrix(nrow = m)
V_orth <- stiefel_orthogonalize(V, U)

# if part of the tangent space then V^T X + X^T V should equal 0
print("inner-product before orthogonalization: ")
print(t(V) %*% U + t(U) %*% V)
print("inner-product after orthogonalization: ")
print(t(V_orth) %*% U + t(U) %*% V_orth)

print("~~Testing stiefel flow")
stiefel_flow <- compute_stiefel_flow(U, V_orth, 0.005)
X <- stiefel_flow$X
V <- stiefel_flow$V

print("X^T X: ")
print(t(X) %*% X)

print("V^T X + X^T V: ")
print(t(V) %*% X + t(X) %*% V)

print("~~Testing sampler")
sample <- emhmc_onesample_eigennetwork(U, L, c, lp, grad_lp, h = c(0.005, 0.1, 0.001), nsteps = 1)
U <- sample$U
L <- sample$L
c <- sample$c

print("U^T U: ")
print(t(U) %*% U)

print(paste("lambda: ", L))
print(paste("c: ", c))