require(tidyverse)
require(expm)

phi <- function(x) dnorm(x)
Phi <- function(x) pnorm(x)

generate_synthetic_data <- function(U, L, c) {
  
  eta <- U %*% diag(L) %*% t(U) + c
  m <- nrow(U)
  p <- ncol(U)
  
  Y <- rbinom(m*m, 1, Phi(eta)) %>% matrix(nrow = m)
  
  Y[Y == 0] <- -1
  diag(Y) <- 0
  Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)] #symmetrize
  
  Y
}

# Y should be m x m w/ 0 on the diagonal. 1 when
# event occured and 0 when event didn't occur. We're taking
# lower diaganol only so it doesn't matter if it's not symmetric
generate_lp_eigennetwork <- function(Y) {
  
  lp <- function(U, L, c) {
      
    eta <- U %*% diag(L) %*% t(U) + c
    m <- nrow(U)
    
    # only get lower diagonal because the matrix is symmetric
    # and we don't want to count connections twice in the likelihood
    Y_eta <- (Y*eta)[lower.tri(Y)]
    sum(log(Phi(Y*eta))[lower.tri(Y)]) - sum(L^2)/(2*m) - c^2/200
  }
}

generate_grad_lp_eigennetwork <- function(Y) {
  
  grad_lp <- function(U, L, c) {
    
    eta <- U %*% diag(L) %*% t(U) + c
    m <- nrow(U)
    p <- ncol(U)
    
    dlp_deta <- Y*phi(Y*eta)/Phi(Y*eta)
    dlp_dU <- matrix(0, nrow = m, ncol = p)
    for(r in 1:p) {
      for(i in 1:m) {
        dlp_dU[i,r] <- sum(dlp_deta[i,]*U[,r]*L[r])
      }
    }
    
    dlp_dL <- rep(1, p)
    for(r in 1:p) {
      U_times_U <- matrix(U[,r], ncol=1) %*% matrix(U[,r], nrow=1)
      dlp_dL[r] <- sum((dlp_deta*U_times_U)[lower.tri(Y)]) - L[r]/m
    }
    
    dlp_dc <- sum(dlp_deta[lower.tri(Y)]) - c/100
      
    list(dlp_dU = dlp_dU, dlp_dL = dlp_dL, dlp_dc = dlp_dc)
  }
}

stiefel_orthogonalize <- function(U, X) U - 0.5*X %*% (t(X) %*% U + t(U) %*% X)

compute_stiefel_flow <- function(X, V, h) {
 
  p <- ncol(X)
  Id <- diag(p)
  Z <- matrix(0, nrow = p, ncol = p)
    
  A <- t(X) %*% V
  S0 <- t(V) %*% V
  
  X_V <- cbind(X, V) %*% expm(h * rbind(cbind(A, -S0), cbind(Id, A))) %*% rbind(cbind(expm(-h*A), Z), cbind(Z, expm(-h*A)))
  X <- X_V[,1:p]
  V <- X_V[,(p+1):(2*p)]
  
  # if p is 1 we need to convert back to matrix
  if(p == 1) {
    X <- matrix(X, ncol = 1)
    V <- matrix(V, ncol = 1)
  }
  
  list(X=X, V=V)
}

emhmc_onesample_eigennetwork <- function(U0, L0, c0, lp, grad_lp, h, nsteps) {
  
  m <- nrow(U0)
  p <- ncol(U0)
  
  # draw momentum for each parameter
  v_U <- rnorm(m*p) %>% matrix(nrow = m)
  v_L <- rnorm(p)
  v_c <- rnorm(1)
  
  # orthogonalize velocity for U
  v_U <- stiefel_orthogonalize(v_U, U0)
  
  # compute total energy
  v <- c(v_U, v_L, v_c)
  H0 <- lp(U0, L0, c0) - 0.5*sum(v^2)

  # copy parameters
  U <- U0
  L <- L0
  c <- c0
  
  # simulate Hamiltonian Dynamics
  for(t in 1:nsteps) {
    
    # update momentum and reorthoganolize
    grad_lp_q <- unlist(grad_lp(U, L, c))
    v_U <- v_U + (h[1]/2)*matrix(grad_lp_q[1:(m*p)], nrow=m)
    v_U <- stiefel_orthogonalize(v_U, U)
    v_L <- v_L + (h[2]/2)*grad_lp_q[(m*p+1):(m*p+p)]
    v_c <- v_c + (h[3]/2)*grad_lp_q[m*p+p+1]
    
    # update position for U
    stiefel_flow <- compute_stiefel_flow(U, v_U, h[1])
    U <- stiefel_flow$X
    v_U <- stiefel_flow$V
    
    # update position for L and c
    L <- L + h[2]*v_L
    c <- c + h[3]*v_c
    
    # update momentum and reorthoganolize
    grad_lp_q <- unlist(grad_lp(U, L, c))
    v_U <- v_U + (h[1]/2)*matrix(grad_lp_q[1:(m*p)], nrow=m)
    v_U <- stiefel_orthogonalize(v_U, U)
    v_L <- v_L + (h[2]/2)*grad_lp_q[(m*p+1):(m*p+p)]
    v_c <- v_c + (h[3]/2)*grad_lp_q[m*p+p+1]
  }
  
  v <- c(v_U, v_L, v_c)
  H1 <- lp(U, L, c) - 0.5*sum(v^2)
  u <- runif(1)
  
  if(u < exp(H1 - H0)) {
    return(list(U = U, L = L, c = c))
  } else {
    return(list(U = U, L = L, c = c))
  }
}

lp_givens_uniform <- function(U) 1.0
grad_lp_givens_uniform <- function(U) matrix(0.0, nrow = nrow(U), ncol = ncol(U))

emhmc_onesample_givens_uniform <- function(U0, lp, grad_lp, h, nsteps) {
  
  m <- nrow(U0)
  p <- ncol(U0)
  
  # draw momentum
  v_U <- rnorm(m*p) %>% matrix(nrow = m)
  
  # orthogonalize velocity for U
  v_U <- stiefel_orthogonalize(v_U, U0)
  
  # compute total energy
  v <- c(v_U)
  H0 <- lp(U0) - 0.5*sum(v^2)
  
  # copy parameters
  U <- U0
  
  # simulate Hamiltonian Dynamics
  for(t in 1:nsteps) {
    
    # update momentum and reorthoganolize
    grad_lp_q <- unlist(grad_lp(U))
    v_U <- v_U + (h[1]/2)*matrix(grad_lp_q[1:(m*p)], nrow=m)
    v_U <- stiefel_orthogonalize(v_U, U)
    
    # update position for U
    stiefel_flow <- compute_stiefel_flow(U, v_U, h[1])
    U <- stiefel_flow$X
    v_U <- stiefel_flow$V
    
    # update momentum and reorthoganolize
    grad_lp_q <- unlist(grad_lp(U))
    v_U <- v_U + (h[1]/2)*matrix(grad_lp_q[1:(m*p)], nrow=m)
    v_U <- stiefel_orthogonalize(v_U, U)
  }
  
  v <- c(v_U)
  H1 <- lp(U) - 0.5*sum(v^2)
  u <- runif(1)
  
  if(u < exp(H1 - H0)) {
    return(list(U = U))
  } else {
    return(list(U = U))
  }
}


emhmc_sample_eigennetwork <- function(nsamples, U0, L0, c0, lp, grad_lp, h, nsteps) {
  
  m <- nrow(U0)
  p <- ncol(U0)
  
  # copy parameters
  U <- U0
  L <- L0
  c <- c0
  
  samples <- matrix(0.0, nrow = nsamples+1, ncol = m*p + p + 1)
  samples[1, 1:(m*p)] <- as.vector(U)
  samples[1, (m*p+1):(m*p+p)] <- L
  samples[1, m*p+p+1] <- c
  
  for(s in 1:nsamples) {
    sample <- emhmc_onesample_eigennetwork(U, L, c, lp, grad_lp, h, nsteps)
    U <- sample$U
    L <- sample$L
    c <- sample$c
    
    samples[s+1, 1:(m*p)] <- as.vector(U)
    samples[s+1, (m*p+1):(m*p+p)] <- L
    samples[s+1, m*p+p+1] <- c
  }
  
  as_tibble(samples) %>%
    set_names(c(paste0("U", rep(1:m, p), "_", rep(1:p, each=m)),
                paste0("L", 1:p),
                "c"))
}

emhmc_sample_givens_uniform <- function(nsamples, U0, lp, grad_lp, h, nsteps) {
  
  m <- nrow(U0)
  p <- ncol(U0)
  
  # copy parameters
  U <- U0
  
  samples <- matrix(0.0, nrow = nsamples+1, ncol = m*p)
  samples[1, 1:(m*p)] <- as.vector(U)
  
  for(s in 1:nsamples) {
    sample <- emhmc_onesample_givens_uniform(U, lp, grad_lp, h, nsteps)
    U <- sample$U
    
    samples[s+1, 1:(m*p)] <- as.vector(U)
  }
  
  as_tibble(samples) %>%
    set_names(paste0("U", rep(1:m, p), "_", rep(1:p, each=m)))
}