LeftRotateCounterClockwise <- function(A, angle, i, j){
  
  #rotating will keep all rows the same except for row i and j
  RA <- A
  RA[i,] <- cos(angle)*A[i,] - sin(angle)*A[j,]
  RA[j,] <- sin(angle)*A[i,] + cos(angle)*A[j,]
  
  return(RA)
}

RightRotateCounterClockwise <- function(A, angle, i, j){
  
  #rotating will keep all rows the same except for row i and j
  AR <- A
  AR[,i] <- cos(angle)*A[,i] + sin(angle)*A[,j]
  AR[,j] <- -sin(angle)*A[,i] + cos(angle)*A[,j]
  
  return(AR)
}

InverseGivensTransform <- function(angles, n, p) {
  G <- diag(n)
  idx <- 1
  pp <- p
  if(n == p) pp <- p-1
  for(i in 1:pp) {
    for(j in (i+1):n) {
      G <- RightRotateCounterClockwise(G, angles[idx], i, j)
      idx <- idx + 1
    }
  }
  
  return(G[,1:p])
}

GivensTransform <- function(W) {
  n <- nrow(W)
  p <- ncol(W)
  A <- W
  angles <- rep(0, n*p-p*(p+1)/2)
  idx <- 1
  if(n == p) p <- p-1
  for(i in 1:p) {
    for(j in (i+1):n) {
      #angle should be a "negative" rotation counter clockwise
      angle <- atan2(A[j,i], A[i,i])
      A <- LeftRotateCounterClockwise(A, -angle, i, j)
      angles[idx] <- angle
      idx <- idx + 1
    }
  }
  
  return(angles)
}

# not finished/working
mirror <- function(W) {
  n <- nrow(W)
  p <- ncol(W)
  theta <- GivensTransform(W)
  theta_new <- theta
  
  idx_lat <- 1
  for(i in 1:p) {
    theta_lat <- theta[idx_lat]
    if(abs(theta_lat) > pi/2) {
      if(theta_lat > 0.0) {
        theta_new[idx_lat] <- -pi/2 + (theta[idx_lat] - 1)
        theta_new[idx_lat:(idx_lat + (n-i-1))] <-
          theta[idx_lat:(idx_lat + (n-i-1))]
      } else{
        
      }
    }
    idx_lat <- idx_lat + (n-i)
  }
}