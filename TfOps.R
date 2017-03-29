library(tensorflow)
library(purrr)

#will expect angle to be a tensor, n, i_in, and j_in will be
#R variables with zero indexing
CreateRotationMatrix <- function(angle, n, i_in, j_in){
  OneIndices <- lapply((1:n)[!(1:n %in% c(i_in+1,j_in+1))], function(x) list(as.integer(x-1), as.integer(x-1)))
  i <- as.integer(i_in)
  j <- as.integer(j_in)
  TrigIndices <- list(list(i,i), list(i,j), list(j,i), list(j,j))
  
  OneElements <- rep(list(1), n-2)
  TrigElements <- list(cos(angle), -sin(angle), sin(angle), cos(angle))
  
  tf$scatter_nd(c(OneIndices, TrigIndices),
                tf$stack(c(OneElements, TrigElements)),
                shape(as.integer(n), as.integer(n)))
}

#angles should be a shape (np-p(p+1)/2,) tensor i.e.
# Theta01 <- tf$placeholder(tf$float32, shape = shape())
# Theta02 <- tf$placeholder(tf$float32, shape = shape())
# Theta12 <- tf$placeholder(tf$float32, shape = shape())
# Theta <- tf$stack(list(Theta01, Theta02, Theta12))
CreateGivensMatrix <- function(angles, n, p) {
  
  idx <- 0
  #initialize G as the identity then repeatedly
  #mult. on the right by rotations
  G <- tf$constant(diag(n),dtype = tf$float32)
  for(i in 0:(p-1)) {
    for(j in (i+1):(n-1)) {
      R <- CreateRotationMatrix(angles[idx], n, i, j)
      G <- tf$matmul(G, R)
      idx <- idx + 1
    }
  }
  return(G)
}

#take in a vector and a list of rank 0 tensors with which we will take the
#derivative w.r.t.
GetJacobian <- function(v, x){
  n <- v$shape$as_list()[[1]]
  
  #get elements as an R list of TF scalars then apply gradients to each
  Elements <- lapply(0:(n-1), function(i) v[i])
  JacobianList <- lapply(Elements, function(e) tf$gradients(e,x))
  
  Jacobian <- tf$concat(JacobianList, axis = 0L)
  
  return(Jacobian)
}

GetStiefelAreaForm <- function(G, angles, n, p) {
  
  FList <- list()
  
  for(i in 0:(p-1)) {
    #Rows <- tf$transpose(G[,(i+1):(n-1)])
    GTransposeRows <- tf$slice(tf$transpose(G), list(as.integer(i+1), 0L), list(as.integer(n-i-1), as.integer(n)))
    Jacobian <- GetJacobian(G[,i], angles)
    OneForms <- tf$transpose(tf$matmul(GTransposeRows, Jacobian))
    
    for(j in 0:(n-i-1-1)) {
      OneForm <- tf$slice(OneForms, list(0L, as.integer(j)), list(as.integer(n), 1L))
      FList <- c(FList, OneForm)
    }
  }
  
  FMatrix <- tf$concat(FList, axis = 1L)
  det <- tf$matrix_determinant(FMatrix)
  
  return(det)
}


