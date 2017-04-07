library(tensorflow)
library(purrr)

#will expect angle to be a tensor, n, i_in, and j_in will be
#R variables with zero indexing
CreateRotationMatrix <- function(angle, n, i_in, j_in){
  with(tf$name_scope(paste0('RotationMatrix', i_in, j_in)), {
    OneIndices <- lapply((1:n)[!(1:n %in% c(i_in+1,j_in+1))], function(x) list(as.integer(x-1), as.integer(x-1)))
    i <- as.integer(i_in)
    j <- as.integer(j_in)
    TrigIndices <- list(list(i,i), list(i,j), list(j,i), list(j,j))
    
    OneElements <- rep(list(1), n-2)
    TrigElements <- list(cos(angle), -sin(angle), sin(angle), cos(angle))
    
    tf$scatter_nd(c(OneIndices, TrigIndices),
                  tf$stack(c(OneElements, TrigElements)),
                  shape(as.integer(n), as.integer(n)))
  })
}

CreateDerivativeOfRotationMatrix <- function(angle, n, i_in, j_in){
  with(tf$name_scope(paste0('DerivativeOfRotationMatrix', i_in, j_in)), {
    i <- as.integer(i_in)
    j <- as.integer(j_in)
    TrigIndices <- list(list(i,i), list(i,j), list(j,i), list(j,j))
    
    TrigElements <- list(-sin(angle), -cos(angle), cos(angle), -sin(angle))
    
    tf$scatter_nd(TrigIndices,TrigElements,shape(as.integer(n), as.integer(n)))
  })
}


#applying this op to a matrix A is equivalent to right multiplying
#that matrix by a rotation matrix that's a rotation of angle
#in the i,j plane.
#the resulting matrix will only differ from A in the ith and jth columns.
#we can thus split A in to the following five parts i, 1L, j-i-1, 1L, n-j-1
# RotateOp <- function(A, angle, n_in, i_in, j_in){
#   n <- as.integer(n_in)
#   i <- as.integer(i_in)
#   j <- as.integer(j_in)
# 
#   with(tf$name_scope(paste0('RotateOp', i_in, j_in)), {
#     SplitCols <- tf$split(A, list(i, 1L, j-i-1L, 1L, n-j-1L), axis = 1L)
#     AColI <- SplitCols[[2]]
#     AColJ <- SplitCols[[4]]
#     
#     SplitCols[[2]] <- cos(angle)*AColI + sin(angle)*AColJ
#     SplitCols[[4]] <- -sin(angle)*AColI + cos(angle)*AColJ
#     
#     Result <- tf$concat(SplitCols, 1L)
#   })
#   
#   return(Result)
# }

#angles should be a shape (np-p(p+1)/2,) tensor i.e.
# Theta01 <- tf$placeholder(tf$float32, shape = shape())
# Theta02 <- tf$placeholder(tf$float32, shape = shape())
# Theta12 <- tf$placeholder(tf$float32, shape = shape())
# Theta <- tf$stack(list(Theta01, Theta02, Theta12))
# CreateGivensMatrix <- function(angles, n, p) {
# 
#   idx <- 0
#   #initialize G as the identity then repeatedly
#   #mult. on the right by rotations
#   G <- tf$constant(diag(n),dtype = tf$float32)
#   for(i in 0:(p-1)) {
#     for(j in (i+1):(n-1)) {
# 
#       #old way building rotation matrices then multiplying
#       R <- CreateRotationMatrix(angles[idx], n, i, j)
#       G <- tf$matmul(G, R, b_is_sparse=TRUE)
# 
#       #new way to apply a rotation cred. Ben
#       #G <- RotateOp(G, angles[idx], n, i, j)
# 
#       idx <- idx + 1
#     }
#   }
#   return(G)
# }

CreateGivensMatrix <- function(PartialRotationsAB, n, p) {
  d <- n*p - p*(p+1)/2
  
  with(tf$name_scope(paste0('GivensMatrix')), {
    GivensMatrix <- tf$identity(PartialRotationsAB$A[[d+1]], name = 'GivensMatrix')
  })
    
  return(GivensMatrix)
}

#take in a vector and a list of rank 0 tensors with which we will take the
#derivative w.r.t.
# GetJacobian <- function(v, x){
#   n <- v$shape$as_list()[[1]]
# 
#   #get elements as an R list of TF scalars then apply gradients to each
#   Elements <- lapply(0:(n-1), function(i) v[i])
#   JacobianList <- lapply(Elements, function(e) tf$gradients(e,x))
# 
#   Jacobian <- tf$concat(JacobianList, axis = 0L)
# 
#   return(Jacobian)
# }

# GetStiefelAreaForm <- function(G, angles, n, p) {
#   
#   FList <- list()
#   StiefelDim <- as.integer(n*p - p*(p+1)/2)
#   
#   for(i in 0:(p-1)) {
#     #Rows <- tf$transpose(G[,(i+1):(n-1)])
#     GTransposeRows <- tf$slice(tf$transpose(G), list(as.integer(i+1), 0L), list(as.integer(n-i-1), as.integer(n)))
#     Jacobian <- GetJacobian(G[,i], angles)
#     OneForms <- tf$transpose(tf$matmul(GTransposeRows, Jacobian))
#     
#     for(j in 0:(n-i-1-1)) {
#       OneForm <- tf$slice(OneForms, list(0L, as.integer(j)), list(as.integer(StiefelDim), 1L))
#       FList <- c(FList, OneForm)
#     }
#   }
#   
#   FMatrix <- tf$concat(FList, axis = 1L)
#   det <- tf$matrix_determinant(FMatrix)
#   
#   return(det)
# }

#need a way to A and B which should be lists of size d=np-p(p+1)/2
#go through these and mult ADB giving an nxp matrix. Stack these
#to get an nxpxd tensor the back slices are the nxd slices of Jacobian
#we need
CreatePartialGivens <- function(angles, n, p) {
  
  with(tf$name_scope(paste0('PartialGivensA')), {
    G <- tf$constant(diag(n),dtype = tf$float32, name = 'Id')
    A <- list(G)
    idx <- 0
    for(i in 0:(p-1)) {
      for(j in (i+1):(n-1)) {
        R <- CreateRotationMatrix(angles[idx], n, i, j)
        G <- tf$matmul(G, R, b_is_sparse=TRUE, name = paste0('P_',0,1,'_',i,j))
        A <- c(A, G)
        idx <- idx + 1
      }
    }
  })
  
  StiefelDim <- as.integer(n*p - p*(p+1)/2)
  with(tf$name_scope(paste0('PartialGivensB')), {
    G <- tf$constant(diag(n),dtype = tf$float32, name = 'Id')[,0:(p-1)]
    B <- list(G)
    idx <- StiefelDim - 1
    for(i in (p-1):0) {
      for(j in (n-1):(i+1)) {
        R <- CreateRotationMatrix(angles[idx], n, i, j)
        G <- tf$matmul(R, G, a_is_sparse=TRUE, name = paste0('P_',i,j,'_',p-1,n-1))
        B <- c(G, B)
        idx <- idx - 1
      }
    }
  })
  
  return(list(A = A,B = B))
}

#For all d := np-p(p+1)/2 angles, get the derivative of each
#of the first nxp elements of the Givens matrix.
#The derivative w.r.t. a single angle, theta_ij, is computed very similarly
#to how the Givens matrix is computed. We multiply all the same
#rotation matrices as in the Givens matrix, except we replace
#the R_ij rotation matrix with its derivative w.r.t theta_ij, D.
#In other words,if A is the result of the multiplications before R_ij
#and if B is the result of multiplications after R_ij, we simply compute ADB.
#We will have the partial Givens matrices, A and B be precomputed
#and passed in.
#Again note that this is the derivative of the first nxp elements
#of the Givens matrix w.r.t a single angle. we of course need to do
#this for all d angles at which point we can stack the d nxp matrices
#to get an nxpxd tensor that we can slice to get p nxd Jacobians, which
#are the Jacobians of the first p columns of G w.r.t. all the angles. 
#Note we only need the Jacobians of the first p columns of G
#so the partial Givens, B on the right, need only consist of the
#first p columns of the partial Givens
GetGivensJacobians <- function(PartialRotationsAB, angles, n, p) {
  
  DerivativesList <- list()
  
  with(tf$name_scope(paste0('GivensJacobians')), {
    idx <- 0
    for(i in 0:(p-1)) {
      for(j in (i+1):(n-1)) {
        with(tf$name_scope(paste0('dG_dTheta_',i,j)), {
          D <- CreateDerivativeOfRotationMatrix(angles[idx], n, i, j)
          A <- PartialRotationsAB$A[[idx + 1]]
          B <- PartialRotationsAB$B[[idx + 1 + 1]]
          DB <- tf$matmul(D,B,a_is_sparse=TRUE, name = paste0('DB'))
          ADB <- tf$matmul(A, DB, name = paste0('ADB'))
          dG_dTheta_ij <- tf$identity(ADB, name = paste0('dG_dTheta'))
          DerivativesList <- c(DerivativesList, dG_dTheta_ij)
        })
        
        idx <- idx + 1
      }
    }
    
    GivensJacobianTensor <- tf$stack(DerivativesList, name = 'GivensJacobianTensor')
    
    StiefelDim <- as.integer(n*p - p*(p+1)/2)
    GivensJacobians <- list()
    for(i in 0:(p-1)) {
      #the much simpler syntax GivensJacobianTensor[,,i] returns an error from RTensorflow
      #so we must resort to manually doing the slice
      J <- tf$transpose(tf$reshape(
        tf$slice(GivensJacobianTensor, list(0L, 0L, as.integer(i)), list(StiefelDim, as.integer(n), 1L)),
        list(StiefelDim, as.integer(n))),
        name = paste0('JacobianColumn', i))
      GivensJacobians <- c(GivensJacobians, J)
    }
  })
  
  return(GivensJacobians)
}

GetStiefelAreaForm <- function(G, GivensJacobians, n, p) {
  
  FList <- list()
  StiefelDim <- as.integer(n*p - p*(p+1)/2)
  
  with(tf$name_scope('StiefelAreaForm'), {
    for(i in 0:(p-1)) {
      with(tf$name_scope(paste0('OneFormsColumn', i)), {
        #Rows <- tf$transpose(G[,(i+1):(n-1)])
        GTransposeRows <- tf$slice(tf$transpose(G), list(as.integer(i+1), 0L), list(as.integer(n-i-1), as.integer(n)),name = 'GTransposeRows')
        #Jacobian <- GetJacobian(G[,i], angles)
        OneForms <- tf$transpose(tf$matmul(GTransposeRows, GivensJacobians[[i+1]]), name = 'OneForms')
        
        for(j in 0:(n-i-1-1)) {
          OneForm <- tf$slice(OneForms, list(0L, as.integer(j)), list(as.integer(StiefelDim), 1L), name = paste0('OneForm',j))
          FList <- c(FList, OneForm)
        }
      })
    }
    
    FMatrix <- tf$concat(FList, axis = 1L, name = 'F')
    det <- tf$matrix_determinant(FMatrix, name = 'detF')
  })
  
  return(det)
}

CreateThetaConstrained <- function(ThetaUnconstrained, n, p) {
  
  ThetaConstrainedList <- list()
  ThetaConstrainedDerivativeList <- list()
  
  with(tf$name_scope('ThetaConstrained'), {
    Pi <- tf$constant(pi, name = 'Pi')
    Pi2 <- tf$constant(pi/2, name = 'Pi/2')
    idx <- 0
    for(i in 0:(p-1)) {
      for(j in (i+1):(n-1)) {
          #first rotation of the column should go from -pi to pi
          # if(j == i+1) {
          #   a <- -Pi
          #   b <- Pi
          # }
          # else {
          #   a <- -Pi2
          #   b <- Pi2
          # }
        a <- -Pi2
        b <- Pi2
        
        with(tf$name_scope(paste0('ThetaConstrained',i,j)), {
          ThetaConstrained <- a + (b-a)*tf$sigmoid(ThetaUnconstrained[idx])
          ThetaConstrainedList <- c(ThetaConstrainedList, ThetaConstrained)
        })
          
        with(tf$name_scope(paste0('ThetaConstrainedDerivative',i,j)), {
          ThetaConstrainedDerivative <- (b-a)*tf$sigmoid(ThetaUnconstrained[idx])*(1-tf$sigmoid(ThetaUnconstrained[idx]))
          ThetaConstrainedDerivativeList <- c(ThetaConstrainedDerivativeList, ThetaConstrainedDerivative)
            
          idx <- idx + 1
        })
      }
    }
    
    ThetaConstrained <- tf$stack(ThetaConstrainedList, name = 'ThetaConstrained')
    ThetaConstrainedDerivative <- tf$stack(ThetaConstrainedDerivativeList, name = 'ThetaConstrainedDerivative')
  })
  
  return(list(ThetaConstrained = ThetaConstrained, ThetaConstrainedDerivative = ThetaConstrainedDerivative))
}

CreateLambdaConstrained <- function(LambdaUnconstrained, p) {
  
  LambdaConstrainedList <- list()

  #transform for ordered vector in section VI of Stan manual
  with(tf$name_scope('LambdaConstrained'), {

    ExpLambdaUnconstrained <- tf$exp(LambdaUnconstrained)
    
    for(i in (p-1):0) {
      with(tf$name_scope(paste0('LambdaConstrained',i)), {
        if(i == p-1) LambdaConstrained <- LambdaUnconstrained[i]
        #i index should be i -1 but R does from 1 indexing
        else LambdaConstrained <- LambdaConstrainedList[[1]] + ExpLambdaUnconstrained[i]
        
        #add to the front of the list because last Lambda will be highest
        LambdaConstrainedList <- c(LambdaConstrained, LambdaConstrainedList)
      })
    }

    LambdaConstrained <- tf$stack(LambdaConstrainedList, name = 'LambdaConstrained')
    LambdaConstrainedDerivative <- tf$reduce_prod(ExpLambdaUnconstrained, name = 'LambdaConstrainedDerivative')
  })
  
  return(list(LambdaConstrained = LambdaConstrained, LambdaConstrainedDerivative = LambdaConstrainedDerivative))
}

