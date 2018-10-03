library(tensorflow)

# n <- 10
# p <- 8
# d <- n*p-p*(p+1)/2

#generate quick synthetic data (for testing)
# N <- 1000
# W <- InverseGivensTransform(rep(0,d),n,p)
# X <- GenerateHighDimData(n, p, W, rep(1,p), 1, 1000)$x %>% as.matrix
# SigmaHat_ <- (1/N)*t(X) %*% X

sess <- tf$Session()
ThetaUnconstrained <- tf$placeholder(tf$float32, shape = shape(d), name = 'ThetaUnconstrained')
#p-1 rather than p when n==p because of a bug in CreateTheta Unconstrained where when p==n
#extra angle variables are created in the loop that creates the angles
ThetaConstrained <- CreateThetaConstrained(ThetaUnconstrained, n, ifelse(p == n, p-1, p))
ThetaConstrainedDerivative <- ThetaConstrained$ThetaConstrainedDerivative
ThetaConstrained <- ThetaConstrained$ThetaConstrained

LambdaUnconstrainedVec <- tf$placeholder(tf$float32, shape(p), 'LambdaUnconstrained')
#LambdaVec <- tf$placeholder(tf$float32, shape(p))
LambdaConstrained <- CreateLambdaConstrained(LambdaUnconstrainedVec, p)
LambdaConstrainedDerivative <- LambdaConstrained$LambdaConstrainedDerivative
LambdaConstrainedVec <- LambdaConstrained$LambdaConstrained
Lambda <- tf$diag(LambdaConstrainedVec, name = 'Lambda')
LambdaSq <- Lambda^2

LogSigmaSq <- tf$placeholder(tf$float32, name = 'LogSigmaSq')
Id <- tf$constant(diag(n),dtype = tf$float32, name = 'Id_n')

PartialRotations <- CreatePartialGivens(ThetaConstrained, n, ifelse(p == n, p-1, p))
G <- CreateGivensMatrix(PartialRotations, n, p)
W <- G[,0:(p-1)]
GivensJacobians <- GetGivensJacobians(PartialRotations, ThetaConstrained, n, ifelse(p == n, p-1, p))
LogStiefelAreaForm <- tf$log(GetStiefelAreaForm(G, GivensJacobians, n, ifelse(p == n, p-1, p)))

C <- tf$matmul(W, tf$matmul(LambdaSq, tf$transpose(W))) + tf$exp(LogSigmaSq)*Id
SigmaHat <- tf$constant(SigmaHat_, dtype = tf$float32)

U <- (N/2)*log(tf$matrix_determinant(C)) + (N/2)*tf$trace(tf$matrix_solve(C, SigmaHat)) - LogStiefelAreaForm - tf$reduce_sum(tf$log(ThetaConstrainedDerivative)) - tf$log(LambdaConstrainedDerivative) - LogSigmaSq
GradU <- tf$gradients(U, list(ThetaUnconstrained, LambdaUnconstrainedVec, LogSigmaSq))

#for testing uniform on the Stiefel
#U <- - LogStiefelAreaForm - tf$reduce_sum(tf$log(ThetaConstrainedDerivative))
#GradU <- tf$gradients(U, ThetaUnconstrained)

#writer <- tf$summary$FileWriter("./TfLogs", sess$graph)
#sess$run(GradU, feed_dict = dict(ThetaUnconstrained = c(pi/2,0,0), LambdaVec = c(5,3), SigmaSq = 1))

# system.time(
# for(i in 1:100) sess$run(grad, feed_dict = dict(Theta = rep(0, d)))
# )


# sess <- tf$Session()
# 
# Theta01 <- tf$placeholder(tf$float32, shape = shape())
# Theta02 <- tf$placeholder(tf$float32, shape = shape())
# Theta12 <- tf$placeholder(tf$float32, shape = shape())
# Theta <- tf$stack(list(Theta01, Theta02, Theta12))
# 
# n <- 4
# p <- 2
# StiefelDim <- n*p - p*(p+1)/2
# Theta <- tf$placeholder((tf$float32), shape = shape(StiefelDim))
# 
# LambdaVec <- tf$placeholder(tf$float32, shape(p))
# Lambda <- tf$diag(LambdaVec)
# LambdaSq <- Lambda^2
# 
# SigmaSq <- tf$placeholder(tf$float32)
# Id <- tf$constant(diag(n),dtype = tf$float32)
# 
# G <- CreateGivensMatrix(Theta, n = n, p = p)
# W <- G[,0:(p-1)]
# 
# C <- tf$matmul(W, tf$matmul(LambdaSq, tf$transpose(W))) + SigmaSq*Id
# SigmaHat <- tf$constant(SigmaHat_, dtype = tf$float32)
# Area <- GetStiefelAreaForm(G, Theta, n = n, p = p)
# 
# U <- (1000/2)*log(tf$matrix_determinant(C)) + (1000/2)*tf$trace(tf$matrix_solve(C, SigmaHat)) - log(Area)
# GradU <- tf$gradients(U, list(Theta, LambdaVec))
# 
# GradLogArea <- tf$gradients(log(Area), list(Theta))
# system.time(for(i in 1:10)sess$run(GradLogArea, feed_dict = dict(Theta = list(pi/4, pi/4, 0, 0, 0), LambdaVec = list(5,3), SigmaSq = 1)))
