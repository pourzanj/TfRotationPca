library(tensorflow)

sess <- tf$Session()

Theta01 <- tf$placeholder(tf$float32, shape = shape())
Theta02 <- tf$placeholder(tf$float32, shape = shape())
Theta12 <- tf$placeholder(tf$float32, shape = shape())
Theta <- tf$stack(list(Theta01, Theta02, Theta12))

n <- 4
p <- 2
StiefelDim <- n*p - p*(p+1)/2
Theta <- tf$placeholder((tf$float32), shape = shape(StiefelDim))

LambdaVec <- tf$placeholder(tf$float32, shape(p))
Lambda <- tf$diag(LambdaVec)
LambdaSq <- Lambda^2

SigmaSq <- tf$placeholder(tf$float32)
Id <- tf$constant(diag(n),dtype = tf$float32)

G <- CreateGivensMatrix(Theta, n = n, p = p)
W <- G[,0:(p-1)]

C <- tf$matmul(W, tf$matmul(LambdaSq, tf$transpose(W))) + SigmaSq*Id
SigmaHat <- tf$constant(SigmaHat_, dtype = tf$float32)
Area <- GetStiefelAreaForm(G, Theta, n = n, p = p)

U <- (1000/2)*log(tf$matrix_determinant(C)) + (1000/2)*tf$trace(tf$matrix_solve(C, SigmaHat)) - log(Area)
GradU <- tf$gradients(U, list(Theta, LambdaVec))

GradLogArea <- tf$gradients(log(Area), list(Theta))
system.time(for(i in 1:10)sess$run(GradLogArea, feed_dict = dict(Theta = list(pi/4, pi/4, 0, 0, 0), LambdaVec = list(5,3), SigmaSq = 1)))
