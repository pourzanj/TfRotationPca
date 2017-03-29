library(tensorflow)

sess <- tf$Session()

Theta01 <- tf$placeholder(tf$float32, shape = shape())
Theta02 <- tf$placeholder(tf$float32, shape = shape())
Theta12 <- tf$placeholder(tf$float32, shape = shape())
Theta <- tf$stack(list(Theta01, Theta02, Theta12))

LambdaVec <- tf$placeholder(tf$float32, shape(2))
Lambda <- tf$diag(LambdaVec)
LambdaSq <- Lambda^2

SigmaSq <- tf$placeholder(tf$float32)
Id <- tf$constant(diag(3),dtype = tf$float32)

G <- CreateGivensMatrix(Theta, n = 3, p = 2)
W <- G[,0:1]

C <- tf$matmul(W, tf$matmul(LambdaSq, tf$transpose(W))) + SigmaSq*Id
SigmaHat <- tf$constant(SigmaHat_, dtype = tf$float32)
Area <- GetStiefelAreaForm(G, Theta, n = 3, p = 2)

U <- (1000/2)*log(tf$matrix_determinant(C)) + (1000/2)*tf$trace(tf$matrix_solve(C, SigmaHat)) - log(Area)
GradU <- tf$gradients(U, list(Theta01, Theta02, Theta12, LambdaVec))

#system.time(for(i in 1:1000)sess$run(GradU, feed_dict = dict(Theta01 = pi/4, Theta02 = pi/4, Theta12 = 0, LambdaVec = list(5,3), SigmaSq = 1)))

