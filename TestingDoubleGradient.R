library(tensorflow)

sess <- tf$Session()

Theta01 <- tf$placeholder(tf$float32)
Theta02 <- tf$placeholder(tf$float32)
Theta12 <- tf$placeholder(tf$float32)

R01 <- CreateRotationMatrix(Theta01, 3, 0, 1)
R02 <- CreateRotationMatrix(Theta02, 3, 0, 2)
R12 <- CreateRotationMatrix(Theta12, 3, 1, 2)
G <- tf$matmul(R01, R02)
G <- tf$matmul(G, R12)

Grad <- tf$gradients(G[0,0], Theta01)
GradStacked <- tf$stack(Grad)
Grad2 <- tf$gradients(GradStacked, Theta01)

sess$run(Grad2, feed_dict = dict(Theta01 = 0, Theta02 = 0, Theta12 = 0))
############################

G <- CreateGivensMatrix(list(Theta01, Theta02, Theta12), n = 3, p = 2)
Grad <- tf$gradients(G[0,0], Theta01)
GradStacked <- tf$stack(Grad)
Grad2 <- tf$gradients(GradStacked, Theta01)

sess$run(Grad2, feed_dict = dict(Theta01 = 0, Theta02 = 0, Theta12 = 0))
