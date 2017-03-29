library(tensorflow)
sess <- tf$Session()

#assume q is 5 dimensional
U_R <- function(q) {
  return(sess$run(U, feed_dict = dict(Theta01 = q[1], Theta02 = q[2], Theta12 = q[3], LambdaVec = list(q[4],q[5]), SigmaSq = 1)))
}

GradU_R <- function(q) {
  return(unlist(sess$run(GradU, feed_dict = dict(Theta01 = q[1], Theta02 = q[2], Theta12 = q[3], LambdaVec = list(q[4],q[5]), SigmaSq = 1))))
} 


HMC <- function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q),0,1)
  current_p = p
  
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if(runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
    return (q) # accept
  else
    return (current_q) # reject
}

HMC(U_R, GradU_R, epsilon = 0.001, L = 10, current_q = c(pi/4, pi/4, 0, 5, 3))

NumSamples <- 1000
q <- c(pi/4, pi/4, 0, 5, 3)
Samples2 <- matrix(nrow = NumSamples+1, ncol = 5)
Samples2[1,] <- q
system.time(
for(s in 1:NumSamples) {
  q <- HMC(U_R, GradU_R, epsilon = 0.001, L = 50, current_q = q)
  Samples2[s+1,] <- q
}
)
Samples2 <- data.frame(Samples2)