library(rstan)

n_0 <-
p_0 <-
d_W_0 <-
d_B_0 <-

n_1 <-
p_1 <-
d_W_1 <-
d_B_1 <-

n_2 <-
p_2 <-
d_W_2 <-
d_B_2 <-

n_3 <-
p_3 <-
d_W_3 <-
d_B_3 <-

n_4 <-
p_4 <-
d_W_4 <-
d_B_4 <-

p_common <-

N <-
SigmaHat <- 

data <- list(n_0, p_0, d_W_0, d_B_0, n_1, p_1, d_W_1, d_W_1, 
             n_2, p_2, d_W_2, d_B_2, n_3, p_3, d_W_3, d_B_3, 
             n_4, p_4, d_W_4, d_B_4, N, SigmaHat = SigmaHat)
fit <- stan(file = "Stan/CcaGen.stan", data = data, chains = 4, iter = 10000)

