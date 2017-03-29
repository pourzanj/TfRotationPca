library(plotly)

N <- 1000
p <- 2
z <- matrix(rnorm(N*p), nrow = N)

W <- matrix(c(1/2, -1/sqrt(2),
              1/2, 1/sqrt(2),
              1/sqrt(2), 0), nrow = 3, byrow = TRUE)

x <- t(W %*% matrix(c(5, 0, 0, 3), nrow = 2) %*% t(z)) + matrix(rnorm(N*3, sd = 1), nrow = N)
df <- data.frame(x[1:1000,])

qplot(z[,1], z[,2])

plot_ly() %>%
  add_trace(data = df, x = X1, y = X2, z = X3, mode = "markers", type = "scatter3d", 
            marker = list(size = 5, color = "red", symbol = 104))

SigmaHat_ <- (1/N)*t(x) %*% x
