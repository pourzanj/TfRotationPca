library(dplyr)
library(plotly)

GenerateHighDimData <- function(n, p, W, LambdaVec, sd = 1, N) {
  z <- matrix(rnorm(N*p), nrow = N)
  W_Lambda_z <- t(W %*% diag(LambdaVec, nrow = p, ncol = p) %*% t(z))
  eps <- matrix(rnorm(N*n, sd = sd), nrow = N)
  
  x <- W_Lambda_z + eps %>% tbl_df()
  
  return(list(z = tbl_df(z), x = x))
}

# W <- matrix(c(1/2, -1/sqrt(2),
#               1/2, 1/sqrt(2),
#               1/sqrt(2), 0), nrow = 3, byrow = TRUE)
# 
# ThreeTwo100 <- GenerateHighDimData(n = 3, p = 2, W = W, LambdaVec = c(5, 3), sd = 1, N = 100)
# z <- ThreeTwo100$z %>% tbl_df
# x <- ThreeTwo100$x %>% tbl_df
# qplot(V1, V2, data = z)
# plot_ly(x, x = ~V1, y = ~V2, z = ~V3) %>% add_markers()

# WA <- InverseGivensTransform(c(-pi/4,pi/4,0), 3, 2)
# WB <- InverseGivensTransform(c(-pi/4,-pi/4,0), 3, 2)
# WC <- InverseGivensTransform(c(pi/4,pi/4,0), 3, 2)
# 
# N <- 20
# ThreeTwoA <- GenerateHighDimData(n = 3, p = 2, W = WA, LambdaVec = c(5, 3), sd = 1.0, N = N)
# ThreeTwoB <- GenerateHighDimData(n = 3, p = 2, W = WB, LambdaVec = c(5, 3), sd = 1.0, N = N)
# ThreeTwoC <- GenerateHighDimData(n = 3, p = 2, W = WC, LambdaVec = c(5, 3), sd = 0.5, N = 100)
# 
# xA <- ThreeTwoA$x %>% tbl_df %>% mutate(Group = "A")
# xB <- ThreeTwoB$x %>% tbl_df %>% mutate(Group = "B")
# xC <- ThreeTwoC$x %>% tbl_df %>% mutate(Group = "C")
# 
# x <- rbind(xA, xB) %>% mutate(Group = as.factor(Group))
# N <- 40
# 
# plot_ly(x, x = ~V1, y = ~V2, z = ~V3, color = ~Group, colors = c('Red', 'Orange')) %>%
#   add_markers()

#   
# plot_ly(showscale = FALSE) %>% add_surface(z = ~z, opacity = 0.9) %>%
#   add_trace(data = x, x = ~V1, y = ~V2, z = ~V3, mode = "markers", type = "scatter3d", 
#             marker = list(size = 5, color = "red", symbol = 104))
# 
# plot_ly() %>%
#   add_trace(data = ThreeDim100, x = X1, y = X2, z = X3, mode = "markers", type = "scatter3d", 
#             marker = list(size = 5, color = "red", symbol = 104))

# xMat <- as.matrix(xB[,1:3])
# SigmaHat_ <- (1/N)*t(xMat) %*% xMat
