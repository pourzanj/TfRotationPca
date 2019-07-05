#requires x and y coordinates to find max an min
GetSurfacePoints <- function(x, y, W) {
  
  xMin <- floor(min(x)) - 1
  xMax <- ceiling(max(x)) + 1
  yMin <- floor(min(y)) - 1
  yMax <- ceiling(max(y)) + 1
  
  x <- c(xMin, xMax, xMin, xMax)
  y <- c(yMin, yMin, yMax, yMax)  
  xyPoints <- matrix(c(x,y), 2, 4, byrow = TRUE)
  ab <- solve(W[1:2,], xyPoints)
  z <- (W %*% ab)[3,]
  
  dim(x) <- c(2,2)
  dim(y) <- c(2,2)
  dim(z) <- c(2,2)
  
  return(list(x = x, y = y, z = z))
}