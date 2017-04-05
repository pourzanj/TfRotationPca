library(purrr)

SampleStiefelUniformQr <- function(n, p) {
  RandomNormalMatrix <- matrix(rnorm(n*p), nrow = n, ncol = p)
  QR <- qr(RandomNormalMatrix)
  
  #Q will be the same sign for v and -v, but R will handles the sign
  Q <- qr.Q(QR)*sign(diag(qr.R(QR)))
  
  return(RandomNormalMatrix)
}

GetSamplesStiefelUniormQr <- function(n, p, NumSamples) {
  Samples <- map(1:NumSamples, function(x) SampleStiefelUniformQr(n,p)) %>%
    map(GivensTransform)
  
  SamplesDf <- do.call(rbind, Samples) %>% tbl_df
  
  return(SamplesDf)
}