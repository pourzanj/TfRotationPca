data {
  int n;
  int p;
  
  int N;  
  matrix[n,n] SigmaHat;
}
parameters {
  matrix[n,p] Y;
  positive_ordered[p] lambdaReversed;
  real<lower = 0> sigmaSq;
}
transformed parameters{
  vector<lower=0>[p] lambdaSq;

  // reverse lambda
  for(i in 1:p) lambdaSq[i] = lambdaReversed[p - i + 1]^2;
}
model {
  matrix[n,n] C = Y*diag_matrix(lambdaSq)*Y' + diag_matrix(rep_vector(sigmaSq, n));
  target += -(N/2)*log(determinant(C)) - (N/2)*trace(C\SigmaHat);
}
