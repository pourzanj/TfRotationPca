functions {
  
  real mirror_atan2(real theta) {
    
    real eps;
    
    // upper plane
    if(theta > 0.0) {
      eps = theta - pi()/2;
      if(eps <= 0.0)
        return(theta);
      else
        return(-pi()/2 + eps);
    }
    // lower plane
    else {
      eps = theta + pi()/2;
      if(eps >= 0.0)
        return(theta);
      else
        return(pi()/2 + eps);
    }
  }
  
  // return Y = R_12 ... R_1n R_23 ... R_2n ... R_pp+1 ... Rpn I_np
  // where I_np is the first p columns of the n x n identity matrix
  matrix givens_lp(int n, int p, vector theta) {

    int d = n*p - p*(p+1)/2;

    // set Y to I_np so we can apply left rotations
    matrix[n,p] Y = diag_matrix(rep_vector(1, n))[,1:p];

    int idx = d; // index to keep track of which angle we're on
    
    // fill in reverse order from right to left. This way we only have to
    // store p colums at a time instead of n
    for(i in 1:p) {
      
      int i_rev = p-i+1; // create a reverse index that starts from p and goes down to 1
      
      for(j in (i_rev+1):n) {
        
        // index goes from p+1 to n, so reverse it to go from n to p+1
        int j_rev = n-j+(i_rev+1);
        
        // apply left rotation which affects rows of Y. we must keep updated 
        // rows in temp variable as oppose to updating Y in place 
        real theta_ij = theta[idx];
        row_vector[p] Y_i_temp = cos(theta_ij)*Y[i_rev,] - sin(theta_ij)*Y[j_rev,];
        row_vector[p] Y_j_temp = sin(theta_ij)*Y[i_rev,] + cos(theta_ij)*Y[j_rev,];
        Y[i_rev,] = Y_i_temp;
        Y[j_rev,] = Y_j_temp;
        
        // update jacobian determinant with this angle
        target += (j_rev-i_rev-1)*log(cos(theta_ij));
        
        // go to the next angle to the left
        idx = idx - 1;
      }
    }

    return(Y);
  }
}

data {
  int n;
  int p;
  
  int N;  
  matrix[n,n] SigmaHat;
  real sigmaSqHyperPrior;
}
transformed data {
  int d = n*p - p*(p+1)/2;
}
parameters {
  vector[d] x;
  vector[d] y;
  
  positive_ordered[p] lambdaReversed;
  real<lower = 0> sigmaSq;
}
transformed parameters{
  vector[d] theta;
  vector[d] r;
  matrix[n,p] Y;
  
  vector<lower=0>[p] lambdaSq;

  // transform donut
  for(i in 1:d) {
    theta[i] = mirror_atan2(atan2(y[i], x[i]));
    r[i] = hypot(x[i],y[i]); 
  }
  
  // reverse lambda
  for(i in 1:p) lambdaSq[i] = lambdaReversed[p - i + 1]^2;
  
  Y = givens_lp(n,p,theta);
}
model {
  matrix[n,n] C;
  r ~ normal(1.0, 0.1);
  
  //sigmaSq ~ normal(1, sigmaSqHyperPrior);
  C = Y*diag_matrix(lambdaSq)*Y' + diag_matrix(rep_vector(sigmaSq, n));
  target += -(N/2)*log(determinant(C)) - (N/2)*trace(C\SigmaHat);
}
