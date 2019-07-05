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
  
  matrix givens_lp(int N, int P, vector theta) {

    int D = N*P - P*(P+1)/2;
    int d = 1; // index to track which angle we're on
    
    matrix[N,P] Y = diag_matrix(rep_vector(1, N))[,1:P];
    vector[D] cos_theta_power;
    
    int idx = 1;
    
    // fill in reverse order
    for(j in 1:P) {
      int p_rev = P-j+1;
      for(i in (p_rev+1):N) {
        
        int n_rev = N-i+p_rev+1;
        row_vector[P] Y_i_temp = cos(theta[d])*Y[p_rev,] - sin(theta[d])*Y[n_rev,];
        row_vector[P] Y_j_temp = sin(theta[d])*Y[p_rev,] + cos(theta[d])*Y[n_rev,];
        Y[p_rev,] = Y_i_temp;
        Y[n_rev,] = Y_j_temp;
        
        d = d + 1;
      }
    }
    
    for(i in 1:P) {
      for(j in (i+1):N) {
        cos_theta_power[idx] = cos(theta[idx])^(j-i-1);
        idx = idx + 1;
    }
  }
    
    target += sum(log(cos_theta_power));

    return(Y);
  }
}

data {
  int N;
  int P;
}
transformed data {
  int D = N*P - P*(P+1)/2;
}
parameters {
  vector[D] x;
  vector[D] y;
}
transformed parameters{
  vector[D] theta_raw;
  vector[D] theta;
  vector[D] r;
  matrix[N,P] Y;

  // transform donut
  for(d in 1:D) {
    theta_raw[d] = atan2(y[d], x[d]);
    theta[d] = mirror_atan2(theta_raw[d]);
    r[d] = hypot(x[d],y[d]); 
  }
  
  Y = givens_lp(N,P,theta);
}
model {
  r ~ normal(1.0, 0.1);

}
