functions {
    
  vector atan2_vec(vector x, vector y) {
    
    int d = num_elements(x);
    vector[d] theta;
    
    for(i in 1:d) {
      theta[i] = atan2(y[i], x[i]);  
    }
    
    return(theta);
  }
  
  vector hypot_vec(vector x, vector y) {
    
    int d = num_elements(x);
    vector[d] r;
    
    for(i in 1:d) {
      r[i] = hypot(x[i], y[i]);  
    }
    
    return(r);
  }
  
  vector mirror(int n, int p, vector theta_lon, vector theta_lat) {

    int idx = 1;
    int idx_lat = 1;
    int d = n*p - (p*(p+1))/2;
    vector[d] theta;
    
    // If n==p then we shouldn't do the last column will have no degrees of freedom
    // thus there is no angle corresponding to htat las column.
    int pp = p;
    if(p == n) pp = p-1;
    for(i in 1:pp) {
      
      // check the longitudinal angle to see if we should mirror
      int mirror;  
      real eps;
      
      if(theta_lon[i] > 0.0) {
        // upper plane
        eps = theta_lon[i] - pi()/2;
        if(eps <= 0.0) {
          theta[idx] = theta_lon[i];
          mirror = 0;
        }
        else {
          theta[idx] = -pi()/2 + eps;
          mirror = 1;
        }
      }
      
      else {
        // lower plane
        eps = theta_lon[i] + pi()/2;
        if(eps >= 0.0) {
          theta[idx] = theta_lon[i];
          mirror = 0;
        }
        else {
          theta[idx] = pi()/2 + eps;
          mirror = 1;
        }
      }
     
     if(i <= n-2) { 
      // set lat angles according to whether they're mirrored or not
      if(mirror == 1) {
        theta[(idx+1):(idx+n-i-1)] = -theta_lat[idx_lat:(idx_lat+n-i-1-1)];
      } else {
        theta[(idx+1):(idx+n-i-1)] = theta_lat[idx_lat:(idx_lat+n-i-1-1)];
      }
     }
      
      idx = idx + n-i;
      idx_lat = idx_lat + n-i-1;
    }
    
    return(theta);
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
    // If n==p then we shouldn't do the last column will have no degrees of freedom
    // thus there is no angle corresponding to htat las column.
    int pp = p;
    if(p == n) pp = p-1;
    for(i in 1:pp) {
      
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
        // only update target for latiduninal angles. otherwise we might
        // get a cosine that is negative and can't take log
        if(j_rev > (i_rev+1) ) {
          target += (j_rev-i_rev-1)*log(cos(theta_ij));
        }
        
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
    int N; // number of data points
    
    int<lower=1> K; // num categories
    int<lower=0> interactions[N,n,n];
}

transformed data {
  int d = n*p - (p*(p+1))/2;
}

parameters {
  
  //transit probs
  simplex[K] phi[K];
  
  // emission parameters
  vector[p] x_lon[K]; vector[p] y_lon[K];
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d-p] theta_lat[K];
    
  positive_ordered[p] lambda_rev[K];
  
  // hierarchical parameter on transition probs
  vector<lower=0>[K] beta;
}

transformed parameters{
  
  vector[p] theta_lon[K];
  vector[d] theta[K];
  vector[p] r[K];
  
  matrix[n, p] W[K];
  matrix[n, n] R[K];
    
  vector<lower=0>[p] lambda_sq[K];
    
  for (k in 1:K) {
    
    // set W
    theta_lon[k] = atan2_vec(x_lon[k], y_lon[k]);
    r[k] = hypot_vec(x_lon[k], y_lon[k]);
    theta[k] = mirror(n, p, theta_lon[k], theta_lat[k]);
    
    W[k] = givens_lp(n,p,theta[k]);
    
    // reverse lambda
    for(i in 1:p) lambda_sq[k][i] = lambda_rev[k][p - i + 1]^2;
    
    // set rates of interactions
    R[k] = W[k]*diag_matrix(lambda_sq[k])*W[k]';
  }
}

model {
  real acc[K];
  real gamma[N, K];
  
  // put dirichlet priors on transition probs
  for (k in 1:K) phi[k] ~ dirichlet(beta);
    
  // likelihood of interactions during first window
  for (k in 1:K) {
    gamma[1,k] = 0;
    for(x in 1:n) {
      for(y in x:n) {
        gamma[1,k] = gamma[1,k] + poisson_lpmf(interactions[1,x,y] | exp(R[k][x,y]));
      }
    }
  }
    
  // likelihood for rest of time windows
  for (t in 2:N) {
    for (k in 1:K) {
      for (j in 1:K) {
        
        acc[j] = gamma[t-1,j] + log(phi[j,k]);
        for(x in 1:n) {
          for(y in x:n) {
            acc[j] = acc[j] + poisson_lpmf(interactions[t,x,y] | exp(R[k][x,y]));
          }
        }
      }
      gamma[t,k] = log_sum_exp(acc);
    }
  }
    
  // accumulate likelihood
  target += log_sum_exp(gamma[N]);
}

generated quantities {
    int<lower=1, upper=K> y_star[N];
    real log_p_y_star;
    {
        int back_ptr[N, K];
        real best_logp[N, K];
        for (k in 1:K) {
            best_logp[1,k] = 0;
            for(x in 1:n) {
                for(y in x:n) {
                    best_logp[1,k] = best_logp[1,k] + poisson_lpmf(interactions[1,x,y] | exp(R[k][x,y]));
                }
            }
        }
    
        for (t in 2:N) {
            for (k in 1:K) {
                best_logp[t,k] = negative_infinity();
                for(j in 1:K) {
                    real logp;
                    logp = best_logp[t-1,j] + log(phi[j,k]);
                    for(x in 1:n) {
                        for(y in x:n) {
                            logp = logp + poisson_lpmf(interactions[t,x,y] | exp(R[k][x,y]));
                        }
                    }
                    if (logp > best_logp[t,k]) {
                        back_ptr[t,k] = j;
                        best_logp[t,k] = logp;
                    }

                }
            }
        }

        log_p_y_star = max(best_logp[N]);
        for(k in 1:K) {
            if(best_logp[N,k] == log_p_y_star) {
                y_star[N] = k;
            }
        }

        for(t in 1:(N-1)){
            y_star[N - t] = back_ptr[N - t + 1, y_star[N - t + 1]];
        }
    } 
}
                
