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
  
  vector set_theta(int n, int p, vector theta_lon, vector theta_lat) {

    int idx = 1;
    int idx_lat = 1;
    int d = n*p - (p*(p+1))/2;
    vector[d] theta;
    
    // If n==p then we shouldn't do the last column will have no degrees of freedom
    // thus there is no angle corresponding to htat las column.
    int pp = p;
    if(p == n) pp = p-1;
    for(i in 1:pp) {
      
      theta[idx] = theta_lon[i];
      
      if(i <= n-2) { 
        theta[(idx+1):(idx+n-i-1)] = theta_lat[idx_lat:(idx_lat+n-i-1-1)];
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
      
      int i_rev = pp-i+1; // create a reverse index that starts from p and goes down to 1
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
          //target += (j_rev-i_rev-1)*log(cos(theta_ij));
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
  
  int N;  
  matrix[n,n] SigmaHat;
  
  matrix[n, p] Y_ml;
  
  real tau0;
  real nu;
  real s;
}
transformed data {
  real slab_scale = s;    // Scale for large slopes
  real slab_scale2 = square(slab_scale);
  real slab_df = nu;      // Effective degrees of freedom for large slopes
  real half_slab_df = 0.5 * slab_df;
  
  int d = n*p - (p*(p+1))/2;
  
  // there is one longitudinal angle for every column
  // except in the case where n == p in which case the
  // last column has no degrees of freedom
  int pp = p;
  if(n == p) pp = p - 1;
}
parameters {
  vector[pp] x_lon;
  vector[pp] y_lon;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d-pp] theta_lat;
  
  positive_ordered[p] lambda_sq_rev;
  real<lower = 0> sigmaSq;
  
  vector<lower=0>[d] lambda;
  real<lower=0> c2_tilde;
  real<lower=0> tau_tilde;
}
transformed parameters{
  vector[pp] theta_lon;
  vector[d] theta_tilde;
  vector[pp] r;
  matrix[n,p] Y;
  vector<lower=0>[p] lambda_sq;

  vector[d] theta;

  theta_lon = atan2_vec(x_lon, y_lon);
  r = hypot_vec(x_lon, y_lon);
  theta_tilde = set_theta(n, p, theta_lon, theta_lat);
  
  {
    real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)

    // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
    // Implies that marginally beta ~ student_t(slab_df, 0, slab_scale)
    real c2 = slab_scale2 * c2_tilde;

    vector[d] lambda_tilde =
      sqrt( c2 * square(lambda) ./ (c2 + square(tau) * square(lambda)) );

    // beta ~ normal(0, tau * lambda_tilde)
    theta = tau * lambda_tilde .* theta_tilde;
  }
  
  // reverse lambda
  for(i in 1:p) lambda_sq[i] = lambda_sq_rev[p - i + 1];
  
  Y = givens_lp(n,p,theta);
}
model {
  matrix[n,n] C;
  r ~ normal(1.0, 0.1);
  
  theta_tilde ~ normal(0, 1);
  lambda ~ cauchy(0, 1);
  tau_tilde ~ cauchy(0, 1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
  
  C = Y*diag_matrix(lambda_sq)*Y' + diag_matrix(rep_vector(sigmaSq, n));
  target += -(N/2)*log(determinant(C)) - (N/2)*trace(C\SigmaHat);
}
generated quantities {
  // compute principal angles between columns of ML estimate
  vector[p] theta_princ;
  real Y_sparsity = 0.0;
  for(j in 1:p) {
    real qTv = abs(dot_product(Y[,j], Y_ml[,j]));
    real q = sqrt(dot_self(Y[,j]));
    real v = sqrt(dot_self(Y_ml[,j]));
    theta_princ[j] = acos(qTv/(q*v));
  }
  
  for(j in 1:p) {
    for(i in 1:n) {
      if(Y[i,j] > 0.01) {
        Y_sparsity += 1.0;
      }
    }
  }
  Y_sparsity = Y_sparsity/(n*p);
}
