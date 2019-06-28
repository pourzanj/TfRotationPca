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
    int d = n*p - p*(p+1)/2;
    vector[d] theta;
    
    for(i in 1:p) {
      
      // check the longitudinal angle to see if we should mirror
      int mirror;  
      real eps;
      
      if(theta_lon[p] > 0.0) {
        // upper plane
        eps = theta_lon[p] - pi()/2;
        if(eps <= 0.0) {
          theta[idx] = theta_lon[p];
          mirror = 0;
        }
        else {
          theta[idx] = -pi()/2 + eps;
          mirror = 1;
        }
      }
      
      else {
        // lower plane
        eps = theta_lon[p] + pi()/2;
        if(eps >= 0.0) {
          theta[idx] = theta_lon[p];
          mirror = 0;
        }
        else {
          theta[idx] = pi()/2 + eps;
          mirror = 1;
        }
      }
      
      // set lat angles according to whether they're mirrored or not
      if(mirror == 1) {
        theta[(idx+1):(idx+n-i-1)] = -theta_lat[idx_lat:(idx_lat+n-i-1-1)];
      } else {
        theta[(idx+1):(idx+n-i-1)] = theta_lat[idx_lat:(idx_lat+n-i-1-1)];
      }
      
      idx = idx + n-i;
      idx_lat = idx_lat + n-i-1;
    }
    
    return(theta);
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
}
transformed data {
  int d = n*p - p*(p+1)/2;
  
  // if p ==n the pth column does not have an angle in it
  int pp = p;
  if(p == n) pp = p-1;
}
parameters {
  vector<lower=-pi(),upper=pi()>[pp] theta_lon;
  vector<lower=-pi()/2.0,upper=pi()/2.0>[d-pp] theta_lat;
}
transformed parameters{
  vector<lower=-pi(),upper=pi()>[d] theta = set_theta(n,p,theta_lon,theta_lat);
  matrix[n,p] Y = givens_lp(n,p,theta);
  real abs_theta = abs(theta[1]);
}
model {
  abs_theta ~ normal(pi(), 0.2);
}
