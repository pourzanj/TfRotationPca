functions {
  
  matrix fill_W(matrix W_prot, vector Lambda_prot, matrix W_teg, vector Lambda_teg, matrix B_prot, vector Gamma_prot, matrix B_teg, vector Gamma_teg) {
    
    int n_prot = rows(W_prot);
    int n_teg = rows(W_teg);
    int p_common = cols(W_prot);
    int p_B_prot = cols(B_prot);
    int p_B_teg = cols(B_teg);
    
    matrix[n_prot+n_teg, p_common + p_B_prot + p_B_teg] W = rep_matrix(0, n_prot+n_teg, p_common + p_B_prot + p_B_teg);
   
    // Construct column-wise
    W[1:n_prot, 1:p_common] = diag_post_multiply(W_prot, Lambda_prot);
    W[(n_prot+1):(n_prot+n_teg), 1:p_common] = diag_post_multiply(W_teg, Lambda_teg);
    W[1:n_prot,(p_common+1):(p_common+p_B_prot)] = diag_post_multiply(B_prot, Gamma_prot);
    W[(n_prot+1):(n_prot+n_teg),(p_common+p_B_prot+1):(p_common+p_B_prot+p_B_teg)] = diag_post_multiply(B_teg, Gamma_teg);
    
    return(W);
  }
  
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
  int<lower=0> n_prot;
  int<lower=0> p_W_prot;
  int<lower=0> p_B_prot;

  int<lower=0> n_teg;
  int<lower=0> p_W_teg;
  int<lower=0> p_B_teg;

  int<lower=0> p_common;

  int N;//number of patients
  matrix[n_prot+n_teg,n_prot+n_teg] SigmaHat;
}
transformed data {
  int<lower=0> d_W_prot = n_prot*p_W_prot - p_W_prot*(p_W_prot+1)/2;
  int<lower=0> d_B_prot = n_prot*p_B_prot - p_B_prot*(p_B_prot+1)/2;
  int<lower=0> d_W_teg = n_teg*p_W_teg - p_W_teg*(p_W_teg+1)/2;
  int<lower=0> d_B_teg = n_teg*p_B_teg - p_B_teg*(p_B_teg+1)/2;
}
parameters {
  
  // Givens parameters for orthogonal matrices W and B
  vector[p_W_prot] x_lon_W_prot; vector[p_W_prot] y_lon_W_prot;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_W_prot-p_W_prot] theta_lat_W_prot;
  
  vector[p_B_prot] x_lon_B_prot; vector[p_B_prot] y_lon_B_prot;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_B_prot-p_B_prot] theta_lat_B_prot;
  
  vector[p_W_teg] x_lon_W_teg; vector[p_W_teg] y_lon_W_teg;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_W_teg-p_W_teg] theta_lat_W_teg;
  
  vector[p_B_teg] x_lon_B_teg; vector[p_B_teg] y_lon_B_teg;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_B_teg-p_B_teg] theta_lat_B_teg;
  
  // diagonal terms Lambda and Gamma
  // ordered vector in Stan goes from smallest to largest so we must reverse them
  positive_ordered[p_W_prot] Lambda_prot_rev;
  positive_ordered[p_B_prot] Gamma_prot_rev;
  positive_ordered[p_W_teg] Lambda_teg_rev;
  positive_ordered[p_B_teg] Gamma_teg_rev;
  
  // noise distribution
  real<lower=0> sigma;
}

transformed parameters{
  
  vector[p_W_prot] theta_lon_W_prot; vector[d_W_prot] theta_W_prot; vector[p_W_prot] r_W_prot; matrix[n_prot,p_W_prot] W_prot;
  vector[p_B_prot] theta_lon_B_prot; vector[d_B_prot] theta_B_prot; vector[p_B_prot] r_B_prot; matrix[n_prot,p_B_prot] B_prot;
  
  vector[p_W_teg] theta_lon_W_teg; vector[d_W_teg] theta_W_teg; vector[p_W_teg] r_W_teg; matrix[n_teg,p_W_teg] W_teg;
  vector[p_B_teg] theta_lon_B_teg; vector[d_B_teg] theta_B_teg; vector[p_B_teg] r_B_teg; matrix[n_teg,p_B_teg] B_teg;
  
  vector[p_W_prot] Lambda_prot;
  vector[p_B_prot] Gamma_prot;
  vector[p_W_teg] Lambda_teg;
  vector[p_B_teg] Gamma_teg;
  
  // get orthogonal matrices from Givens representation
  theta_lon_W_prot = atan2_vec(x_lon_W_prot, y_lon_W_prot); r_W_prot = hypot_vec(x_lon_W_prot, y_lon_W_prot);
  theta_W_prot = mirror(n_prot, p_W_prot, theta_lon_W_prot, theta_lat_W_prot); W_prot = givens_lp(n_prot,p_W_prot,theta_W_prot);
  
  theta_lon_B_prot = atan2_vec(x_lon_B_prot, y_lon_B_prot); r_B_prot = hypot_vec(x_lon_B_prot, y_lon_B_prot);
  theta_B_prot = mirror(n_prot, p_B_prot, theta_lon_B_prot, theta_lat_B_prot); B_prot = givens_lp(n_prot,p_B_prot,theta_B_prot);
  
  theta_lon_W_teg = atan2_vec(x_lon_W_teg, y_lon_W_teg); r_W_teg = hypot_vec(x_lon_W_teg, y_lon_W_teg);
  theta_W_teg = mirror(n_teg, p_W_teg, theta_lon_W_teg, theta_lat_W_teg); W_teg = givens_lp(n_teg,p_W_teg,theta_W_teg);
  
  theta_lon_B_teg = atan2_vec(x_lon_B_teg, y_lon_B_teg); r_B_teg = hypot_vec(x_lon_B_teg, y_lon_B_teg);
  theta_B_teg = mirror(n_teg, p_B_teg, theta_lon_B_teg, theta_lat_B_teg); B_teg = givens_lp(n_teg,p_B_teg,theta_B_teg);
  
  // reverse and square lambda
  for(i in 1:p_W_prot) Lambda_prot[i] = Lambda_prot_rev[p_W_prot - i + 1];
  for(i in 1:p_B_prot) Gamma_prot[i] = Gamma_prot_rev[p_B_prot - i + 1];
  for(i in 1:p_W_teg) Lambda_teg[i] = Lambda_teg_rev[p_W_teg - i + 1];
  for(i in 1:p_B_teg) Gamma_teg[i] = Gamma_teg_rev[p_B_teg - i + 1];
}

model {

  int n_total;
  int p_total;

  matrix[n_prot+n_teg, p_common + p_B_prot + p_B_teg] W = fill_W(W_prot, Lambda_prot, W_teg, Lambda_teg, B_prot, Gamma_prot, B_teg, Gamma_teg);
  matrix[n_prot+n_teg, n_prot+n_teg] C = W*W' + diag_matrix(rep_vector(sigma^2, n_prot+n_teg));
  
  // likelihood
  target += -(N/2)*log(determinant(C)) - (N/2)*trace(C\SigmaHat);
}
