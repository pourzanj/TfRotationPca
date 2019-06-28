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
  int<lower=0> J;//number of groups
  
  int<lower=0> n_prot;
  int<lower=0> p_W_prot;
  int<lower=0> p_B_prot;

  int<lower=0> n_teg;
  int<lower=0> p_W_teg;
  int<lower=0> p_B_teg;

  int<lower=0> p_common;

  int N[J];//number of patients
  matrix[n_prot+n_teg,n_prot+n_teg] SigmaHat[J];
}
transformed data {
  int<lower=0> d_W_prot = n_prot*p_W_prot - p_W_prot*(p_W_prot+1)/2;
  int<lower=0> d_B_prot = n_prot*p_B_prot - p_B_prot*(p_B_prot+1)/2;
  int<lower=0> d_W_teg = n_teg*p_W_teg - p_W_teg*(p_W_teg+1)/2;
  int<lower=0> d_B_teg = n_teg*p_B_teg - p_B_teg*(p_B_teg+1)/2;
}
parameters {

  // hyperparameters on Givens representations
  // location parameter mu should be implemented like the angle coordinates themselves
  vector[p_W_prot] x_lon_mu_W_prot; vector[p_W_prot] y_lon_mu_W_prot;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_W_prot-p_W_prot] theta_lat_mu_W_prot;
  
  vector[p_B_prot] x_lon_mu_B_prot; vector[p_B_prot] y_lon_mu_B_prot;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_B_prot-p_B_prot] theta_lat_mu_B_prot;
  
  vector[p_W_teg] x_lon_mu_W_teg; vector[p_W_teg] y_lon_mu_W_teg;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_W_teg-p_W_teg] theta_lat_mu_W_teg;
  
  vector[p_B_teg] x_lon_mu_B_teg; vector[p_B_teg] y_lon_mu_B_teg;
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_B_teg-p_B_teg] theta_lat_mu_B_teg;
  
  // variance hyperparameters for Givens representations
  vector<lower=0>[d_W_prot] sigma_theta_W_prot;
  vector<lower=0>[d_B_prot] sigma_theta_B_prot;
  vector<lower=0>[d_W_teg] sigma_theta_W_teg;
  vector<lower=0>[d_B_teg] sigma_theta_B_teg;
  
  // hyperparameters on Lambda and Gamma
  vector[p_W_prot] mu_Lambda_prot; vector<lower=0>[p_W_prot] sigma_Lambda_prot;
  vector[p_B_prot] mu_Gamma_prot; vector<lower=0>[p_B_prot] sigma_Gamma_prot;
  vector[p_W_teg] mu_Lambda_teg; vector<lower=0>[p_W_teg] sigma_Lambda_teg;
  vector[p_B_teg] mu_Gamma_teg; vector<lower=0>[p_B_teg] sigma_Gamma_teg;
  
  // Givens parameters for orthogonal matrices W and B
  vector[p_W_prot] x_lon_W_prot[J]; vector[p_W_prot] y_lon_W_prot[J];
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_W_prot-p_W_prot] theta_lat_W_prot[J];
  
  vector[p_B_prot] x_lon_B_prot[J]; vector[p_B_prot] y_lon_B_prot[J];
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_B_prot-p_B_prot] theta_lat_B_prot[J];
  
  vector[p_W_teg] x_lon_W_teg[J]; vector[p_W_teg] y_lon_W_teg[J];
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_W_teg-p_W_teg] theta_lat_W_teg[J];
  
  vector[p_B_teg] x_lon_B_teg[J]; vector[p_B_teg] y_lon_B_teg[J];
  vector<lower=-pi()/2.0 + 1e-5,upper=pi()/2.0 - 1e-5>[d_B_teg-p_B_teg] theta_lat_B_teg[J];
  
  // diagonal terms Lambda and Gamma
  // ordered vector in Stan goes from smallest to largest so we must reverse them
  positive_ordered[p_W_prot] Lambda_prot_rev[J];
  positive_ordered[p_B_prot] Gamma_prot_rev[J];
  positive_ordered[p_W_teg] Lambda_teg_rev[J];
  positive_ordered[p_B_teg] Gamma_teg_rev[J];
  
  // noise distribution
  real<lower=0> sigma;
}

transformed parameters{
  
  // declare hierarchical orthogonal matrices
  vector[p_W_prot] theta_lon_mu_W_prot; vector[d_W_prot] theta_mu_W_prot; vector[p_W_prot] r_mu_W_prot; matrix[n_prot,p_W_prot] mu_W_prot;
  vector[p_B_prot] theta_lon_mu_B_prot; vector[d_B_prot] theta_mu_B_prot; vector[p_B_prot] r_mu_B_prot; matrix[n_prot,p_B_prot] mu_B_prot;
  
  vector[p_W_teg] theta_lon_mu_W_teg; vector[d_W_teg] theta_mu_W_teg; vector[p_W_teg] r_mu_W_teg; matrix[n_teg,p_W_teg] mu_W_teg;
  vector[p_B_teg] theta_lon_mu_B_teg; vector[d_B_teg] theta_mu_B_teg; vector[p_B_teg] r_mu_B_teg; matrix[n_teg,p_B_teg] mu_B_teg;
  
  // declare orthogonal matrices for each group
  vector[p_W_prot] theta_lon_W_prot[J]; vector[d_W_prot] theta_W_prot[J]; vector[p_W_prot] r_W_prot[J]; matrix[n_prot,p_W_prot] W_prot[J];
  vector[p_B_prot] theta_lon_B_prot[J]; vector[d_B_prot] theta_B_prot[J]; vector[p_B_prot] r_B_prot[J]; matrix[n_prot,p_B_prot] B_prot[J];
  
  vector[p_W_teg] theta_lon_W_teg[J]; vector[d_W_teg] theta_W_teg[J]; vector[p_W_teg] r_W_teg[J]; matrix[n_teg,p_W_teg] W_teg[J];
  vector[p_B_teg] theta_lon_B_teg[J]; vector[d_B_teg] theta_B_teg[J]; vector[p_B_teg] r_B_teg[J]; matrix[n_teg,p_B_teg] B_teg[J];
  
  // declare Lambda and Gamma for each group  
  vector[p_W_prot] Lambda_prot[J];
  vector[p_B_prot] Gamma_prot[J];
  vector[p_W_teg] Lambda_teg[J];
  vector[p_B_teg] Gamma_teg[J];
  
  // set hierarchical orthogonal matrices
  theta_lon_mu_W_prot = atan2_vec(x_lon_mu_W_prot, y_lon_mu_W_prot); r_mu_W_prot = hypot_vec(x_lon_mu_W_prot, y_lon_mu_W_prot);
  theta_mu_W_prot = mirror(n_prot, p_W_prot, theta_lon_mu_W_prot, theta_lat_mu_W_prot); mu_W_prot = givens_lp(n_prot,p_W_prot,theta_mu_W_prot);
    
  theta_lon_mu_B_prot = atan2_vec(x_lon_mu_B_prot, y_lon_mu_B_prot); r_mu_B_prot = hypot_vec(x_lon_mu_B_prot, y_lon_mu_B_prot);
  theta_mu_B_prot = mirror(n_prot, p_B_prot, theta_lon_mu_B_prot, theta_lat_mu_B_prot); mu_B_prot = givens_lp(n_prot,p_B_prot,theta_mu_B_prot);
    
  theta_lon_mu_W_teg = atan2_vec(x_lon_mu_W_teg, y_lon_mu_W_teg); r_mu_W_teg = hypot_vec(x_lon_mu_W_teg, y_lon_mu_W_teg);
  theta_mu_W_teg = mirror(n_teg, p_W_teg, theta_lon_mu_W_teg, theta_lat_mu_W_teg); mu_W_teg = givens_lp(n_teg,p_W_teg,theta_mu_W_teg);
    
  theta_lon_mu_B_teg = atan2_vec(x_lon_mu_B_teg, y_lon_mu_B_teg); r_mu_B_teg = hypot_vec(x_lon_mu_B_teg, y_lon_mu_B_teg);
  theta_mu_B_teg = mirror(n_teg, p_B_teg, theta_lon_mu_B_teg, theta_lat_mu_B_teg); mu_B_teg = givens_lp(n_teg,p_B_teg,theta_mu_B_teg);
  
  // set orthogonal matrices for each group
  for(j in 1:J) {
    // get orthogonal matrices from Givens representation
    theta_lon_W_prot[j] = atan2_vec(x_lon_W_prot[j], y_lon_W_prot[j]); r_W_prot[j] = hypot_vec(x_lon_W_prot[j], y_lon_W_prot[j]);
    theta_W_prot[j] = mirror(n_prot, p_W_prot, theta_lon_W_prot[j], theta_lat_W_prot[j]); W_prot[j] = givens_lp(n_prot,p_W_prot,theta_W_prot[j]);
    
    theta_lon_B_prot[j] = atan2_vec(x_lon_B_prot[j], y_lon_B_prot[j]); r_B_prot[j] = hypot_vec(x_lon_B_prot[j], y_lon_B_prot[j]);
    theta_B_prot[j] = mirror(n_prot, p_B_prot, theta_lon_B_prot[j], theta_lat_B_prot[j]); B_prot[j] = givens_lp(n_prot,p_B_prot,theta_B_prot[j]);
    
    theta_lon_W_teg[j] = atan2_vec(x_lon_W_teg[j], y_lon_W_teg[j]); r_W_teg[j] = hypot_vec(x_lon_W_teg[j], y_lon_W_teg[j]);
    theta_W_teg[j] = mirror(n_teg, p_W_teg, theta_lon_W_teg[j], theta_lat_W_teg[j]); W_teg[j] = givens_lp(n_teg,p_W_teg,theta_W_teg[j]);
    
    theta_lon_B_teg[j] = atan2_vec(x_lon_B_teg[j], y_lon_B_teg[j]); r_B_teg[j] = hypot_vec(x_lon_B_teg[j], y_lon_B_teg[j]);
    theta_B_teg[j] = mirror(n_teg, p_B_teg, theta_lon_B_teg[j], theta_lat_B_teg[j]); B_teg[j] = givens_lp(n_teg,p_B_teg,theta_B_teg[j]);
    
    // reverse and square lambda
    for(i in 1:p_W_prot) Lambda_prot[j][i] = Lambda_prot_rev[j][p_W_prot - i + 1];
    for(i in 1:p_B_prot) Gamma_prot[j][i] = Gamma_prot_rev[j][p_B_prot - i + 1];
    for(i in 1:p_W_teg) Lambda_teg[j][i] = Lambda_teg_rev[j][p_W_teg - i + 1];
    for(i in 1:p_B_teg) Gamma_teg[j][i] = Gamma_teg_rev[j][p_B_teg - i + 1];
  }
}

model {

  int n_total;
  int p_total;

  for(j in 1:J) {
    matrix[n_prot+n_teg, p_common + p_B_prot + p_B_teg] W = fill_W(W_prot[j], Lambda_prot[j], W_teg[j], Lambda_teg[j], B_prot[j], Gamma_prot[j], B_teg[j], Gamma_teg[j]);
    matrix[n_prot+n_teg, n_prot+n_teg] C = W*W' + diag_matrix(rep_vector(sigma^2, n_prot+n_teg));
    
    //prior
    for(i in 1:d_W_prot) theta_W_prot[j][i] ~ normal(theta_mu_W_prot[i], sigma_theta_W_prot[i]) T[-pi()/2, pi()/2];
    for(i in 1:d_B_prot) theta_B_prot[j][i] ~ normal(theta_mu_B_prot[i], sigma_theta_B_prot[i]) T[-pi()/2, pi()/2];
    for(i in 1:d_W_teg) theta_W_teg[j][i] ~ normal(theta_mu_W_teg[i], sigma_theta_W_teg[i]) T[-pi()/2, pi()/2];
    for(i in 1:d_B_teg) theta_B_teg[j][i] ~ normal(theta_mu_B_teg[i], sigma_theta_B_teg[i]) T[-pi()/2, pi()/2];
    
    for(i in 1:p_W_prot) Lambda_prot[j][i] ~ normal(mu_Lambda_prot[i], sigma_Lambda_prot[i]) T[0.0,];
    for(i in 1:p_B_prot) Gamma_prot[j][i] ~ normal(mu_Gamma_prot[i], sigma_Gamma_prot[i]) T[0.0,];
    for(i in 1:p_W_teg) Lambda_teg[j][i] ~ normal(mu_Lambda_teg[i], sigma_Lambda_teg[i]) T[0.0,];
    for(i in 1:p_B_teg) Gamma_teg[j][i] ~ normal(mu_Gamma_teg[i], sigma_Gamma_teg[i]) T[0.0,];
    
    // likelihood
    target += -(N[j]/2)*log(determinant(C)) - (N[j]/2)*trace(C\SigmaHat[j]);
  }
 
}
