functions {
  
    matrix left_rotation(matrix A, real angle, int n, int p, int i, int j) {
        matrix[n, p] RA;
        RA = A;
        
        RA[i,] = cos(angle)*A[i,] - sin(angle)*A[j,];
        RA[j,] = sin(angle)*A[i,] + cos(angle)*A[j,];
        
        return RA;
    }
    
    matrix right_rotation(matrix A, real angle, int n, int i, int j) {
        matrix[n, n] AR;
        AR = A;
        
        AR[,i] = cos(angle)*A[,i] + sin(angle)*A[,j];
        AR[,j] = -sin(angle)*A[,i] + cos(angle)*A[,j];
        
        return AR;
    }
    
    matrix d_rotation_matrix(real angle, int n, int i, int j) {
        matrix[n, n] dR = diag_matrix(rep_vector(0, n));
        
        dR[i, i] = -sin(angle);
        dR[i, j] = -cos(angle);
        dR[j, i] = cos(angle);
        dR[j, j] = -sin(angle);
        
        return dR;
    }
    
    matrix[] generate_forward_pgivens(vector angles, int n, int p) {
        int d = n*p - p*(p+1)/2;
        int idx;
        matrix[n, n] G;
        matrix[n, n] partial_givens[d+1];
        matrix[n, n] R;
        int pp;
        if(p == n) pp = p - 1;
        else pp = p;
        
        G = diag_matrix(rep_vector(1, n));
        idx = 1;
        partial_givens[1] = G;
        for(i in 1:pp){
            for(j in i+1:n){
                G = right_rotation(G, angles[idx], n, i, j);
                partial_givens[idx + 1] = G;
                idx = idx + 1;
            }
        }
        
        return partial_givens;      
    }
    
    matrix[] generate_reverse_pgivens(vector angles, int n, int p) {
        int d = n*p - p*(p+1)/2;
        int idx;
        matrix[n, n] G_eye;
        matrix[n, p] G;
        matrix[n, p] partial_givens[d+1];
        matrix[n, n] R;
        int pp;
        if(p == n) pp = p - 1;
        else pp = p;

        G_eye = diag_matrix(rep_vector(1, n));
        G = G_eye[,1:p];
        
        partial_givens[d+1] = G;
        idx = d;
        for(i in 1:pp){
            int i_st = pp - i + 1;
            for(j in i_st+1:n){
                G = left_rotation(G, angles[idx], n, p, i_st, n - j + i_st + 1);
                partial_givens[idx] = G;
                idx = idx - 1;
            }
        }
        
        return partial_givens;      
    }
    
    matrix[] generate_givens_jacobians(matrix[] partial_givens_forward, matrix[] partial_givens_reverse, vector angles, int n, int p) {
        int d = n*p - p*(p+1)/2;
        matrix[n,p] derivative_list[d];
        matrix[n,d] givens_jacobian[p];
        int idx = 1;
        int pp;
        if(p == n) pp = p - 1;
        else pp = p;
        
        for(i in 1:p){
            for(j in i+1:n){
                matrix[n,n] dR = d_rotation_matrix(angles[idx], n, i, j);
                matrix[n,n] a = partial_givens_forward[idx];
                matrix[n,p] b = partial_givens_reverse[idx + 1];
                
                derivative_list[idx] = a * dR * b;
                idx = idx + 1;
            }
        }
        
        for(i in 1:pp) {
            for(j in 1:d) {
                vector[n] t = derivative_list[j][,i];
                matrix[n,d] z = givens_jacobian[i];
                z[,j] = t;
                givens_jacobian[i] = z;
            }
        }
        
        return givens_jacobian;
        
    }
    
    matrix area_form_lp(vector angles, int n, int p) {
        int d = n*p - p*(p+1)/2;
        int idx;
        matrix[n, n] givens;
        matrix[n, n] partial_givens_forward[d+1];
        matrix[n, p] partial_givens_reverse[d+1];
        matrix[n, d] givens_jacobians[p];
        matrix[d, d] area_mat;
        int pp;
        if(p == n) pp = p - 1;
        else pp = p;
        

        partial_givens_forward = generate_forward_pgivens(angles, n, p);
        partial_givens_reverse = generate_reverse_pgivens(angles, n, p);
        givens = partial_givens_forward[d+1];
        
        givens_jacobians = generate_givens_jacobians(partial_givens_forward, partial_givens_reverse, angles, n, p);
        
        idx = 1;
        for(i in 1:pp){
            matrix[d, n-i] one_forms;
            one_forms = (givens'[i+1:n,] * givens_jacobians[i])';
            for(j in 1:n-i) {
                area_mat[,idx] = one_forms[,j]; 
                idx = idx + 1;
            }
        }
        
        target += log_determinant(area_mat); 
        return givens[,1:p];
    }
}

data {
  int n_0;
  int p_0;

  int n_1;
  int p_1;

  int n_2;
  int p_2;

  int n_3;
  int p_3;

  int n_4;
  int p_4;

  int p_common;
  
  int d_W_0;
  int d_B_0;

  int d_W_1;
  int d_B_1;

  int d_W_2;
  int d_B_2;

  int d_W_3;
  int d_B_3;

  int d_W_4;
  int d_B_4;
  
  real sigmaSqHyperPrior;
  
  real ardHyperHyperPrior;
  real<lower =0> sparseHyperHyperPrior;
  
  int J; //number groups
  int N[J]; //numer of patients for each group
  
  matrix[n_0+n_1+n_2+n_3+n_4,n_0+n_1+n_2+n_3+n_4] SigmaHat[J]; //estimated covariance matrice for each group
}

parameters {
  
  //TEG
  vector<lower = -pi()/2, upper = pi()/2>[d_W_0] theta_W_0[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_W_0] mu_theta_W_0;
  vector<lower = 0>[d_W_0] Sigma_theta_W_0;
  
  positive_ordered[p_common] lambda_0_reversed[J];
  positive_ordered[p_common] mu_lambda_0_reversed;
  vector<lower = 0>[p_common] Sigma_lambda_0_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_0] theta_B_0[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_B_0] mu_theta_B_0;
  vector<lower = 0>[d_B_0] Sigma_theta_B_0;
  
  positive_ordered[p_0] gamma_0_reversed[J];
  positive_ordered[p_0] mu_gamma_0_reverse;
  vector<lower = 0>[p_0] Sigma_gamma_0_reverse;
  
  //PT
  vector<lower = -pi()/2, upper = pi()/2>[d_W_1] theta_W_1[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_W_1] mu_theta_W_1;
  vector<lower = 0>[d_W_1] Sigma_theta_W_1;
  
  positive_ordered[p_common] lambda_1_reversed[J];
  positive_ordered[p_common] mu_lambda_1_reversed;
  vector<lower = 0>[p_common] Sigma_lambda_1_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_1] theta_B_1[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_B_1] mu_theta_B_1;
  vector<lower = 0>[d_B_1] Sigma_theta_B_1;
  
  positive_ordered[p_1] gamma_1_reversed[J];
  positive_ordered[p_1] mu_gamma_1_reverse;
  vector<lower = 0>[p_1] Sigma_gamma_1_reverse;

  //2
  vector<lower = -pi()/2, upper = pi()/2>[d_W_2] theta_W_2[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_W_2] mu_theta_W_2;
  vector<lower = 0>[d_W_2] Sigma_theta_W_2;
  
  positive_ordered[p_common] lambda_2_reversed[J];
  positive_ordered[p_common] mu_lambda_2_reversed;
  vector<lower = 0>[p_common] Sigma_lambda_2_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_2] theta_B_2[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_B_2] mu_theta_B_2;
  vector<lower = 0>[d_B_2] Sigma_theta_B_2;
  
  positive_ordered[p_2] gamma_2_reversed[J];
  positive_ordered[p_2] mu_gamma_2_reverse;
  vector<lower = 0>[p_2] Sigma_gamma_2_reverse;

  //3
  vector<lower = -pi()/2, upper = pi()/2>[d_W_3] theta_W_3[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_W_3] mu_theta_W_3;
  vector<lower = 0>[d_W_3] Sigma_theta_W_3;
  
  positive_ordered[p_common] lambda_3_reversed[J];
  positive_ordered[p_common] mu_lambda_3_reversed;
  vector<lower = 0>[p_common] Sigma_lambda_3_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_3] theta_B_3[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_B_3] mu_theta_B_3;
  vector<lower = 0>[d_B_3] Sigma_theta_B_3;
  
  positive_ordered[p_3] gamma_3_reversed[J];
  positive_ordered[p_3] mu_gamma_3_reverse;
  vector<lower = 0>[p_3] Sigma_gamma_3_reverse;

  //4
  vector<lower = -pi()/2, upper = pi()/2>[d_W_4] theta_W_4[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_W_4] mu_theta_W_4;
  vector<lower = 0>[d_W_4] Sigma_theta_W_4;
  
  positive_ordered[p_common] lambda_4_reversed[J];
  positive_ordered[p_common] mu_lambda_4_reversed;
  vector<lower = 0>[p_common] Sigma_lambda_4_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_4] theta_B_4[J];
  vector<lower = -pi()/2, upper = pi()/2>[d_B_4] mu_theta_B_4;
  vector<lower = 0>[d_B_4] Sigma_theta_B_4;
  
  //positive_ordered[p_4] gamma_4_reversed[J];
  //positive_ordered[p_4] mu_gamma_4_reverse;
  //vector<lower = 0>[p_4] Sigma_gamma_4_reverse;
  
  //PT doesn't have private latent variables
  
  //hyper-priors
  //real<lower = 0> sparseHyperPrior;
  real<lower = 0> sigmaSq;
  
  real<lower = 0> ardHyperPrior;
  real<lower = 0> sparseHyperPrior;
}

transformed parameters{
  //TEG
  vector<lower=0>[p_common] lambda_0[J];
  vector<lower=0>[p_0] gamma_0[J];

  vector<lower=0>[p_common] lambda_1[J];
  vector<lower=0>[p_1] gamma_1[J];

  vector<lower=0>[p_common] lambda_2[J];
  vector<lower=0>[p_2] gamma_2[J];

  vector<lower=0>[p_common] lambda_3[J];
  vector<lower=0>[p_3] gamma_3[J];

  vector<lower=0>[p_common] lambda_4[J];
  //vector<lower=0>[p_4] gamma_4[J];
  
  for(j in 1:J) {
    for (i in 1:p_common) lambda_0[j][i] = lambda_0_reversed[j][p_common - i + 1];
    for (i in 1:p_0) gamma_0[j][i] = gamma_0_reversed[j][p_0 - i + 1];
  
    for (i in 1:p_common) lambda_1[j][i] = lambda_1_reversed[j][p_common - i + 1];
    for (i in 1:p_1) gamma_1[j][i] = gamma_1_reversed[j][p_1 - i + 1];
  
    for (i in 1:p_common) lambda_2[j][i] = lambda_2_reversed[j][p_common - i + 1];
    for (i in 1:p_2) gamma_2[j][i] = gamma_2_reversed[j][p_2 - i + 1];
  
    for (i in 1:p_common) lambda_3[j][i] = lambda_3_reversed[j][p_common - i + 1];
    for (i in 1:p_3) gamma_3[j][i] = gamma_3_reversed[j][p_3 - i + 1];
  
    for (i in 1:p_common) lambda_4[j][i] = lambda_4_reversed[j][p_common - i + 1];
    //for (i in 1:p_4) gamma_4[j][i] = gamma_4_reversed[j][p_4 - i + 1];
  }
}

model {

  int n_total;
  int p_total;
  
  matrix[n_0, p_common] W_0[J];
  matrix[n_0, p_0] B_0[J];
  
  matrix[n_1, p_common] W_1[J];
  matrix[n_1, p_1] B_1[J];

  matrix[n_2, p_common] W_2[J];
  matrix[n_2, p_2] B_2[J];

  matrix[n_3, p_common] W_3[J];
  matrix[n_3, p_3] B_3[J];

  matrix[n_4, p_common] W_4[J];
  //matrix[n_4, p_4] B_4[J];
  
  matrix[n_0 + n_1 + n_2 + n_3 + n_4, p_common + p_0 + p_1 + p_2 + p_3 + p_4] W[J];
  matrix[n_0 + n_1 + n_2 + n_3 + n_4, n_0 + n_1 + n_2 + n_3 + n_4] Id_n;
  matrix[n_0 + n_1 + n_2 + n_3 + n_4, n_0 + n_1 + n_2 + n_3 + n_4] C[J];
  
  n_total = n_0 + n_1 + n_2 + n_3 + n_4;
  p_total = p_common + p_0 + p_1 + p_2 + p_3 + p_4;
  
  //Hyper Priors
  Id_n = diag_matrix(rep_vector(1, n_total));
  sigmaSq ~ normal(1, sigmaSqHyperPrior);
  
  //sparse Pca priors on W
  sparseHyperPrior ~ normal(0.1, sparseHyperHyperPrior);
  
  //control on lambda size for regularization
  ardHyperPrior ~ normal(0.1, ardHyperHyperPrior);
  
  //construct W and B matrices (CCA model) for each of the J sub-cohorts
  for(j in 1:J) {
    
    theta_W_0[j] ~ multi_normal(mu_theta_W_0, diag_matrix(Sigma_theta_W_0));
    W_0[j] = area_form_lp(theta_W_0[j], n_0, p_common);
    theta_B_0[j] ~ multi_normal(mu_theta_B_0, diag_matrix(Sigma_theta_B_0));
    B_0[j] = area_form_lp(theta_B_0[j], n_0, p_0);
  
    theta_W_1[j] ~ multi_normal(mu_theta_W_1, diag_matrix(Sigma_theta_W_1));
    W_1[j] = area_form_lp(theta_W_1[j], n_1, p_common);
    theta_B_1[j] ~ multi_normal(mu_theta_B_1, diag_matrix(Sigma_theta_B_1));
    B_1[j] = area_form_lp(theta_B_1[j], n_1, p_1);
  
    theta_W_2[j] ~ multi_normal(mu_theta_W_2, diag_matrix(Sigma_theta_W_2));
    W_2[j] = area_form_lp(theta_W_2[j], n_2, p_common);
    theta_B_2[j] ~ multi_normal(mu_theta_B_2, diag_matrix(Sigma_theta_B_2));
    B_2[j] = area_form_lp(theta_B_2[j], n_2, p_2);
  
    theta_W_3[j] ~ multi_normal(mu_theta_W_3, diag_matrix(Sigma_theta_W_3));
    W_3[j] = area_form_lp(theta_W_3[j], n_3, p_common);
    theta_B_3[j] ~ multi_normal(mu_theta_B_3, diag_matrix(Sigma_theta_B_3));
    B_3[j] = area_form_lp(theta_B_3[j], n_3, p_3);
  
    theta_W_4[j] ~ multi_normal(mu_theta_W_4, diag_matrix(Sigma_theta_W_4));
    W_4[j] = area_form_lp(theta_W_4[j], n_4, p_common);
    //theta_B_4[j] ~ multi_normal(mu_theta_B_4, diag_matrix(Sigma_theta_B_4));
    //B_4[j] = area_form_lp(theta_B_4[j], n_4, p_4);
    
    //create empty grand W for sub-cohort
    W[j] = rep_matrix(0, n_total, p_total);
    
    //fill in the W for the sub-cohort row wise
    W[j,1:n_0, 1:p_common] = W_0[j] * diag_matrix(lambda_0[j]);
    W[j,n_0 + 1:n_0 + n_1, 1:p_common] = W_1[j] * diag_matrix(lambda_1[j]);
    W[j,n_0 + n_1 + 1:n_0 + n_1 + n_2, 1:p_common] = W_2[j] * diag_matrix(lambda_2[j]);
    W[j,n_0 + n_1 + n_2 + 1:n_0 + n_1 + n_2 + n_3, 1:p_common] = W_3[j] * diag_matrix(lambda_3[j]);
    W[j,n_0 + n_1 + n_2 + n_3 + 1:n_total, 1:p_common] = W_4[j] * diag_matrix(lambda_4[j]);
  
    W[j,1:n_0, p_common + 1:p_common + p_0] = B_0[j] * diag_matrix(gamma_0[j]);
    W[j,n_0 + 1:n_0 + n_1, p_common + p_0 + 1:p_common + p_0 + p_1] = B_1[j] * diag_matrix(gamma_1[j]);
    W[j,n_0 + n_1 + 1:n_0 + n_1 + n_2, p_common + p_0 + p_1 + 1:p_common + p_0 + p_1 + p_2] = B_2[j] * diag_matrix(gamma_2[j]);
    W[j,n_0 + n_1 + n_2 + 1:n_0 + n_1 + n_2 + n_3, p_common + p_0 + p_1 + p_2 + 1:p_common + p_0 + p_1 + p_2 + p_3] = B_3[j] * diag_matrix(gamma_3[j]);
    //W[j,n_0 + n_1 + n_2 + n_3 + 1:n_total, p_common + p_0 + p_1 + p_2 + p_3 + 1:p_total] = B_4[j] * diag_matrix(gamma_4[j]);
    
    // sparsity prior on entries of W and B (sparse PCA)
    for(i in 1:p_common) {
      W_0[j][,i] ~ cauchy(0, sparseHyperPrior);
      W_1[j][,i] ~ cauchy(0, sparseHyperPrior);
      W_2[j][,i] ~ cauchy(0, sparseHyperPrior);
      W_3[j][,i] ~ cauchy(0, sparseHyperPrior);
      W_4[j][,i] ~ cauchy(0, sparseHyperPrior);
    }
    
    for(i in 1:p_0) B_0[j][,i] ~ cauchy(0, sparseHyperPrior);
    for(i in 1:p_1) B_1[j][,i] ~ cauchy(0, sparseHyperPrior);
    for(i in 1:p_2) B_2[j][,i] ~ cauchy(0, sparseHyperPrior);
    for(i in 1:p_3) B_3[j][,i] ~ cauchy(0, sparseHyperPrior);
    //for(i in 1:p_4) B_4[j][,i] ~ cauchy(0, sparseHyperPrior);

    //control on lambda size for regularization
    lambda_0_reversed[j] ~ cauchy(0, ardHyperPrior);
    lambda_0_reversed[j]~ multi_normal(mu_lambda_0_reversed, diag_matrix(Sigma_lambda_0_reversed));
    
    lambda_1_reversed[j] ~ cauchy(0, ardHyperPrior);
    lambda_1_reversed[j] ~ multi_normal(mu_lambda_1_reversed, diag_matrix(Sigma_lambda_1_reversed));
    
    lambda_2_reversed[j] ~ cauchy(0, ardHyperPrior);
    lambda_2_reversed[j] ~ multi_normal(mu_lambda_2_reversed, diag_matrix(Sigma_lambda_2_reversed));
    
    lambda_3_reversed[j] ~ cauchy(0, ardHyperPrior);
    lambda_3_reversed[j] ~ multi_normal(mu_lambda_3_reversed, diag_matrix(Sigma_lambda_3_reversed));
    
    lambda_4_reversed[j] ~ cauchy(0, ardHyperPrior);
    lambda_4_reversed[j] ~ multi_normal(mu_lambda_4_reversed, diag_matrix(Sigma_lambda_4_reversed));
    
    gamma_0_reversed[j] ~ cauchy(0, ardHyperPrior);
    gamma_0_reversed[j] ~ multi_normal(mu_gamma_0_reverse, diag_matrix(Sigma_gamma_0_reverse));
    
    gamma_1_reversed[j] ~ cauchy(0, ardHyperPrior);
    gamma_1_reversed[j] ~ multi_normal(mu_gamma_1_reverse, diag_matrix(Sigma_gamma_1_reverse));
    
    gamma_2_reversed[j] ~ cauchy(0, ardHyperPrior);
    gamma_2_reversed[j] ~ multi_normal(mu_gamma_2_reverse, diag_matrix(Sigma_gamma_2_reverse));
    
    gamma_3_reversed[j] ~ cauchy(0, ardHyperPrior);
    gamma_3_reversed[j] ~ multi_normal(mu_gamma_3_reverse, diag_matrix(Sigma_gamma_3_reverse));
    
    //gamma_4_reversed[j] ~ cauchy(0, ardHyperPrior);
    //gamma_4_reversed[j] ~ multi_normal(mu_gamma_4_reverse, diag_matrix(Sigma_gamma_4_reverse));
    
    //gamma_4_reversed ~ cauchy(0, ardHyperPrior);
    
    //PPCA likelihood from Ch. 12 of Kevin Murphy
    C[j] = W[j]*W[j]' + sigmaSq*Id_n;
  
    target += -(N[j]/2)*log(determinant(C[j])) -(N[j]/2)*trace(C[j]\SigmaHat[j]);
  }
  
}
