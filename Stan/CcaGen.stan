functions {
  matrix rotation_matrix(real angle, int n, int i, int j) {
    matrix[n, n] R = diag_matrix(rep_vector(1, n));
    
    R[i, i] = cos(angle);
    R[i, j] = -sin(angle);
    R[j, i] = sin(angle);
    R[j, j] = cos(angle);
    
    return R;
  }
  
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
    if(p == n) pp = p-1;
    else pp = p;
    
    G = diag_matrix(rep_vector(1, n));
    idx = 1;
    partial_givens[1] = G;
    for(i in 1:pp){
      for(j in i+1:n){
        //R = rotation_matrix(angles[idx], n, i, j);
        //G = G * R;
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
    if(p == n) pp = p-1;
    else pp = p;
    
    G_eye = diag_matrix(rep_vector(1, n));
    G = G_eye[,1:p];
    
    partial_givens[d+1] = G;
    idx = d;
    for(i in 1:pp){
      int i_st = pp - i + 1;
      for(j in i_st+1:n){
        //R = rotation_matrix(angles[idx], n, i_st, n - j + i_st + 1);
        //G = R * G;
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
    if(p == n) pp = p-1;
    else pp = p;
    
    for(i in 1:pp){
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
  
  real area_form(matrix[] partial_givens_forward, matrix[] partial_givens_reverse, vector angles, int n, int p) {

    int d = n*p - p*(p+1)/2;
    int idx;
    matrix[n, n] givens;
    //matrix[n, n] partial_givens_forward[d+1];
    //matrix[n, p] partial_givens_reverse[d+1];
    matrix[n, d] givens_jacobians[p];
    matrix[d, d] area_mat;
    int pp;
    if(p == n) pp = p-1;
    else pp = p;
    /**
      * Create Partial Givens
    */
      
    //partial_givens_forward = generate_forward_pgivens(angles, n, p);
    //partial_givens_reverse = generate_reverse_pgivens(angles, n, p);
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
      
    return log(determinant(area_mat));
  }
}

data {
  int n_teg;
  int p_teg;
  int n_pt;
  int p_pt;
  int p_common;
  
  int d_W_teg;
  int d_B_teg;
  int d_W_pt;
  int d_B_pt;
  
  real sigmaSqHyperPrior;
  //real<lower =0> sparseHyperHyperPrior;
  int N; //num of observations
  matrix[n_teg+n_pt,n_teg+n_pt] SigmaHat;
}
parameters {
  //TEG
  vector<lower = -pi()/2, upper = pi()/2>[d_W_teg] theta_W_teg;
  positive_ordered[p_common] lambda_teg_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_teg] theta_B_teg;
  positive_ordered[p_teg] gamma_teg_reversed;
  
  //PT
  vector<lower = -pi()/2, upper = pi()/2>[d_W_pt] theta_W_pt;
  positive_ordered[p_common] lambda_pt_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_pt] theta_B_pt;
  positive_ordered[p_pt] gamma_pt_reversed;
  
  //hyper-priors
  //real<lower = 0> sparseHyperPrior;
  real<lower = 0> sigmaSq;
}
transformed parameters{
  //TEG
  vector<lower=0>[p_common] lambda_teg;
  vector<lower=0>[p_teg] gamma_teg;
  vector<lower=0>[p_common] lambda_pt;
  vector<lower=0>[p_pt] gamma_pt;
  
  for (i in 1:p_common) lambda_teg[i] = lambda_teg_reversed[p_common - i + 1];
  for (i in 1:p_teg) gamma_teg[i] = gamma_teg_reversed[p_teg - i + 1];
  for (i in 1:p_common) lambda_pt[i] = lambda_pt_reversed[p_common - i + 1];
  for (i in 1:p_pt) gamma_pt[i] = gamma_pt_reversed[p_pt - i + 1];
}
model {
  matrix[n_teg, n_teg] W_teg_partial_givens_forward[d_W_teg+1];
  matrix[n_teg, p_common] W_teg_partial_givens_reverse[d_W_teg+1];
  matrix[n_teg, n_teg] B_teg_partial_givens_forward[d_B_teg+1];
  matrix[n_teg, p_teg] B_teg_partial_givens_reverse[d_B_teg+1];
  
  matrix[n_pt, n_pt] W_pt_partial_givens_forward[d_W_pt+1];
  matrix[n_pt, p_common] W_pt_partial_givens_reverse[d_W_pt+1];
  matrix[n_pt, n_pt] B_pt_partial_givens_forward[d_B_pt+1];
  matrix[n_pt, p_pt] B_pt_partial_givens_reverse[d_B_pt+1];
  
  matrix[n_teg, n_teg] G_W_teg;
  matrix[n_teg, p_common] W_teg;
  matrix[n_teg, n_teg] G_B_teg;
  matrix[n_teg, p_teg] B_teg;
  
  matrix[n_pt, n_pt] G_W_pt;
  matrix[n_pt, p_common] W_pt;
  matrix[n_pt, n_pt] G_B_pt;
  matrix[n_pt, p_pt] B_pt;
  
  matrix[n_teg + n_pt, p_common + p_teg + p_pt] W;
  matrix[n_teg + n_pt, n_teg + n_pt] Id_n;
  matrix[n_teg + n_pt, n_teg + n_pt] C;
  
  //TEG Area-Forms
  W_teg_partial_givens_forward = generate_forward_pgivens(theta_W_teg, n_teg, p_common);
  W_teg_partial_givens_reverse = generate_reverse_pgivens(theta_W_teg, n_teg, p_common);
  target += area_form(W_teg_partial_givens_forward, W_teg_partial_givens_reverse, theta_W_teg, n_teg, p_common);
  
  B_teg_partial_givens_forward = generate_forward_pgivens(theta_B_teg, n_teg, p_teg);
  B_teg_partial_givens_reverse = generate_reverse_pgivens(theta_B_teg, n_teg, p_teg);
  target += area_form(B_teg_partial_givens_forward, B_teg_partial_givens_reverse, theta_B_teg, n_teg, p_teg);
  
  G_W_teg = W_teg_partial_givens_forward[d_W_teg+1];
  W_teg = G_W_teg[,1:p_common];
  
  G_B_teg = B_teg_partial_givens_forward[d_B_teg+1];
  B_teg = G_B_teg[,1:p_teg];
  
  //PT Area-Forms
  W_pt_partial_givens_forward = generate_forward_pgivens(theta_W_pt, n_pt, p_common);
  W_pt_partial_givens_reverse = generate_reverse_pgivens(theta_W_pt, n_pt, p_common);
  target += area_form(W_pt_partial_givens_forward, W_pt_partial_givens_reverse, theta_W_pt, n_pt, p_common);
  
  B_pt_partial_givens_forward = generate_forward_pgivens(theta_B_pt, n_pt, p_pt);
  B_pt_partial_givens_reverse = generate_reverse_pgivens(theta_B_pt, n_pt, p_pt);
  target += area_form(B_pt_partial_givens_forward, B_pt_partial_givens_reverse, theta_B_pt, n_pt, p_pt);
  
  G_W_pt = W_pt_partial_givens_forward[d_W_pt+1];
  W_pt = G_W_pt[,1:p_common];
  
  G_B_pt = B_pt_partial_givens_forward[d_B_pt+1];
  B_pt = G_B_pt[,1:p_pt];
  
  //create W matrix
  W = rep_matrix(0, n_teg + n_pt, p_common + p_teg + p_pt);
  W[1:n_teg, 1:p_common] = W_teg * diag_matrix(lambda_teg);
  W[(n_teg+1):(n_teg+n_pt),1:p_common] = W_pt * diag_matrix(lambda_pt);
  W[1:n_teg, (p_common+1):(p_common+p_teg)] = B_teg * diag_matrix(gamma_teg);
  W[(n_teg+1):(n_teg+n_pt),(p_common+p_teg+1):(p_common+p_teg+p_pt)] = B_pt * diag_matrix(gamma_pt);
  
  //priors
  sigmaSq ~ normal(1, sigmaSqHyperPrior);
  
  //sparse Pca priors on W
  //sparseHyperPrior ~ normal(0.01, sparseHyperHyperPrior);
  //for(i in 1:p) W[,i] ~ cauchy(0, sparseHyperPrior);
  
  //PPCA likelihood from Ch. 12 of Kevin Murphy
  Id_n = diag_matrix(rep_vector(1, n_teg + n_pt));
  C = W*W' + sigmaSq*Id_n;
  
  target += -(N/2)*log(determinant(C)) -(N/2)*trace(C\SigmaHat);
}