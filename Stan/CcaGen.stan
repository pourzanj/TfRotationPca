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
  int N; //num of observations
  matrix[n_0+n_1+n_2+n_3+n_4,n_0+n_1+n_2+n_3+n_4] SigmaHat;
}

parameters {
  //TEG
  vector<lower = -pi()/2, upper = pi()/2>[d_W_0] theta_W_0;
  positive_ordered[p_common] lambda_0_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_0] theta_B_0;
  positive_ordered[p_0] gamma_0_reversed;
  
  //PT
  vector<lower = -pi()/2, upper = pi()/2>[d_W_1] theta_W_1;
  positive_ordered[p_common] lambda_1_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_1] theta_B_1;
  positive_ordered[p_1] gamma_1_reversed;

  //2
  vector<lower = -pi()/2, upper = pi()/2>[d_W_2] theta_W_2;
  positive_ordered[p_common] lambda_2_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_2] theta_B_2;
  positive_ordered[p_1] gamma_2_reversed;

  //3
  vector<lower = -pi()/2, upper = pi()/2>[d_W_3] theta_W_3;
  positive_ordered[p_common] lambda_3_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_3] theta_B_3;
  positive_ordered[p_1] gamma_3_reversed;

  //4
  vector<lower = -pi()/2, upper = pi()/2>[d_W_4] theta_W_4;
  positive_ordered[p_common] lambda_4_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[d_B_4] theta_B_4;
  positive_ordered[p_1] gamma_4_reversed;
  
  //hyper-priors
  //real<lower = 0> sparseHyperPrior;
  real<lower = 0> sigmaSq;
}

transformed parameters{
  //TEG
  vector<lower=0>[p_common] lambda_0;
  vector<lower=0>[p_0] gamma_0;

  vector<lower=0>[p_common] lambda_1;
  vector<lower=0>[p_1] gamma_1;

  vector<lower=0>[p_common] lambda_2;
  vector<lower=0>[p_2] gamma_2;

  vector<lower=0>[p_common] lambda_3;
  vector<lower=0>[p_3] gamma_3;

  vector<lower=0>[p_common] lambda_4;
  vector<lower=0>[p_4] gamma_4;
  
  for (i in 1:p_common) lambda_0[i] = lambda_0_reversed[p_common - i + 1];
  for (i in 1:p_0) gamma_0[i] = gamma_0_reversed[p_0 - i + 1];

  for (i in 1:p_common) lambda_1[i] = lambda_1_reversed[p_common - i + 1];
  for (i in 1:p_1) gamma_1[i] = gamma_1_reversed[p_1 - i + 1];

  for (i in 1:p_common) lambda_2[i] = lambda_2_reversed[p_common - i + 1];
  for (i in 1:p_2) gamma_2[i] = gamma_p2_reversed[p_2 - i + 1];

  for (i in 1:p_common) lambda_3[i] = lambda_3_reversed[p_common - i + 1];
  for (i in 1:p_3) gamma_3[i] = gamma_3_reversed[p_3 - i + 1];

  for (i in 1:p_common) lambda_4[i] = lambda_4_reversed[p_common - i + 1];
  for (i in 1:p_4) gamma_4[i] = gamma_4_reversed[p_4 - i + 1];
}

model {

  int n_total;
  n_total = n_0 + n_1 + n_2 + n_3 + n_4;
  p_total = p_common + p_0 + p_1 + p_2 + p_3 + p_4;
  
  matrix[n_0, p_common] W_0;
  matrix[n_0, p_0] B_0;
  
  matrix[n_1, p_common] W_1;
  matrix[n_1, p_1] B_1;

  matrix[n_2, p_common] W_2;
  matrix[n_2, p_2] B_2;

  matrix[n_3, p_common] W_3;
  matrix[n_3, p_3] B_3;

  matrix[n_4, p_common] W_4;
  matrix[n_4, p_4] B_4;
  
  matrix[n_total, p_total] W;
  matrix[n_total, n_total] Id_n;
  matrix[n_total, n_total] C;

  W_0 = area_form_lp(theta_W_0, n_0, p_common);
  B_0 = area_form_lp(theta_B_0, n_0, p_0);

  W_1 = area_form_lp(theta_W_1, n_1, p_common);
  B_1 = area_form_lp(theta_W_1, n_1, p_1);

  W_2 = area_form_lp(theta_W_2, n_2, p_common);
  B_2 = area_form_lp(theta_W_2, n_2, p_2);

  W_3 = area_form_lp(theta_W_3, n_3, p_common);
  B_3 = area_form_lp(theta_W_3, n_3, p_3);

  W_4 = area_form_lp(theta_W_4, n_4, p_common);
  B_4 = area_form_lp(theta_W_4, n_4, p_4);
  
  W = rep_matrix(0, n_total, p_total);

  // Construct Row Wise
  W[1:n_0, 1:p_common] = W_0 * diag_matrix(lambda_0);
  W[n_0 + 1:n_0 + n_1, 1:p_common] = W_1 * diag_matrix(lambda_1);
  W[n_0 + n_1 + 1:n_0 + n_1 + n_2, 1:p_common] = W_2 * diag_matrix(lambda_2);
  W[n_0 + n_1 + n_2 + 1:n_0 + n_1 + n_2 + n_3, 1:p_common] = W_3 * diag_matrix(lambda_3);
  W[n_0 + n_1 + n_2 + n_3 + 1:n_total, 1:p_common] = W_4 * diag_matrix(lambda_4);

  W[1:n_0, p_common + 1:p_common + p_0] = B_0 * diag_matrix(gamma_0);
  W[n_0 + 1:n_0 + n_1, p_common + p_0 + 1:p_common + p_0 + p_1] = B_1 * diag_matrix(gamma_1);
  W[n_0 + n_1 + 1:n_0 + n_1 + n_2, p_common + p_0 + p_1 + 1:p_common + p_0 + p_1 + p_2] = B_2 * diag_matrix(gamma_2);
  W[n_0 + n_1 + n_2 + 1:n_0 + n_1 + n_2 + n_3, p_common + p_0 + p_1 + p_2 + 1:p_common + p_0 + p_1 + p_2 + p_3] = B_3 * diag_matrix(gamma_3);
  W[n_0 + n_1 + n_2 + n_3 + 1:n_total, p_common + p_0 + p_1 + p_2 + p_3 + 1:p_total] = B_4 * diag_matrix(gamma_4);
  
  //priors
  sigmaSq ~ normal(1, sigmaSqHyperPrior);
  
  //sparse Pca priors on W
  //sparseHyperPrior ~ normal(0.01, sparseHyperHyperPrior);
  //for(i in 1:p) W[,i] ~ cauchy(0, sparseHyperPrior);
  
  //PPCA likelihood from Ch. 12 of Kevin Murphy
  Id_n = diag_matrix(rep_vector(1, n_total));
  C = W*W' + sigmaSq*Id_n;
  
  target += -(N/2)*log(determinant(C)) -(N/2)*trace(C\SigmaHat);
}