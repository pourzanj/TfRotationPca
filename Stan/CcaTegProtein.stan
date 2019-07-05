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
    
    matrix[] generate_forward_pgivens(vector principal_angles, vector lower_angles, int n, int p) {
        int d = n*p - p*(p+1)/2;
        int partial_idx;
        int lower_idx;
        matrix[n, n] G;
        matrix[n, n] partial_givens[d+1];
        matrix[n, n] R;
        
        //if we're doing full rank i.e. n = p. The last column has no
        //degrees of freedom and therefore no angle so we shouldn't include
        //it in the string of rotation matrix multiplications
        int pp;
        if(p == n) pp = p - 1;
        else pp = p;
        
        G = diag_matrix(rep_vector(1, n));
        partial_idx = 1;
        lower_idx = 1;
        partial_givens[1] = G;
        
        for(i in 1:pp){
          for(j in i+1:n){
            //the first angle of the column is a principal angle
            //and is separate from the lower angles because we might need it
            //to go from (-pi,pi) rather than (-pi/2,pi/2)
            if(j == i+1) {
              G = right_rotation(G, principal_angles[i], n, i, j);
            }
            else {
              G = right_rotation(G, lower_angles[lower_idx], n, i, j);
              lower_idx = lower_idx + 1;
            }
              
            partial_givens[partial_idx + 1] = G;
            partial_idx = partial_idx + 1;  
          }
        }
        
        return partial_givens;      
    }
    
    matrix[] generate_reverse_pgivens(vector principal_angles, vector lower_angles, int n, int p) {
        int d = n*p - p*(p+1)/2;
        int partial_idx;
        int lower_idx;
        matrix[n, n] G_eye;
        matrix[n, p] G;
        matrix[n, p] partial_givens[d+1];
        matrix[n, n] R;
        
        //if we're doing full rank i.e. n = p. The last column has no
        //degrees of freedom and therefore no angle so we shouldn't include
        //it in the string of rotation matrix multiplications
        int pp;
        if(p == n) pp = p - 1;
        else pp = p;

        G_eye = diag_matrix(rep_vector(1, n));
        G = G_eye[,1:p];
        
        partial_givens[d+1] = G;
        partial_idx = d;
        lower_idx = d-pp;
        
        for(i in 1:pp){
          //stan can't do for loops backward so use i_st to go backwards through columns
          int i_st = pp - i + 1;
          
          for(j in i_st+1:n){
            //the first angle of the column is a principal angle
            //and is separate from the lower angles because we might need it
            //to go from (-pi,pi) rather than (-pi/2,pi/2)
            if((n - j + i_st + 1) == i_st+1) {
              G = left_rotation(G, principal_angles[i_st], n, p, i_st, n - j + i_st + 1);
            }
            else {
              G = left_rotation(G, lower_angles[lower_idx], n, p, i_st, n - j + i_st + 1);
              lower_idx = lower_idx - 1;
            }
              
            partial_givens[partial_idx] = G;
            partial_idx = partial_idx - 1; 
          }
        }
        
        return partial_givens;      
    }
    
    matrix[] generate_givens_jacobians(matrix[] partial_givens_forward, matrix[] partial_givens_reverse, vector principal_angles, vector lower_angles, int n, int p) {
        int d = n*p - p*(p+1)/2;
        matrix[n,p] derivative_list[d];
        matrix[n,d] givens_jacobian[p];
        int d_idx = 1;
        int lower_idx = 1;
        matrix[n,n] a;
        matrix[n,n] dR;
        matrix[n,p] b;
        
        //if we're doing full rank i.e. n = p. The last column has no
        //degrees of freedom and therefore no angle so we shouldn't include
        //it in the string of rotation matrix multiplications
        int pp;
        if(p == n) pp = p - 1;
        else pp = p;
        
        //for each angle get derivatve by replacing it's associated matrix
        //by the derivative matrix and carrying out string of multiplies
        for(i in 1:pp){
            for(j in i+1:n){
              
                //the first angle of the column is a principal angle
                //and is separate from the lower angles because we might need it
                //to go from (-pi,pi) rather than (-pi/2,pi/2)
                if(j == i+1) {
                  dR = d_rotation_matrix(principal_angles[i], n, i, j);
                }
                else {
                  dR = d_rotation_matrix(lower_angles[lower_idx], n, i, j);
                  lower_idx = lower_idx + 1;
                }
                
                a = partial_givens_forward[d_idx];
                b = partial_givens_reverse[d_idx + 1];
                
                derivative_list[d_idx] = a * dR * b;
                d_idx = d_idx + 1;
            }
        }
        
        //slice appropriately to get a givens_jacobian which is partial deriv.
        //of each column w.r.t. all angles
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
    
    matrix area_form_lp(vector principal_angles, vector lower_angles, int n, int p) {
        int d = n*p - p*(p+1)/2;
        int idx;
        matrix[n, n] givens;
        matrix[n, n] partial_givens_forward[d+1];
        matrix[n, p] partial_givens_reverse[d+1];
        matrix[n, d] givens_jacobians[p];
        matrix[d, d] area_mat;
        
        //if we're doing full rank i.e. n = p. The last column has no
        //degrees of freedom and therefore no angle so we shouldn't include
        //it in the string of rotation matrix multiplications
        int pp;
        if(p == n) pp = p - 1;
        else pp = p;

        partial_givens_forward = generate_forward_pgivens(principal_angles, lower_angles, n, p);
        partial_givens_reverse = generate_reverse_pgivens(principal_angles, lower_angles, n, p);
        givens = partial_givens_forward[d+1];
        
        givens_jacobians = generate_givens_jacobians(partial_givens_forward, partial_givens_reverse, principal_angles, lower_angles, n, p);
        
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

  int p_common;
  
  int d_W_0;
  int d_B_0;

  int d_W_1;
  int d_B_1;
  
  real sigmaSqHyperPrior;
  
  //real<lower =0> ardHyperHyperPrior;
  //real<lower =0> sparseHyperHyperPrior;
  
  int N; //num of observations
  matrix[n_0+n_1,n_0+n_1] SigmaHat;
}

parameters {
  //Proteins
  vector<lower = -pi()/2, upper = pi()/2>[min(p_common, n_0-1)] theta_principal_W_0;
  vector<lower = -pi()/2, upper = pi()/2>[d_W_0-min(p_common, n_0-1)] theta_lower_W_0;
  positive_ordered[p_common] lambda_0_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[min(p_0, n_0-1)] theta_principal_B_0;
  vector<lower = -pi()/2, upper = pi()/2>[d_B_0-min(p_0, n_0-1)] theta_lower_B_0;
  positive_ordered[p_0] gamma_0_reversed;
  
  //Teg
  vector<lower = -pi(), upper = pi()>[min(p_common, n_1-1)] theta_principal_W_1;
  vector<lower = -pi()/2, upper = pi()/2>[d_W_1-min(p_common, n_1-1)] theta_lower_W_1;
  positive_ordered[p_common] lambda_1_reversed;
  
  vector<lower = -pi()/2, upper = pi()/2>[min(p_1, n_1-1)] theta_principal_B_1;
  vector<lower = -pi()/2, upper = pi()/2>[d_B_1-min(p_1, n_1-1)] theta_lower_B_1;
  positive_ordered[p_1] gamma_1_reversed;
  
  //hyper-priors
  real<lower = 0> sigmaSq;
  
  real<lower = 0> ardHyperPrior;
  real<lower = 0> sparseHyperPrior;
}

transformed parameters{
  
  matrix[n_0, p_common] W_0;
  matrix[n_0, p_0] B_0;
  
  matrix[n_1, p_common] W_1;
  matrix[n_1, p_1] B_1;
  
  vector<lower=0>[p_common] lambda_0;
  vector<lower=0>[p_0] gamma_0;

  vector<lower=0>[p_common] lambda_1;
  vector<lower=0>[p_1] gamma_1;
  
  for (i in 1:p_common) lambda_0[i] = lambda_0_reversed[p_common - i + 1];
  for (i in 1:p_0) gamma_0[i] = gamma_0_reversed[p_0 - i + 1];

  for (i in 1:p_common) lambda_1[i] = lambda_1_reversed[p_common - i + 1];
  for (i in 1:p_1) gamma_1[i] = gamma_1_reversed[p_1 - i + 1];

  W_0 = area_form_lp(theta_principal_W_0, theta_lower_W_0, n_0, p_common);
  B_0 = area_form_lp(theta_principal_B_0, theta_lower_B_0, n_0, p_0);

  W_1 = area_form_lp(theta_principal_W_1, theta_lower_W_1, n_1, p_common);
  B_1 = area_form_lp(theta_principal_B_1, theta_lower_W_1, n_1, p_1);

}

model {
  matrix[n_0 + n_1, p_common + p_0 + p_1] W;
  matrix[n_0 + n_1, n_0 + n_1] Id_n;
  matrix[n_0 + n_1, n_0 + n_1] C;
  
  W = rep_matrix(0, n_0 + n_1, p_common + p_0 + p_1);

  // Construct Row Wise
  W[1:n_0, 1:p_common] = W_0 * diag_matrix(lambda_0);
  W[n_0 + 1:n_0 + n_1, 1:p_common] = W_1 * diag_matrix(lambda_1);

  W[1:n_0, p_common + 1:p_common + p_0] = B_0 * diag_matrix(gamma_0);
  W[n_0 + 1:n_0 + n_1, p_common + p_0 + 1:p_common + p_0 + p_1] = B_1 * diag_matrix(gamma_1);
  
  //priors
  sigmaSq ~ normal(1, sigmaSqHyperPrior);
  
  //control on lambda size
  //ardHyperPrior ~ normal(0.1, 10);
  
  lambda_0_reversed ~ cauchy(0, ardHyperPrior);
  lambda_1_reversed ~ cauchy(0, ardHyperPrior);
  
  gamma_0_reversed ~ cauchy(0, ardHyperPrior);
  gamma_1_reversed ~ cauchy(0, ardHyperPrior);
  
  //sparse Pca priors on W
  //sparseHyperPrior ~ normal(0.1, 10);

  for(i in 1:p_common) {
    W_0[,i] ~ cauchy(0, sparseHyperPrior);
    W_1[,i] ~ cauchy(0, sparseHyperPrior);
  }
  
  for(i in 1:p_0) B_0[,i] ~ cauchy(0, sparseHyperPrior);
  for(i in 1:p_1) B_1[,i] ~ cauchy(0, sparseHyperPrior);
  
  //PPCA likelihood from Ch. 12 of Kevin Murphy
  Id_n = diag_matrix(rep_vector(1, n_0 + n_1));
  C = W*W' + sigmaSq*Id_n;
  
  target += -(N/2)*log(determinant(C)) -(N/2)*trace(C\SigmaHat);
}
