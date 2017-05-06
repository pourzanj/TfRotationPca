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
  int n;
  int p;
  
  int d;
  
  real sigmaSqHyperPrior;
  //real<lower=0> ardHyperHyperPrior;
  //real<lower=0> sparseHyperHyperPrior;
  
  int N; //num of observations
  
  real lowerAngle2; real upperAngle2;
  real lowerAngle3; real upperAngle3;
  real lowerAngle64; real upperAngle64;
  real sigma_hier_hyper;
  
  matrix[n,n] SigmaHat2;
  matrix[n,n] SigmaHat3;
  matrix[n,n] SigmaHat64;
}

parameters {
  vector<lower = lowerAngle2, upper = upperAngle2>[min(p, n-1)] theta_principal2;
  vector<lower = -pi()/2, upper = pi()/2>[d-min(p, n-1)] theta_lower2;
  positive_ordered[p] lambda_reversed2;
  real<lower = 0> sigmaSq2;
  
  vector<lower = lowerAngle3, upper = upperAngle3>[min(p, n-1)] theta_principal3;
  vector<lower = -pi()/2, upper = pi()/2>[d-min(p, n-1)] theta_lower3;
  positive_ordered[p] lambda_reversed3;
  real<lower = 0> sigmaSq3;
  
  vector<lower = lowerAngle64, upper = upperAngle64>[min(p, n-1)] theta_principal64;
  vector<lower = -pi()/2, upper = pi()/2>[d-min(p, n-1)] theta_lower64;
  positive_ordered[p] lambda_reversed64;
  real<lower = 0> sigmaSq64;
  
  real<lower = 0, upper = pi()/2> mu_hier;
  real<lower = 0> sigma_hier;
}

transformed parameters{
  
  matrix[n, p] W2; matrix[n, p] W3; matrix[n, p] W64;
  vector<lower=0>[p] lambdaSq2; vector<lower=0>[p] lambdaSq3; vector<lower=0>[p] lambdaSq64;
  
  for (i in 1:p) lambdaSq2[i] = pow(lambda_reversed2[p - i + 1], 2);
  for (i in 1:p) lambdaSq3[i] = pow(lambda_reversed3[p - i + 1], 2);
  for (i in 1:p) lambdaSq64[i] = pow(lambda_reversed64[p - i + 1], 2);
  W2 = area_form_lp(theta_principal2, theta_lower2, n, p);
  W3 = area_form_lp(theta_principal3, theta_lower3, n, p);
  W64 = area_form_lp(theta_principal64, theta_lower64, n, p);

}

model {

  matrix[n, n] Id_n;
  matrix[n, n] C2; matrix[n, n] C3; matrix[n, n] C64;
  
  sigmaSq2 ~ normal(1, sigmaSqHyperPrior);
  sigmaSq3 ~ normal(1, sigmaSqHyperPrior);
  sigmaSq64 ~ normal(1, sigmaSqHyperPrior);
  
  Id_n = diag_matrix(rep_vector(1, n));
  C2 = W2*diag_matrix(lambdaSq2)*W2' + sigmaSq2*Id_n;
  C3 = W3*diag_matrix(lambdaSq3)*W3' + sigmaSq3*Id_n;
  C64 = W64*diag_matrix(lambdaSq64)*W64' + sigmaSq64*Id_n;
  
  target += -(N/2)*log(determinant(C2)) -(N/2)*trace(C2\SigmaHat2);
  target += -(N/2)*log(determinant(C3)) -(N/2)*trace(C3\SigmaHat3);
  target += -(N/2)*log(determinant(C64)) -(N/2)*trace(C64\SigmaHat64);
  
  sigma_hier ~ normal(0.1, sigma_hier_hyper);
  theta_principal2 ~ normal(mu_hier, sigma_hier);
  theta_principal3 ~ normal(mu_hier, sigma_hier);
  theta_principal64 ~ normal(mu_hier, sigma_hier);
}
