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
    int n;
    int p;
    int d;
}

parameters {
  real<lower = -pi(), upper = pi()> theta01;
  real<lower = -pi()/2, upper = pi()/2> theta02;
  real<lower = -pi(), upper = pi()> theta12;
}

model {
  matrix[n,p] W;

  vector[3] theta;
  theta[1] = theta01;
  theta[2] = theta02;
  theta[3] = theta12;
  
  W = area_form_lp(theta, n, p);
}
