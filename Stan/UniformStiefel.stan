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
        
        G = diag_matrix(rep_vector(1, n));
        idx = 1;
        partial_givens[1] = G;
        for(i in 1:p){
            for(j in i+1:n){
                R = rotation_matrix(angles[idx], n, i, j);
                G = G * R;
                //G = right_rotation(G, angles[idx], n, i, j);
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

        G_eye = diag_matrix(rep_vector(1, n));
        G = G_eye[,1:p];
        
        partial_givens[d+1] = G;
        idx = d;
        for(i in 1:p){
            int i_st = p - i + 1;
            for(j in i_st+1:n){
                R = rotation_matrix(angles[idx], n, i_st, n - j + i_st + 1);
                G = R * G;
                //G = left_rotation(G, angles[idx], n, p, i_st, n - j + i_st + 1);
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
        
        for(i in 1:p){
            for(j in i+1:n){
                matrix[n,n] dR = d_rotation_matrix(angles[idx], n, i, j);
                matrix[n,n] a = partial_givens_forward[idx];
                matrix[n,p] b = partial_givens_reverse[idx + 1];
                
                derivative_list[idx] = a * dR * b;
                idx = idx + 1;
            }
        }
        
        for(i in 1:p) {
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
        /**
         * Create Partial Givens
        */
        
        //partial_givens_forward = generate_forward_pgivens(angles, n, p);
        //partial_givens_reverse = generate_reverse_pgivens(angles, n, p);
        givens = partial_givens_forward[d+1];
        
        givens_jacobians = generate_givens_jacobians(partial_givens_forward, partial_givens_reverse, angles, n, p);
        
        idx = 1;
        for(i in 1:p){
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
    int n;
    int p;
    int d;
}
parameters {
  vector<lower = -pi()/2, upper = pi()/2>[d] theta;
}
model {
  matrix[n, n] partial_givens_forward[d+1];
  matrix[n, p] partial_givens_reverse[d+1];
  matrix[n, n] G;
  matrix[n, p] W;

  //add Stiefel "area form" to log probability since we are sampling on angle space
  //this requires partial matrix multiplications of rotation matrices
  partial_givens_forward = generate_forward_pgivens(theta, n, p);
  partial_givens_reverse = generate_reverse_pgivens(theta, n, p);
  target += area_form(partial_givens_forward, partial_givens_reverse, theta, n, p);
  
  //last "partial" multiplication is actually full multiplication that gives us
  //nxn orthonormal Givens matrix. From it we can slice out nxp orthornomal W
  G = partial_givens_forward[d+1];
  W = G[,1:p];
}