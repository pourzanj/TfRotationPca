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
    
    matrix area_form(vector angles, int n, int p) {
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
    
        return givens[,1:p];
    }
}

data {
    int n;
    int p;
    int d;
    int N; // number of data points
    
    int<lower=1> K; // num categories
    int<lower=0> interactions[N,n,n];
}

transformed data {
}

parameters {
    simplex[K] phi[K]; //transit probs
    vector<lower = -pi()/2, upper = pi()/2>[d] theta[K]; // emission parameters
    vector[p] Lambda[K];
    vector<lower=0>[K] beta;
}

transformed parameters{
    matrix[n, n] Ws[K];
    matrix[n, p] Wt[K];
    
    for (k in 1:K) {
        matrix[p, p] diag_l;
        diag_l = diag_matrix(Lambda[k]);
        Wt[k] = area_form_lp(theta[k], n, p);
        Ws[k] = Wt[k] * diag_l * Wt[k]';
    }
}

model {
    real acc[K];
    real gamma[N, K];
    
    for (k in 1:K)
        phi[k] ~ dirichlet(beta);
    
    // compute jacobian correction for area representation
    // and combined matrix
    for (k in 1:K) {
        gamma[1,k] = 0;
        for(x in 1:n) {
            for(y in x:n) {
                gamma[1,k] = gamma[1,k] + poisson_lpmf(interactions[1,x,y] | exp(Ws[k][x,y]));
            }
        }
    }
    
    // forwards algorithm to compute likelihood
    for (t in 2:N) {
        for (k in 1:K) {
            for (j in 1:K) {
                acc[j] = gamma[t-1,j] + log(phi[j,k]);
                for(x in 1:n) {
                    for(y in x:n) {
                        acc[j] = acc[j] + poisson_lpmf(interactions[t,x,y] | exp(Ws[k][x,y]));
                    }
                }
            }
            gamma[t,k] = log_sum_exp(acc);
        }
    }
    
    // accumulate likelihood
    target += log_sum_exp(gamma[N]);
}

generated quantities {
    int<lower=1, upper=K> y_star[N];
    real log_p_y_star;
    {
        int back_ptr[N, K];
        real best_logp[N, K];
        for (k in 1:K) {
            best_logp[1,k] = 0;
            for(x in 1:n) {
                for(y in x:n) {
                    best_logp[1,k] = best_logp[1,k] + poisson_lpmf(interactions[1,x,y] | exp(Ws[k][x,y]));
                }
            }
        }
    
        for (t in 2:N) {
            for (k in 1:K) {
                best_logp[t,k] = negative_infinity();
                for(j in 1:K) {
                    real logp;
                    logp = best_logp[t-1,j] + log(phi[j,k]);
                    for(x in 1:n) {
                        for(y in x:n) {
                            logp = logp + poisson_lpmf(interactions[t,x,y] | exp(Ws[k][x,y]));
                        }
                    }
                    if (logp > best_logp[t,k]) {
                        back_ptr[t,k] = j;
                        best_logp[t,k] = logp;
                    }

                }
            }
        }

        log_p_y_star = max(best_logp[N]);
        for(k in 1:K) {
            if(best_logp[N,k] == log_p_y_star) {
                y_star[N] = k;
            }
        }

        for(t in 1:(N-1)){
            y_star[N - t] = back_ptr[N - t + 1, y_star[N - t + 1]];
        }
    } 
}
                
