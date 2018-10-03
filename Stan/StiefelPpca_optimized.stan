functions {
  
  // real get_time();
  
  real mirror_atan2(real theta) {
    
    real eps;
    
    // upper plane
    if(theta > 0.0) {
      eps = theta - pi()/2;
      if(eps <= 0.0)
        return(theta);
      else
        return(-pi()/2 + eps);
    }
    // lower plane
    else {
      eps = theta + pi()/2;
      if(eps >= 0.0)
        return(theta);
      else
        return(pi()/2 + eps);
    }
  }
  
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
  
  // this is equivalent to multiplying the matrix A on the left by the matrix
  // that is all zeros but has dR[i, i] = -sin(angle);dR[i, j] = -cos(angle);dR[j, i] = cos(angle);dR[j, j] = -sin(angle);
  matrix d_rotation_matrix_left_apply(matrix A, real angle, int n, int p, int i, int j) {
    matrix[n,p] AR;

    // dR[i, i] = -sin(angle);
    // dR[i, j] = -cos(angle);
    // dR[j, i] = cos(angle);
    // dR[j, j] = -sin(angle);
    // 
    // AR[,i] = -sin(angle)*A[,i] - cos(angle)*A[,j];
    // AR[,j] = -sin(angle)*A[,i] + cos(angle)*A[,j];
    
    return AR;
  }
  
  matrix[] generate_forward_pgivens(vector angles, int n, int p) {
    int d = n*p - p*(p+1)/2;
    int idx;
    matrix[n, n] G;
    matrix[n, n] partial_givens[d+1];
    // matrix[n, n] partial_givens[1];
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
        // G = right_rotation(G, angles[idx], n, i, j);
        vector[n] G_i_temp = cos(angles[idx])*G[,i] + sin(angles[idx])*G[,j];
        G[,j] = -sin(angles[idx])*G[,i] + cos(angles[idx])*G[,j];
        G[,i] = G_i_temp;
        
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
        // matrix[n,n] dR = d_rotation_matrix(angles[idx], n, i, j);
        matrix[n,n] a = partial_givens_forward[idx];
        matrix[n,p] b = partial_givens_reverse[idx + 1];
        
        // derivative_list[idx] = a * dR * b;
        // derivative_list[idx] = partial_givens_forward[idx] * d_rotation_matrix(angles[idx], n, i, j) * partial_givens_reverse[idx + 1];
        // derivative_list[idx] = partial_givens_forward[idx] * d_rotation_matrix_left_apply(partial_givens_reverse[idx + 1], angles[idx], n, i, j);
        
        // dR is mostly sparse with just 4 elements non-zero at locations, (i,i),(i,j),(j,i),(j,j)
        // thus dR * b is an n x p matrix with only the ith and jth rows non-zero. 
        row_vector[p] dR_b_row_i = -sin(angles[idx])*b[i,] - cos(angles[idx])*b[j,];
        row_vector[p] dR_b_row_j = cos(angles[idx])*b[i,] - sin(angles[idx]*b[j,]);
        
        // dR*b is all zeros with only the ith and jth rows non-zero so multiplication by a is cheaper
        for(k in 1:p) {
          derivative_list[idx,,k] = dR_b_row_i[k]*a[,i] + dR_b_row_j[k]*a[,j];
        }
        
        idx = idx + 1;
      }
    }
    
    for(i in 1:pp) {
      for(j in 1:d) {
        givens_jacobian[i,,j] = derivative_list[j][,i];
      }
    }
    
    // for(i in 1:pp) {
    //   for(j in 1:d) {
    //     vector[n] t = derivative_list[j][,i];
    //     matrix[n,d] z = givens_jacobian[i];
    //     z[,j] = t;
    //     // givens_jacobian[i] = z;
    //   }
    // }
    
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
    // return 1.0;
  }
}

data {
  int n;
  int p;
  int d;
  
  real sigmaSqHyperPrior;
  int N; //num of latent observations
  // int Y[n*n];
}
parameters {
  vector[d] x;
  vector[d] y;
  // vector<lower = -pi()/2, upper = pi()/2>[d] theta;
  // positive_ordered[p] lambdaReversed;
  // real<lower = 0> sigmaSq;
}
transformed parameters{
  vector[d] theta;
  vector[d] theta_mirrored;
  vector[d] r;
  // vector<lower=0>[p] lambdaSq;
  real J;
  // transform donut
  for(i in 1:d) {
    theta[i] = atan2(y[i], x[i]);
    theta_mirrored[i] = mirror_atan2(theta[i]);
    r[i] = hypot(x[i],y[i]); 
  }
  // 
  // vector<lower=0>[p] lambdaSq;
  // for (i in 1:p) lambdaSq[i] = pow(lambdaReversed[p - i + 1], 2);
}
model {
  // int pp;
  // matrix[n, n] partial_givens_forward[d+1];
  matrix[n, p] derivative_list[d];
  matrix[n, d] givens_jacobians[p];
  matrix[n, p] partial_givens_reverse[d+1]; //partial_givens_reverse[i,,j]
  matrix[n, n] G;
  matrix[d,d] area_mat;
  
  matrix[n, p] W;
  matrix[n, n] Id_n;
  matrix[n, n] C;
  
  // avoid center of donut
  r ~ normal(1.0, 0.1);
  // theta[1] ~ normal(pi()/2,0.1);
  
  
  // real begin_reverse_pgivens; real end_reverse_pgivens;
  // real begin_compute_G_and_J; real end_compute_G_and_J;
  // real begin_copy_J; real end_copy_J;
  // real begin_oneforms; real end_oneforms;
  // real begin_determinant; real end_determinant;
  // real begin_model; real end_model;
  
  // begin_reverse_pgivens = get_time();
  partial_givens_reverse = generate_reverse_pgivens(theta_mirrored, n, p);
  // end_reverse_pgivens = get_time();
  // 
  // print("Partial Givens Reverse: ", end_reverse_pgivens - begin_reverse_pgivens);
  // 
  //add Stiefel "area form" to log probability since we are sampling on angle space
  //this requires partial matrix multiplications of rotation matrices
  // partial_givens_forward = generate_forward_pgivens(theta, n, p);
  // begin_compute_G_and_J = get_time();
  {
    int idx;
    matrix[n, n] R;
    
    vector[n] G_i_temp;
    vector[n] G_j_temp;
    
    int pp;
    if(p == n) pp = p-1;
    else pp = p;
    
    G = diag_matrix(rep_vector(1, n));
    idx = 1;

    // built G by right multiplying by the rotation matrices in sequence
    // simultaneously fill the jacobians
    for(i in 1:pp){
      for(j in i+1:n){
        
        // we need to do a * dR * b where b comes from partial givens reverse and a
        // is the successive partial multiplications of G
        
        // dR is mostly sparse with just 4 elements non-zero at locations, (i,i),(i,j),(j,i),(j,j)
        // thus dR * b is an n x p matrix with only the ith and jth rows non-zero.
        row_vector[p] dR_b_row_i = -sin(theta_mirrored[idx])*partial_givens_reverse[idx + 1][i,] - cos(theta_mirrored[idx])*partial_givens_reverse[idx + 1][j,];
        row_vector[p] dR_b_row_j = cos(theta_mirrored[idx])*partial_givens_reverse[idx + 1][i,] - sin(theta_mirrored[idx]*partial_givens_reverse[idx + 1][j,]);

        for(k in 1:p) {
          derivative_list[idx,,k] = dR_b_row_i[k]*G[,i] + dR_b_row_j[k]*G[,j];
        }

        // compute rest of G
        G_i_temp = cos(theta_mirrored[idx])*G[,i] + sin(theta_mirrored[idx])*G[,j];
        G_j_temp = -sin(theta_mirrored[idx])*G[,i] + cos(theta_mirrored[idx])*G[,j];
        G[,i] = G_i_temp;
        G[,j] = G_j_temp;

        idx = idx + 1;
      }
    }
  }
  
  
  // end_compute_G_and_J = get_time();
  // 
  // print("Compute G and J: ", end_compute_G_and_J - begin_compute_G_and_J);
  // 
  // begin_copy_J = get_time();
  // copy derivative list in to a more accessible format
  for(i in 1:p) {
    for(j in 1:d) {
      givens_jacobians[i,,j] = derivative_list[j,,i];
    }
  }
  // end_copy_J = get_time();
  // print("Copy J: ", end_copy_J - begin_copy_J);

  // begin_oneforms = get_time();
  {
    int idx_2 = 1;
    for(i in 1:p){
      matrix[n-i, d] one_forms;
      matrix[d, n-i] one_forms_t;
      matrix[n, n-i] G_partial = G[,i+1:n];
      matrix[n-i, n] G_partial_t = G_partial';
      matrix[n, d] givens_jacobian_temp = givens_jacobians[i];
      one_forms = G_partial_t * givens_jacobian_temp;
      one_forms_t = one_forms';
      // one_forms = ((G[,i+1:n])' * givens_jacobians[i])';
      area_mat[,idx_2:(idx_2+n-i-1)] = one_forms_t;
      idx_2 = idx_2 + n-i;
      // for(j in 1:n-i) {
      //   area_mat[,idx_2] = one_forms[,j];
      //   idx_2 = idx_2 + 1;
      // }
    }
  }

  // end_oneforms = get_time();
  // print("Oneforms: ", end_oneforms - begin_oneforms);

  // begin_determinant = get_time();
  target += log_determinant(area_mat);
  // end_determinant = get_time();
  // print("Determinant: ", end_determinant - begin_determinant);
  
  // begin_model = get_time();
  // target += area_form(partial_givens_forward, partial_givens_reverse, theta, n, p);
  
  //last "partial" multiplication is actually full multiplication that gives us
  //nxn orthonormal Givens matrix. From it we can slice out nxp orthornomal W
  // G = partial_givens_forward[d+1];
  W = G[,1:p];
  //priors
  // sigmaSq ~ normal(1, sigmaSqHyperPrior);
  // 
  // //PPCA likelihood from Ch. 12 of Kevin Murphy
  // Id_n = diag_matrix(rep_vector(1, n));
  // C = Phi_approx(W*diag_matrix(lambdaSq)*W' + sigmaSq*Id_n);
  // 
  // Y ~ bernoulli(to_vector(C));

  // target += -(N/2)*log(determinant(C)) -(N/2)*trace(C\SigmaHat);
  
  // end_model = get_time();
  // print("Model: ", end_determinant - begin_determinant);  
}
