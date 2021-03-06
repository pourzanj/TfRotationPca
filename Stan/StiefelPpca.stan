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
  
  // counter-clockwise rotation of columns
  // A %*% R where R = [[c, -s],
  //                    [s, c]]
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
  int n;
  int p;
  int d;
  
  real sigmaSqHyperPrior;
  int N; //num of latent observations
  matrix[n,n] SigmaHat;
}
parameters {
  vector[d] x;
  vector[d] y;
  
  // positive_ordered[p] lambdaReversed;
  // real<lower = 0> sigmaSq;
}
transformed parameters{
  
  vector[d] theta;
  vector[d] r;
  real J;
  // vector<lower=0>[p] lambdaSq;
  
  // transform donut
  for(i in 1:d) {
    theta[i] = atan2(x[i], y[i])/2.0;
    r[i] = hypot(x[i],y[i]); 
  }
  
  // reverse lambda
  // for (i in 1:p) lambdaSq[i] = pow(lambdaReversed[p - i + 1], 2);
}
model {
  int pp;
  matrix[n, n] partial_givens_forward[d+1];
  matrix[n, p] partial_givens_reverse[d+1];
  matrix[n, n] G;
  matrix[n, p] W;
  matrix[n, n] Id_n;
  matrix[n, n] C;
  
  
  // avoid center of donut
  r ~ normal(1.0, 0.1);
  // theta ~ normal(pi()/2,0.1);
  
  //add Stiefel "area form" to log probability since we are sampling on angle space
  //this requires partial matrix multiplications of rotation matrices
  partial_givens_forward = generate_forward_pgivens(theta, n, p);
  partial_givens_reverse = generate_reverse_pgivens(theta, n, p);
  J = area_form(partial_givens_forward, partial_givens_reverse, theta, n, p);
  target += J;
  
  //last "partial" multiplication is actually full multiplication that gives us
  //nxn orthonormal Givens matrix. From it we can slice out nxp orthornomal W
  G = partial_givens_forward[d+1];
  W = G[,1:p];
  
  //priors
  // sigmaSq ~ normal(1, sigmaSqHyperPrior);
  
  //PPCA likelihood from Ch. 12 of Kevin Murphy
  Id_n = diag_matrix(rep_vector(1, n));
  C = W*diag_matrix(lambdaSq)*W' + sigmaSq*Id_n;

  target += -(N/2)*log(determinant(C)) -(N/2)*trace(C\SigmaHat);
}