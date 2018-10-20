functions {
  
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
  
  matrix[] generate_reverse_partial_givens(int N, int P, vector theta) {
    
    int D = N*P - P*(P+1)/2;
    int d = D; // index to track which angle we're on
    
    matrix[N, P] G = diag_matrix(rep_vector(1, N))[,1:P];
    matrix[N, P] B[D+1];
    
    // if N == P the for loop shouldn't iterate all the way to the last column
    // because there is no angle there.
    int P_;
    if(P == N) P_ = P-1;
    else P_ = P;
    
    // start from the back which is just I_np
    B[D+1] = G;
    for(p in 1:P_){
      for(n in (p+1):N){
        
        int p_rev = P_-p+1;
        int n_rev = N-n+p+1;
        
        // apply left rotation to G
        row_vector[P] G_i_temp = cos(theta[d])*G[p_rev,] - sin(theta[d])*G[n_rev,];
        row_vector[P] G_j_temp = sin(theta[d])*G[p_rev,] + cos(theta[d])*G[n_rev,];
        G[p_rev,] = G_i_temp;
        G[n_rev,] = G_j_temp;
        
        // store in B then go to next rotation matrix on the left
        B[d] = G;
        d = d - 1;
      }
    }
    
    return B;      
  }
  
  
  matrix givens_lp(int N, int P, vector theta) {

    int D = N*P - P*(P+1)/2;
    int d = 1; // index to track which angle we're on
    
    matrix[N,P] B[D+1]; 
    matrix[N,N] G = diag_matrix(rep_vector(1, N)); // the eventual product of rotation matrices
    matrix[N,P] J[D]; // derivative of the final output w.r.t. each of the D angles
    matrix[N,D] J_[P];
    matrix[D,D] A; // matrix of one forms
    
    // if N == P the for loop shouldn't iterate all the way to the last column
    // because there is no angle there.
    int P_;
    if(P == N) P_ = P-1;
    else P_ = P;
    
    //////////////////////////
    // 1. Compute and store partial matrix products of the rotations starting from the
    // last matrix i.e. compute I_np, R_d*I_np, R_d-1*R_d*I_np 
    B = generate_reverse_partial_givens(N,P, theta);
    
    //////////////////////////
    // 2. Simultaneously compute rotation matrix products that comprise G and also
    // compute the gradient of the final N x P matrix w.r.t. to each input and store in a tensor
    // This gradient is simply the product of the rotation matrices with the rotation matrix corresponding
    // the theta_ij swapped out for the corresponding derivative matrix
    for(p in 1:P_) {
      for(n in (p+1):N) {
        
        vector[N] G_i_temp;
        vector[N] G_j_temp;

        // compute gradient term
        // we need to do a * dR * b where b comes from partial givens reverse and a
        // is the successive partial multiplications of G
        // dR is mostly sparse with just 4 elements non-zero at locations, (i,i),(i,j),(j,i),(j,j)
        // thus dR * B is an n x p matrix with only the ith and jth rows non-zero.
        row_vector[P] dR_b_row_i = -sin(theta[d])*B[d + 1][p,] - cos(theta[d])*B[d + 1][n,];
        row_vector[P] dR_b_row_j = cos(theta[d])*B[d + 1][p,] - sin(theta[d]*B[d + 1][n,]);

        // since dR*B is also sparse this simplifies G * dR * B
        for(pp in 1:P_) {
          J[d,,pp] = dR_b_row_i[pp]*G[,p] + dR_b_row_j[pp]*G[,n];
        }
        
        // compute G by multiplying on the left by the rotation matrix. this is equivalent
        // to just updating the columns
        G_i_temp = cos(theta[d])*G[,p] + sin(theta[d])*G[,n];
        G_j_temp = -sin(theta[d])*G[,p] + cos(theta[d])*G[,n];
        G[,p] = G_i_temp;
        G[,n] = G_j_temp;
        
        // update angle index
        d = d + 1;
      }
    }
    
    //////////////////////////
    // 3. Compute one forms using matrix multiplication and use them to fill a D x D matrix
    // J has to be rearranged in to matrix form
    for(p in 1:P) {
      for(dd in 1:D) {
        J_[p,,dd] = J[dd,,p];
      }
    }
    
    d = 1;
    for(p in 1:P){
      matrix[D, N-p] one_forms = ((G[,(p+1):N])' * J_[p])';
      A[,d:(d+N-p-1)] = one_forms;
      d = d + N-p;
    }

    //////////////////////////
    // 4. increment lp with log determinant of one form matrix and return final orthogonal N x P matrix
    target += log_determinant(A);

    return(B[1]);
  }
}

data {
  int N;
  int P;
  
  int Y[N,N];
}
transformed data {
  int D = N*P - P*(P+1)/2;
}
parameters {
  vector[D] x;
  vector[D] y;
  
  positive_ordered[P] L_rev;
  real c;
}
transformed parameters{
  vector[D] theta_raw;
  vector[D] theta;
  vector[D] r;
  vector[P] L;
  
  matrix[N,P] U;
  for(p in 1:P) {
    L[p] = L_rev[P-p+1];
  }

  // transform donut
  for(d in 1:D) {
    theta_raw[d] = atan2(y[d], x[d]);
    theta[d] = mirror_atan2(theta_raw[d]);
    r[d] = hypot(x[d],y[d]); 
  }
  
  U = givens_lp(N,P,theta);
}
model {
  matrix[N,N] Probs = Phi(quad_form(diag_matrix(L), U') + c);
  r ~ normal(1.0, 0.1);

  L ~ normal(0, P);
  c ~ normal(0, 10);
  
  // we only want to count each connection once so only
  // include entries of connecitvity matrix below diagonal
  for(n in 1:N) {
    Y[(n+1):N,n] ~ bernoulli(Probs[(n+1):N,n]);
  }
}
