functions {
  matrix custom_matmul(matrix A, matrix B);
}
data {
  int n;
  int p;
}
transformed data {
  int d =  n*p-p*(p+1)/2;
}
parameters {
  matrix[n, n] G;
  matrix[n, d] J[p];
}
model {
  
  matrix[d,d] A;
  int idx = 1;
  
  to_array_1d(G) ~ normal(0, 1);
  for(i in 1:p) {
    to_array_1d(J[p]) ~ normal(0, 1);
  }
  
  for(i in 1:p){
    
    matrix[n, n-i] G_cols = G[,i+1:n];
    matrix[n-i, n] G_cols_T = G_cols';
      
    matrix[n-i, d] one_forms = G_cols_T * J[p];
    matrix[d, n-i] one_forms_T = one_forms';
    
    A[,idx:(idx+n-i-1)] = one_forms_T;
    idx = idx + n-i;
  }

  target += log_determinant(A);
}
