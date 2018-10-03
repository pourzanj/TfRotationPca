data {
  
}
parameters {
  real x;
  real y;
}
transformed parameters {
  real theta = atan2(x, y)/2.0;
  real r = hypot(x,y);
}
model {
  r ~ normal(1.0, 0.1);
}