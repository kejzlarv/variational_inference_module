//Gamma~Poisson model for STAN
data {
  int N; // number of observations
  int y[N]; // discrete valued observations
}
parameters {
  // The unknown rate parameter
  real lambda_real;
}
transformed parameters{
  real<lower=0> lambda;
  lambda = exp(lambda_real);
}
model {
  // prior is on the transformed parameters, therefore we need to have  jacobian adjustment
  target += gamma_lpdf(lambda | 100, 2);
  target += lambda_real;
  // likelihood
  for (n in 1:N){
    target += poisson_lpmf(y[n]|lambda);
  }
}
