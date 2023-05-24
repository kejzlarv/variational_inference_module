//Logistic regression model for STAN
data {
  int<lower=0> N;
  vector[N] x;
  int<lower=0,upper=1> y[N];
}
parameters {
  real alpha; // intercept
  real beta; // slope
}
model {
  y ~ bernoulli_logit(alpha + beta * x);

  // Priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
}
