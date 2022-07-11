data {
  int<lower = 1> N; // Number of groups
  vector[N] Y; // Group-level regression responses
  vector[N] X; // "Aggregated" covariates (e.g. sample means of X's across groups)
}

parameters {
   real alpha_raw;
   real beta_raw;
   real<lower = 0> sigma_Y_1;
   real<lower = 0> sigma_Y_2;
}

transformed parameters{
  real sigma_Y = sigma_Y_1*inv_sqrt(2.0*sigma_Y_2);
  real alpha = 20.0*alpha_raw;
  real beta = 20.0*beta_raw;
}

model{
  // Likelihood
  Y ~ normal(alpha + beta*X, sigma_Y);

  // Priors
  sigma_Y_1 ~ std_normal(); // Implies sigma_Y ~ student_t(4, 0, sqrt(0.5))
  sigma_Y_2 ~ gamma(2,2);
  beta_raw ~ normal(0, sigma_Y); // Implies beta ~ normal(0, 20*sigma_Y)
  alpha_raw ~ normal(0, sigma_Y); // Implies alpha ~ normal(0, 20*sigma_Y)
}

generated quantities{
  vector[N] fitted_values = alpha + beta*X;
}
