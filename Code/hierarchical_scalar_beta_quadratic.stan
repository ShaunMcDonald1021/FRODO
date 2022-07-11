data {
  int<lower = 1> N; // Number of groups
  int<lower = 1> M; // Total number of covariate observations (sum of n_i's)
  vector[N] Y; // Group-level regression responses
  vector[M] X; // Individual-level covariate observations for all groups
  int n[N]; // Number of individuals per group;
}

parameters {
  real alpha_raw;
  real beta_raw;
  real<lower = 0> sigma_Y_1;
  real<lower = 0> sigma_Y_2;
  vector<lower=0, upper = 3>[N] xi;
}

transformed parameters{
  real sigma_Y = sigma_Y_1*inv_sqrt(2.0*sigma_Y_2);
  real alpha = 10.0*alpha_raw;
  real beta = 10.0*beta_raw;
  vector[N] true_covar = 0.25*(1+1 ./(2*xi + 1));
}

model{
  // Likelihoods
  int pos = 1;
  for(i in 1:N){
    segment(X, pos, n[i]) ~ beta(xi[i], xi[i]);
    pos += n[i];
  }
  
  Y ~ normal(alpha + beta*true_covar, sigma_Y);

  // Priors
  sigma_Y_1 ~ std_normal(); // Implies sigma_Y ~ student_t(0, 4, sqrt(0.5))
  sigma_Y_2 ~ gamma(2,2);
  beta_raw ~ normal(0, sigma_Y); // Implies beta ~ normal(0, 10*sigma_Y)
  alpha_raw ~ normal(0, sigma_Y); // Implies alpha ~ normal(0, 10*sigma_Y)
}

generated quantities{
  vector[N] fitted_values = alpha + beta*true_covar;
}