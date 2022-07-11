data {
  int<lower = 1> N; // Number of groups
  int<lower = 1> M; // Total number of covariate observations (sum of n_i's)
  vector[N] Y; // Group-level regression responses
  vector[M] X; // Individual-level covariate observations for all groups
  int n[N]; // Number of individuals per group
}

parameters {
  real alpha_raw;
  real beta_raw;
  real<lower = 0> sigma_Y_1;
  real<lower = 0> sigma_Y_2;
  real<lower = 0> sigma_xi;
  real<lower = 0> sigma_X;
  vector[N] xi;
  real mu_xi;
}

transformed parameters{
  real sigma_Y = sigma_Y_1*inv_sqrt(2.0*sigma_Y_2);
  real alpha = 20.0*alpha_raw;
  real beta = 20.0*beta_raw;
}

model{
  // Likelihoods
  int pos = 1;
  for(i in 1:N){
    segment(X, pos, n[i]) ~ normal(xi[i], sigma_X);
    pos += n[i];
  }
  Y ~ normal(alpha + beta*xi, sigma_Y);

  // Priors
  xi ~ normal(mu_xi, sigma_xi);
  sigma_Y_1 ~ std_normal(); // Implies sigma_Y ~ student_t(0, 4, sqrt(0.5))
  sigma_Y_2 ~ gamma(2,2);
  sigma_xi ~ std_normal(); 
  sigma_X ~ std_normal();
  beta_raw ~ normal(0, sigma_Y); // Implies beta ~ normal(0, 20*sigma_Y)
  alpha_raw ~ normal(0, sigma_Y); // Implies alpha ~ normal(0, 20*sigma_Y)
  mu_xi ~ std_normal();
}

generated quantities{
  vector[N] fitted_values = alpha + beta*xi;
}
