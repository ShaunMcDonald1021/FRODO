functions{
  matrix make_RW_mat(int r, int K){
    // Output: matrix RW_mat such that, if X is a (K-1) vector of independent Normals,
    // RW_mat * X is a  Gaussian random walk of order r (with first step omitted, since it
    // equals 0 for identifiability). MUCH faster than defining RW's in loops
    
    // NOTE: currently, FRODO only supports r = 1, 2, 3. Shapes corresponding to r > 3
    // are presumably less common/relevant
    
    int tmp_vec[r]; 
    vector[r] sign_vec;
    vector[r] binom_vec;
    matrix[r, r] coef_mat;
    matrix[K, K] RW_mat_full = diag_matrix(rep_vector(1.0, K));
    matrix[K - 1, K - 1] RW_mat;
    
    if(!(r > 0 && K > 0)){
      reject("r and K must be positive");
    }
    
    for(j in 1:r){
      tmp_vec[j] = j; // tmp_vec = (1, 2, ..., r)
      sign_vec[j] = pow(-1.0, r - tmp_vec[j]); // Alternating signs
      binom_vec[j] = choose(r, tmp_vec[j] - 1); // Binomial coefficients
    }
    coef_mat = diag_matrix(sign_vec .* binom_vec);
    
    for(i in (r+1):K){
      RW_mat_full[i, 1:(i-1)] = rep_row_vector(1.0, r)*(coef_mat * RW_mat_full[(i-r):(i-1), 1:(i-1)]);
    }
    
    // Omit first row and column since first step is 0 in all uses here
    RW_mat = RW_mat_full[2:K, 2:K];
    return RW_mat;
  }
}

data{ 
  // Note: we don't directly require the X-samples
  int<lower = 1> r; // Order of random walk prior on densities
  int<lower = 1> N; // Number of groups
  int<lower = 1> K; // Dimensionality of spline basis/nb. of histogram bins
  matrix[K, N] m; // Number of covariate measurments in each bin, for each group
  vector[N] Y; // Group-level regression responses (presumably standardized)
  // NOTE: currently only supports Gaussian responses
  row_vector[K] eval_points; // Bin midpoints
  int<lower=0> non_zero_nb; // Number of nonzero entries in m matrix
  int non_zero_m_inds[non_zero_nb]; // Indices of nonzero elements of m
  vector[N] delta; // Scale factors for density smoothing parameter priors
  real a; // Left endpoint of density support
  real b; // Right endpoint of density support
  real mu_xi_mean; // Global covariate hypermean (for r = 3)
  vector[K] cent_dens; // Empirical marginal covariate density
}

transformed data{
  matrix[K-1, K-1] dens_RW_mat = make_RW_mat(r, K); // Random walk matrix for densities
  matrix[K-1, K-1] coef_RW_mat = make_RW_mat(2, K); // Random walk matrix for beta
  row_vector[N] zeros = rep_row_vector(0, N); // A vector of zeros
  vector[K*N] m_vec = to_vector(m); // "Flattened" version of m
  vector[non_zero_nb] nonzero_m_vec = m_vec[non_zero_m_inds]; //Nonzero elements of m_vec
}

parameters {
  matrix[K - 1, N] theta_raw; // Un-transformed coefficients for log densities
  vector<lower = 0>[N] tau_dens_raw; // Unscaled ensity smoothing parameters

  real alpha; // Regression intercept
  vector[K-1] beta_raw; // Un-transformed coefficients for regression function
  real<lower = 0> tau_beta_raw; // Smoothing parameter for beta
  
  // Regression scale
  real<lower = 0> sigma_y_1;
  real<lower = 0> sigma_y_2;

  // r = 3: Scale of covariate densities when they are Gaussian
  // r = 2: Unscaled shape parameter for latent group-level covariate rates (alpha_lambda/10 in manuscript)
  real<lower = 0> sigma_X[r > 1];
  
  // r = 3: Scale of latent group-level covariate means
  // r = 2: Mean of latent goup-level covariate rates (mu_lambda in manuscript)
  real<lower = 0> sigma_xi[r > 1];
  
  row_vector[N] xi_raw[r == 3]; // For r = 3, untransformed latent group-level covariate means
  real mu_xi_raw[r == 3]; // For r = 3, untransformed mean of latent xi's
  
  row_vector<lower = 0>[N] lambda_raw[r == 2]; // For r = 2, untransformed latent group-level covariate rates
}

transformed parameters{
  vector[N] tau_dens = delta .* tau_dens_raw;
  
  matrix[K-1, N] theta_free_means = rep_matrix(0, K-1, N); // Means for "free parameters"
  matrix[K, N] theta; // Density coefficients on log scale
  
  matrix[K, N] phi; // Density coefficients
  vector[K*N] phi_vec; // Flattened version of phi
  vector[non_zero_nb] log_phi_vec_nz; // log(phi_vec) components corresponding to nonzero elements of m_vec

  matrix[K, N] phi_cent; // Covariate densities minus empirical central density

  // Regression scale
  real sigma_y = sigma_y_1 * inv_sqrt(sigma_y_2);

  // Coefficients of regression function (before centering)
  vector[K-1] beta0 = append_col(coef_RW_mat[,1], tau_beta_raw*coef_RW_mat[,2:(K-1)])*beta_raw;

  if (r == 3){
    {
    // Latent group-level covariate means
    row_vector[N] xi = mu_xi_mean + mu_xi_raw[1]*15.0*inv_square(K) + sigma_xi[1]*xi_raw[1];
    real tau_X = inv_square(sigma_X[1]);
    
    theta_free_means[1,] = (b-a)*(1.0/K)*(xi - (b-a)*1.0/K - a)*tau_X;
    theta_free_means[2,] = (b-a)*(2.0/K)*(xi - (b-a)*3.0/(2.0*K) - a)*tau_X;
    }  
  } else if(r == 2){
    {
    // Latent group-level covariate rates
    row_vector[N] lambda = lambda_raw[1];
    
    theta_free_means[1,] = -(b-a)*lambda/(1.0*K);
    }
  }
  
  theta = append_row(zeros, dens_RW_mat*(theta_free_means +  diag_post_multiply(theta_raw, tau_dens)));
  
  for(i in 1:N){
    phi[,i] = softmax(theta[,i]);
    phi_cent[,i] = phi[,i] - cent_dens; // "Centered" densities
  }

  phi_vec = to_vector(phi);
  log_phi_vec_nz = log(phi_vec[non_zero_m_inds]);
}

model{
  // LIKELIHOOD
  target += dot_product(nonzero_m_vec, log_phi_vec_nz);
  // Implies m[,i] ~ multinomial(phi[,i]) for i = 1,...,N
  
  Y ~ normal_id_glm(phi_cent[2:K]', alpha, beta0, sigma_y);
  
  // PRIORS
  to_vector(theta_raw) ~ std_normal(); // Implies for all i:
  // theta[2:r,i] ~ normal(theta_free[2:(r-1), i], tau_dens[i])
  // theta[k, i] ~ normal(rth_order_RW(theta[1:(k-1),i], tau_dens[i])) for k > r
  tau_dens_raw ~ exponential(1.0); // Implies tau_dens[i] ~ exponential(1/delta[i]) for all i
  
  if(r == 3){
    sigma_xi[1] ~ std_normal();
    xi_raw[1] ~ std_normal(); // Implies xi ~ normal(mu_xi, sigma_xi)
    mu_xi_raw[1] ~ std_normal(); // Implies mu_xi ~ normal(mu_xi_mean, 15/K^2)
    sigma_X[1] ~ std_normal();
  } else if(r == 2){
    lambda_raw[1] ~ gamma(10*sigma_X[1], 10*(sigma_X[1]/sigma_xi[1]));
    sigma_X[1] ~ std_normal();
    sigma_xi[1] ~ std_normal();
  }
  
  sigma_y_1 ~ normal(0, inv_sqrt(2.0));
  sigma_y_2 ~ gamma(2.0, 2.0); // Implies sigma_y ~ student_t(4, 0, 1/sqrt(2))
  
  beta_raw[1] ~ normal(0, (20.0*(b-a)/K)*sigma_y);
  beta_raw[2:(K-1)] ~ std_normal();
  // Implies beta0[k] ~ normal(2*beta0[k-1] - beta0[k], sigma_y*tau_beta) for k > 1,
  // with beta0[0] = 0
  tau_beta_raw ~ exponential(2.0*inv(sigma_y)); // Implies tau_beta ~ exponential(2)
  
  alpha ~ normal(0, 20.0*sigma_y);
}

generated quantities{
  vector[N] fitted_values = phi_cent[2:K]'*beta0 + alpha; // Estimated regression responses
  row_vector[N] phi_means = eval_points * phi; // Expected values of densities
  vector[K] beta = append_row(0, beta0) - cent_dens'*append_row(0,beta0); // True (centered) beta
  real tau_beta = tau_beta_raw * inv(sigma_y); // True regression smoothing param
}
