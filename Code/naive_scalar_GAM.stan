functions{
  matrix make_RW_mat(int r, int K){
    // Output: matrix RW_mat such that, if X is a (K-1) vector of independent Normals,
    // RW_mat * X is a random walk of order r (with first step omitted, since it
    // equals 0 for identifiability). MUCH faster than defining RW's in loops
    
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
  
  // Function to evaluate spline basis functions, taken from this case study:
  // https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}


data {
  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;
  // Above lines also taken from case study
  int<lower = 1> N; // Number of groups
  vector[N] Y; // Group-level regression responses
  real X[N]; // "Aggregated" covariates (e.g. sample means of X's across groups)
  int num_eval_points; // Size of grid on which to evaluate basis functions
  real eval_points[num_eval_points]; // Grid on which to evaluate basis functions
}

transformed data{
  // Taken from case study and modified not to use repeated knots at endpoints
  // so proper P-spline basis is ensured (see "Bayesian P-Splines", Lang & Brezger, 2004)
  int num_basis = num_knots + spline_degree - 1;
  int K = num_basis - 6; 
  matrix[K-1, K-1] coef_RW_mat = make_RW_mat(2, K);
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  matrix[K, N] B;
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 4:(num_basis - 3))
    B[ind-3,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
}

parameters {
 real alpha_raw;
 vector[K-1] beta_raw;
 real<lower = 0> sigma_Y_1;
 real<lower = 0> sigma_Y_2;
  real<lower = 0> tau_beta; // Spline smoothing parameter
}

transformed parameters{
  real sigma_Y = sigma_Y_1*inv_sqrt(2.0*sigma_Y_2);
  real alpha = 20.0*alpha_raw;
  vector[K-1] beta0 = append_col(20.0*coef_RW_mat[,1], tau_beta*coef_RW_mat[,2:(K-1)])*beta_raw;
  // Ensure spline coefficients sum to 0 for identifiability
  vector[K] beta = append_row(0, beta0) - sum(append_row(0, beta0))/K;
  
}

model{
  // Likelihood
  Y ~ normal_id_glm(B', alpha, beta, sigma_Y);
  
  // Priors
  tau_beta ~ exponential(2);
  sigma_Y_1 ~ std_normal(); // Implies sigma_Y ~ student_t(0, 4, sqrt(0.5))
  sigma_Y_2 ~ gamma(2,2);
  beta_raw ~ normal(0, sigma_Y); // Implies beta ~ normal(0, 20*sigma_Y)
  alpha_raw ~ normal(0, sigma_Y); // Implies alpha ~ normal(0, 20*sigma_Y)
}

generated quantities{
  vector[N] fitted_values = alpha + B'*beta;
  
  // Modified from case study code to evaluate regression function on fine grid
  matrix[K, num_eval_points] Bmat;
  for (ind in 4:(num_basis - 3))
    Bmat[ind-3,:] = to_row_vector(build_b_spline(eval_points, to_array_1d(ext_knots), ind, spline_degree + 1));
}
