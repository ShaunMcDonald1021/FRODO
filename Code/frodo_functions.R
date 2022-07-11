# Dependencies: gtools, mgcv, rstan, fda

package_get = function(package_name, install_fun = function(c) install.packages(package_name)){
  # Function to automatically attach packages from library
  # If they're not already installed, it'll install them, THEN attach
  
  # package_name: name of package as character object
  # install_fun: function to install package if not already installed
  
  tryCatch(library(package_name, character.only = TRUE), error = install_fun,
           finally = library(package_name, character.only = TRUE))
}

# Get all dependencies
package_get('fda')
package_get('rstan')
pkgbuild::has_build_tools(debug = TRUE)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
package_get('gtools')
package_get('mgcv')

group_of_pspline_densities = function(X, Xmin, Xmax, degree, lambdaseq, nbasis, eval_grid_size = 500,
                                      pen_order = 3){
  
  # Fit P-spline densities by using Poisson model for bin counts
  # (e.g. "Flexible Smoothing with B-splines and Penalties, Eilers & Marx, 1996, Section 8)
  # Given a group of i.i.d. samples, this function fits P-spline densities over a range
  # of mutual smoothing parameters, then picks the parameter which minimizes the sum of AIC
  # scores across samples
  
  # INPUTS:
  # X: a list of samples (i.e. i^th entry is an i.i.d. sample of size n_i)
  # Xmin (resp. Xmax): the assumed minimum (resp. maximum )of the domain on which the densities are defined
  # degree: the degree to use for the P-spline basis (e.g. degree = 3 corresponds to cubic splines)
  # lambdaseq: a sequence of smoothing parameter values
  # nbasis: number of basis functions
  # eval_grid_size: the size of the grid on which we wish to return density values, minus 1
  # pen_order: the order of roughness penalty to apply. pen_order = r corresponds to penalization
  # of r^th-order differences between coefficients
  
  # OUTPUTS:
  # denslist: a list of the same length as X, each entry of which is the GAM
  # object corresponding to the density estimate for one of the samples in X
  # pvals: a length(X) * (eval_grid_size + 1) matrix. The i^th row contains
  # the values of the i^th density on the evaluation grid
  # xval: the grid used for evaluation, a vector of size eval_grid_size + 1
  # lambda_opt: the AIC_optimal smoothing parameter
  # AICs: a length(lambda_seq) * length(X) matrix of AIC scores, for each
  # density and smoothing parameter
  
  denslist_temp = vector('list', length(X))
  xval = seq(Xmin, Xmax, length.out = eval_grid_size + 1)
  pvals = matrix(0, length(X), eval_grid_size)
  AICs = matrix(0, length(lambdaseq), length(X))
  
  best_AIC = Inf
  
  for(j in seq_along(lambdaseq)){
    for(i in seq_along(denslist)){
      binned = hist(X[[i]], xval, plot = FALSE)
      denslist_temp[[i]] = gam(binned$counts ~ s(binned$mids, k = nbasis, bs = 'ps',
                                            m = c(degree - 1, pen_order), sp = lambdaseq[j]),
                          family = 'poisson', knots = list(`binned$mids` = xval))
      AICs[j,i] = AIC(denslist[[i]])
      print(i)
    }
    
    # Evaluation metric is sum of AIC scores for all densities
    current_AIC = sum(AICs[j,]) 
    if(current_AIC < best_AIC){
      best_AIC = current_AIC
      lambda_opt = lambdaseq[j]
      denslist = denslist_temp
    }
  }

  # Calculate pvals once we have AIC-optimal density estimates
  for(i in seq_along(denslist)){
    pvals[i,] = exp(predict.gam(denslist[[i]], type = 'lpmatrix') %*% 
                      denslist[[i]]$coefficients)/(length(X[[i]])*diff(xval)[1])
  }
  
  return(list(denslist = denslist, pvals = pvals, xval = xval,
              lambda_opt = lambda_opt, AICs = AICs))
}

plot_stepwise_functions = function(fvals, xvals, ...){
  # A function to plot multiple stepwise functions
  
  # INPUTS:
  # fvals: a matrix with nrow = number of functions and ncol = number of "bins"
  # containing the values of the stepwise functions
  # xvals: a vector of length (number of bins) + 1 containing the bin boundaries
  
  matplot(c(xvals[1], rep(xvals[2:(length(xvals)-1)], each = 2), xvals[length(xvals)]),
       apply(fvals, 1, rep, each = 2), type = 'l', ...)
}

make_RW_mat = function(r, K){
  # Ouputs a matrix RW_mat such that, if X is a (K-1) vector of independent Normals,
  # RW_mat * X is a  Gaussian random walk of order r (with first step omitted, since it
  # equals 0 for identifiability). MUCH faster than defining RW's in loops 
  
  tmp_vec = 1:r
  sign_vec = (-1)^(r - tmp_vec)
  binom_vec = choose(r, tmp_vec - 1)
  
  coef_mat = diag(sign_vec * binom_vec)
  RW_mat_full = diag(rep(1, K))
  for(i in (r+1):K){
    RW_mat_full[i, 1:(i-1)] = colSums(coef_mat %*% RW_mat_full[(i-r):(i-1), 1:(i-1)])
  }
  
  # Remove first row and column
  RW_mat = RW_mat_full[2:K, 2:K]
  return(RW_mat)
}

standardize = function(x, a, b, z, w){
  # Take a vector x with values in [a, b], and shift&scale it so its values are
  # in [z, w]
  return(((w*a - z*b) + (z-w)*x)/(a-b))
}

frodo = function(Y, X, r, K, orig_bounds = range(X), scaled_bounds = (orig_bounds - mean(unlist(X)))/sd(unlist(X)),
                 delta = rep(0.1, length(Y)), ...){
  # The main FRODO function. Sets everything up, then initalizes & runs NUTS chains
  # for the model. Currently only supports a single multilevel covariate and Gaussian responses
  
  # INPUTS:
  # Y: a vector of group-level responses
  # X: a list of the same length as Y, each entry of which contains the individual-level
  # covariate responses for a group
  # r: the order of the random walk prior on the densities
  # K: the dimensionality of the densities/regression function (i.e. number of bins)
  # orig_bounds: the assumed domain of X on the original scale (denoted [a', b'] in the manuscript)
  # scaled_bounds: the assumed domain of X after standardizing (denoted [a, b] in the manuscript)
  # delta: a vector of same length as Y containing scale factors for priors on density smoothing parameters
  # ...: additional arguments for Stan sampler
  
  # OUTPUTS:
  # stan_runs: the Stan object for the model after running sampler
  # data_list: the data fed to the Stan program
  
  N = length(Y)
  binwidth = diff(scaled_bounds)/K
  dens_RW_mat = make_RW_mat(r, K)
  beta_RW_mat = make_RW_mat(2, K)
  
  # Standardize covariates to [-1, 1]
  X_stand = lapply(X, standardize, orig_bounds[1], orig_bounds[2], scaled_bounds[1], scaled_bounds[2])
  
  # Create function to use for initialization
  # (Unfortunately, nested function definition is best way to do this)
  init_fun = function(){
    # Preliminary density estimates to initialize f_i's. Have to be careful initializing,
    # because RW prior means default initializations usually result in underflow
    smooth_amt = 10^(runif(1, -1, 2.5))
    X_densities = group_of_pspline_densities(X_stand, scaled_bounds[1], scaled_bounds[2], 0, smooth_amt, K, K, pen_order = r)
    pvals = X_densities$pvals
    # Remove any zeros so we don't have problems when we log-transform
    pvals = t(apply(pvals, 1, function(x) ifelse(x == 0, min(x[x != 0]), x)))
    # Take logs, then ensure that first value is == 0 for each density
    log_pvals = t(log(pvals[,-1]) - log(pvals[,1]))

    # Initializing theta_raw this way ensures that inital values for densities
    # aren't debilitatingly far from initial estimates, but still reasonably diffuse
    tau_dens_raw = rgamma(N, 4, 4)/(sqrt(smooth_amt)*delta)
    tau_dens = delta*tau_dens_raw
    theta_raw_temp = solve(dens_RW_mat)%*%(log_pvals + rnorm(N*(K-1), 0, 1))
    theta_raw = rbind(matrix(rnorm((r-1)*N, 0, 0.5), r-1, N), theta_raw_temp[(r):(K-1),]) %*% diag(1/tau_dens)

    tau_beta_raw = runif(1, 0.05, 2)
    
    # Initialize the parameters involved in random walk strucutres, as these
    # are the ones that create immediate problems w/default initialization
    init_list = list(
      tau_beta_raw = tau_beta_raw,
      tau_dens_raw = tau_dens_raw,
      theta_raw = theta_raw,
      beta_raw = binwidth*as.vector(solve(cbind(beta_RW_mat[,1], tau_beta_raw*beta_RW_mat[,2:(K-1)]))%*% rnorm(K-1, 0, 1))
    )
    
    # Initialize hyperparameters for "free parameters" for third- or second-order density priors
    if(r == 3) {
      init_list$sigma_X = array(1/sqrt(rgamma(1, 0.5, 0.5*mean(sapply(X_stand,var)))), dim = 1)
      xi_hat = sapply(X_stand, mean)
      init_list$xi_raw = matrix((xi_hat-mean(xi_hat))/sd(xi_hat) + rnorm(N, 0, 0.5), 1, N)
      init_list$sigma_xi = array(1/sqrt(rgamma(1, 0.5, 0.5*var(xi_hat))), dim = 1)
    } else if(r == 2){
      lambda_hat = 1/sapply(X_stand, sd)
      init_list$lambda_raw = matrix(1/rgamma(1, 1, sd(lambda_hat)/lambda_hat), 1, N)
    }
    return(init_list)
  }
  
  # Data to use in Stan program
  xval = seq(scaled_bounds[1], scaled_bounds[2], length.out = K+1) # Bin boundaries
  m = sapply(X_stand, function(x) hist(x, xval, plot = FALSE)$counts) # Counts of data in each bin, for each group
  # Explicit indication of nonzero entries of m improves efficiency and stability
  non_zero_m_inds = which(m != 0)
  non_zero_nb = length(non_zero_m_inds)
  
  data_list =  list(K = K, m = m, r = r, N = N, delta = delta, a = scaled_bounds[1], b = scaled_bounds[2],
                    cent_dens = rowSums(m)/sum(m), # "Empirical central density"
                    Y = (Y-mean(Y))/sd(Y), # Standardized regression responses
                    mu_xi_mean = (scaled_bounds[2]*orig_bounds[1] - scaled_bounds[1]*orig_bounds[2])/
                      (orig_bounds[1] - orig_bounds[2]), # Global covariate hypermean for r = 3
                    non_zero_m_inds = non_zero_m_inds, non_zero_nb = non_zero_nb,
                    eval_points = (xval[-length(xval)] + diff(xval)/2) # Bin midpoints
                    )
  
  # Run Stan program
  stan_runs = stan(file = 'frodo.stan', data = data_list, init = init_fun, ...)
  return(list(stan_runs = stan_runs, data_list = data_list))
}

frodo_with_scalar = function(Y, X, r, K, Z, orig_bounds = range(X), scaled_bounds = (orig_bounds - mean(unlist(X)))/sd(unlist(X)),
                             delta = rep(0.1, length(Y)), ...){
  # Mostly the same as the FRODO function, but with an additional (noise-free)
  # scalar group-level covariate Z
  
  # INPUTS:
  # Y: a vector of group-level responses
  # Z: a vector of the same length as Y, containing values for scalar group-level covariate
  # X: a list of the same length as Y, each entry of which contains the individual-level
  # covariate responses for a group
  # r: the order of the random walk prior on the densities
  # K: the dimensionality of the densities/regression function (i.e. number of bins)
  # orig_bounds: the assumed domain of X on the original scale (denoted [a', b'] in the manuscript)
  # scaled_bounds: the assumed domain of X after standardizing (denoted [a, b] in the manuscript)
  # delta: a vector of same length as Y containing scale factors for priors on density smoothing parameters
  # ...: additional arguments for Stan sampler
  
  # OUTPUTS:
  # stan_runs: the Stan object for the model after running sampler
  # data_list: the data fed to the Stan program
  
  N = length(Y)
  binwidth = diff(scaled_bounds)/K
  dens_RW_mat = make_RW_mat(r, K)
  beta_RW_mat = make_RW_mat(2, K)
  
  # Standardize covariates to [-1, 1]
  X_stand = lapply(X, standardize, orig_bounds[1], orig_bounds[2], scaled_bounds[1], scaled_bounds[2])
  
  # Create function to use for initialization
  # (Unfortunately, nested function definition is best way to do this)
  init_fun = function(){
    # Preliminary density estimates to initialize f_i's. Have to be careful initializing,
    # because RW prior means default initializations usually result in underflow
    smooth_amt = 10^(runif(1, -1, 2.5))
    X_densities = group_of_pspline_densities(X_stand, scaled_bounds[1], scaled_bounds[2], 0, smooth_amt, K, K, pen_order = r)
    pvals = X_densities$pvals
    # Remove any zeros so we don't have problems when we log-transform
    pvals = t(apply(pvals, 1, function(x) ifelse(x == 0, min(x[x != 0]), x)))
    # Take logs, then ensure that first value is == 0 for each density
    log_pvals = t(log(pvals[,-1]) - log(pvals[,1]))
    
    # Initializing theta_raw this way ensures that inital values for densities
    # aren't debilitatingly far from initial estimates, but still reasonably diffuse
    tau_dens_raw = rgamma(N, 4, 4)/(sqrt(smooth_amt)*delta)
    tau_dens = delta*tau_dens_raw
    theta_raw_temp = solve(dens_RW_mat)%*%(log_pvals + rnorm(N*(K-1), 0, 1))
    theta_raw = rbind(matrix(rnorm((r-1)*N, 0, 0.5), r-1, N), theta_raw_temp[(r):(K-1),]) %*% diag(1/tau_dens)
    
    tau_beta_raw = runif(1, 0.05, 2)
    
    # Initialize the parameters involved in random walk strucutres, as these
    # are the ones that create immediate problems w/default initialization
    init_list = list(
      tau_beta_raw = tau_beta_raw,
      tau_dens_raw = tau_dens_raw,
      theta_raw = theta_raw,
      beta_raw = binwidth*as.vector(solve(cbind(beta_RW_mat[,1], tau_beta_raw*beta_RW_mat[,2:(K-1)]))%*% rnorm(K-1, 0, 1))
    )
    
    # Initialize hyperparameters for "free parameters" for third- or second-order density priors
    if(r == 3) {
      init_list$sigma_X = array(1/sqrt(rgamma(1, 0.5, 0.5*mean(sapply(X_stand,var)))), dim = 1)
      xi_hat = sapply(X_stand, mean)
      init_list$xi_raw = matrix((xi_hat-mean(xi_hat))/sd(xi_hat) + rnorm(N, 0, 0.5), 1, N)
      init_list$sigma_xi = array(1/sqrt(rgamma(1, 0.5, 0.5*var(xi_hat))), dim = 1)
    } else if(r == 2){
      lambda_hat = 1/sapply(X_stand, sd)
      init_list$lambda_raw = matrix(1/rgamma(1, 1, sd(lambda_hat)/lambda_hat), 1, N)
    }
    return(init_list)
  }
  
  # Data to use in Stan program
  xval = seq(scaled_bounds[1], scaled_bounds[2], length.out = K+1) # Bin boundaries
  m = sapply(X_stand, function(x) hist(x, xval, plot = FALSE)$counts) # Counts of data in each bin, for each group
  # Explicit indication of nonzero entries of m improves efficiency and stability
  non_zero_m_inds = which(m != 0)
  non_zero_nb = length(non_zero_m_inds)
  
  data_list =  list(K = K, m = m, r = r, N = N, Z = Z, delta = delta, a = scaled_bounds[1], b = scaled_bounds[2],
                    cent_dens = rowSums(m)/sum(m), # "Empirical central density"
                    Y = (Y-mean(Y))/sd(Y), # Standardized regression responses
                    mu_xi_mean = (scaled_bounds[2]*orig_bounds[1] - scaled_bounds[1]*orig_bounds[2])/
                      (orig_bounds[1] - orig_bounds[2]), # Global covariate hypermean for r = 3
                    non_zero_m_inds = non_zero_m_inds, non_zero_nb = non_zero_nb,
                    eval_points = (xval[-length(xval)] + diff(xval)/2) # Bin midpoints
  )
  
  # Run Stan program
  stan_runs = stan(file = 'frodo_with_scalar.stan', data = data_list, init = init_fun, ...)
  return(list(stan_runs = stan_runs, data_list = data_list))
}
