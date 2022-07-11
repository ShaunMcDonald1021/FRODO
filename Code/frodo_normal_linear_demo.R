# Simulation study from Section 4.1 of manuscript
# (Additive Gaussian covariate structure, linear regression model)

# Set stuff up for plotting later
library(extrafont)
loadfonts(device = "win")
par(mfrow = c(1,2), cex.axis = 1.25, cex.lab = 1.25,
    cex.main = 1.25, mgp = c(2.1, 1, 0), family = 'LM Roman 10')
source('frodo_functions.R')

# Generate data
set.seed(24601)
N = 275 # Number of groups
n = rep(20, N) # Individuals per group
xi = rnorm(N, 0, 2)
Y = 0.3 + .4*xi + rnorm(N, 0, 0.5) # Group-level responses
X = lapply(1:N, function(i) rnorm(n[i], xi[i], 3)) # Covariate measurements

# Standardize data
X_stand = lapply(X, function(x) (x-mean(unlist(X)))/sd(unlist(X)))
Y_stand = (Y-mean(Y))/sd(Y)

# Assumed domain for covariate densities on original and standardized scales
orig_bounds = c(sum(range(X)) - 11, 11)
scale_bounds = (orig_bounds - mean(unlist(X)))/sd(unlist(X))

K = 10
r = 3

# Run FRODO
frodo_runs = frodo(Y, X, r, K, orig_bounds = orig_bounds, scale_bounds = scale_bounds, delta = rep(0.1, N),
                   verbose = TRUE, refresh = 2, iter = 1750, warmup = 750,
                   control = list(max_treedepth = 12, adapt_delta = 0.985), cores= 4, chains = 4)
stan_runs = frodo_runs$stan_runs

# Check that sampler converged & behaved well, we got enough effective samples, etc.
check_hmc_diagnostics(stan_runs)
summ = summary(stan_runs)$summary
n_eff = summ[,'n_eff']
head(sort(n_eff))
Rhat=  summ[,'Rhat']
tail(sort(Rhat))

# Run naive scalar model
naive_scalar_stan_runs = stan('naive_scalar_linmod.stan', data = list(Y = Y_stand, X = sapply(X_stand, mean), N = N),
                              verbose = TRUE, refresh = 2, iter = 1750, warmup = 750, cores = 4, chains = 4)

# Run hierarchical scalar model
scalar_stan_runs = stan('hierarchical_scalar_normal_linear.stan',
                        data = list(Y = Y_stand, X = unlist(X_stand), N = N, M = length(unlist(X)), n = n),
                        verbose = TRUE, refresh = 2, iter = 1750, warmup = 750, cores = 4, chains = 4)

params = extract(stan_runs)
scalar_params = extract(scalar_stan_runs)
naive_scalar_params = extract(naive_scalar_stan_runs)
# save.image('normal_linear')

# Plot regression function (on original scale)
xval = seq(orig_bounds[1], orig_bounds[2], length.out = K+1)
evalseq = seq(orig_bounds[1], orig_bounds[2], length.out = 1001)
xval_plotseq = c(xval[1], rep(xval[2:K], each = 2), xval[K+1])

plot_stepwise_density(rbind(apply(params$beta,2,mean),apply(params$beta,2,quantile,0.025),
                            apply(params$beta,2,quantile,0.975))*sd(Y), xval,
                      main = 'Regression function', xlab = expression(italic(x)),
                      ylab = expression(italic(beta(x))), lty = c(1, 3, 3), col = rep(1,3),
                      lwd = c(2, 0.75, 0.75))
upper_CI = rep(apply(params$beta,2,quantile,0.975)*sd(Y), each = 2)
lower_CI = rep(apply(params$beta,2,quantile,0.025)*sd(Y), each = 2)
polygon(c(xval_plotseq, rev(xval_plotseq)), c(lower_CI, rev(upper_CI)),
        col = rgb(0.25,0.25,0.25, 0.25), border = NA)

# True regression function (shifted to have integral 0)
abline(-frodo_runs$data_list$cent_dens %*% diff(0.4*xval^2/2)*K/diff(orig_bounds),
       0.4, lty = 4, lwd = 2)

# Regression function from naive scalar model (shifted to have integral 0)
abline(-frodo_runs$data_list$cent_dens %*% (sd(Y)*mean(naive_scalar_params$beta)*diff(xval^2)*K/(2*diff(orig_bounds)*sd(unlist(X)))),
       mean(naive_scalar_params$beta)*sd(Y)/sd(unlist(X)), col = 'red', lwd = 1.75)

# Regression function from hierarchical scalar model (shifted to have integral 0)
abline(-frodo_runs$data_list$cent_dens %*% (sd(Y)*mean(scalar_params$beta)*diff(xval^2)*K/(2*diff(orig_bounds)*sd(unlist(X)))),
       mean(scalar_params$beta)*sd(Y)/sd(unlist(X)), col = 'blue', lwd = 1.75)

legend('topleft', legend = c('FRODO posterior mean', 'FRODO pointwise 95% C.I.', 'True function', 
                             "Scalar posterior mean", "Naive scalar posterior mean"), 
       lty = c(1, 0, 4, 1, 1), lwd = c(2, 0, 2, 1.75, 1.75), pch = c(NA,22, NA, NA, NA), 
       col = c('black',rgb(0.25, 0.25, 0.25, 0.25), 'black', 'blue', 'red'), 
       pt.bg = c(NA,rgb(0.25, 0.25, 0.25, 0.25), NA, NA), pt.cex = 2)

# Plot estimated vs. true responses
mean_Yhat = sd(Y)*apply(params$fitted_values, 2, mean) + mean(Y)
# 95% credible intervals for mean responses (NOT prediction intervals)
lower_CI_Yhat = sd(Y)*apply(params$fitted_values, 2, quantile, 0.025) + mean(Y)
upper_CI_Yhat = sd(Y)*apply(params$fitted_values, 2, quantile, 0.975) + mean(Y)
plot(mean_Yhat, Y, main = 'Response values', pch = 16, cex = 0.5,
     xlab = 'Predicted responses', ylab = 'Actual responses',
     ylim = c(min(lower_CI_Yhat), max(upper_CI_Yhat)))
abline(0, 1, lwd = 1.5)
polygon(c(sort(mean_Yhat), rev(sort(mean_Yhat))),
        c(lower_CI_Yhat[order(mean_Yhat)], rev(upper_CI_Yhat[order(mean_Yhat)])),
        col = rgb(0.25,0.25,0.25, 0.25), border = NA)
legend('topleft', legend = c('Posterior mean responses', "95% C.I.'s"), 
       pch = c(16,22), 
       col = c('black',rgb(0.25, 0.25, 0.25, 0.25)), 
       pt.bg = c(NA,rgb(0.25,0.25,0.25,0.25)), pt.cex = 2)

# Plot estimated densities alongside true densities
par(mfrow = c(1,3))
min_xi = which.min(xi)
max_xi = which.max(xi)
mean_xi = which.min(abs(xi - mean(xi)))

# Group with lowest true mean
plot_stepwise_density(rbind(apply(params$phi[,,min_xi],2,mean),apply(params$phi[,,min_xi],2,quantile,0.025),
                            apply(params$phi[,,min_xi],2,quantile,0.975))*K/diff(orig_bounds), xval,
                      main = bquote(italic(f)[.(min_xi)]), lab = expression(italic(x)),
                      ylab = expression(italic(f(x))), lty = c(1, 3, 3), col = rep(1,3),
                      lwd = c(2, 0.75, 0.75), ylim = c(0, 0.18))
polygon(c(xval_plotseq, rev(xval_plotseq)),
        rep(c(apply(params$phi[,,min_xi],2,quantile,0.025), rev(apply(params$phi[,,min_xi],2,quantile,0.975))),
            each=2)*K/diff(orig_bounds), 
        col = rgb(0.25,0.25,0.25, 0.25), border = NA)
lines(evalseq, dnorm(evalseq, xi[min_xi], 3), col = 'red', lwd = 2)
legend('topright', legend = c('FRODO posterior mean', 'FRODO pointXise 95% C.I.', 'True density'), 
       lty = c(1, 0, 1), lwd = c(2, 0, 2), pch = c(NA,22, NA), 
       col = c('black',rgb(0.25, 0.25, 0.25, 0.25), 'red'), 
       pt.bg = c(NA,rgb(0.25, 0.25, 0.25, 0.25), NA), pt.cex = 2)

# Group with true mean closest to mean(xi)
plot_stepwise_density(rbind(apply(params$phi[,,mean_xi],2,mean),apply(params$phi[,,mean_xi],2,quantile,0.025),
                            apply(params$phi[,,mean_xi],2,quantile,0.975))*K/diff(orig_bounds), xval,
                      main = bquote(italic(f)[.(mean_xi)]), xlab = expression(italic(x)),
                      ylab = '', lty = c(1, 3, 3), col = rep(1,3), lwd = c(2, 0.75, 0.75), ylim = c(0, 0.18))
polygon(c(xval_plotseq, rev(xval_plotseq)),
        rep(c(apply(params$phi[,,mean_xi],2,quantile,0.025),
              rev(apply(params$phi[,,mean_xi],2,quantile,0.975))),each=2)*K/diff(orig_bounds), 
        col = rgb(0.25,0.25,0.25, 0.25), border = NA)
lines(evalseq, dnorm(evalseq, xi[mean_xi], 3), col = 'red', lwd = 2)

# Group with highest true mean
plot_stepwise_density(rbind(apply(params$phi[,,max_xi],2,mean),apply(params$phi[,,max_xi],2,quantile,0.025),
                            apply(params$phi[,,max_xi],2,quantile,0.975))*K/diff(orig_bounds), xval,
                      main = bquote(italic(f)[.(max_xi)]), xlab = expression(italic(x)),
                      ylab = '', lty = c(1, 3, 3), col = rep(1,3), lwd = c(2, 0.75, 0.75), ylim = c(0, 0.18))
polygon(c(xval_plotseq, rev(xval_plotseq)),
        rep(c(apply(params$phi[,,max_xi],2,quantile,0.025),
              rev(apply(params$phi[,,max_xi],2,quantile,0.975))),each=2)*K/diff(orig_bounds), 
        col = rgb(0.25,0.25,0.25, 0.25), border = NA)
lines(evalseq, dnorm(evalseq, xi[max_xi], 3), col = 'red', lwd = 2)
