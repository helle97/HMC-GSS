library(microbenchmark)
library(doParallel)
library(phytools)
library(ape)
library(numDeriv)
library(hmclearn)
library(foreach)
library(mvMORPH)

# we need to compute the reference distribution
# for each parameter based on the MCMC samples.
# we estimate the reference distribution using 
# density in R.

#for cleaning up parallel computing
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()


## Simulated dataset
set.seed(14)
# Generating a random tree with 50 species
tree<-pbtree(n=50)
X <- mvSIM(tree, nsim = 1, error = NULL, model = c("BM1"),
           param = list(theta = c(4,6),
                        sigma = matrix(c(0.4, 0.1, 0.1, 0.6),nrow = 2), 
                        trend = c(-0.2, 0.2)))

dim(X)

params <- c(rnorm(1,5,0.1), rnorm(1,5, 0.1), 0.2, 0.2,0.2, 0.2 )



LogLik_simple_posterior <- function(Y, tree, theta){
  n <- nrow(Y)
  p <- ncol(Y)
  Y <- c(Y) # converting to n*p length vector
  sig2mu1 <- 0.5;  sig2mu2 <- 0.5
  mu1_mean <- 5; mu2_mean <- 5
  # gamma params
  gshape = 35
  gscale = 1/150
  
  
  D <- matrix(rep(NA, n*p*p), byrow = FALSE, ncol = p, nrow =n*p)
  for (i in 1:(n*p)){
    for (j in 1:p){
      if (((j-1)*n < i) & (i <= j*n)){D[i,j] <- 1}
      else {D[i,j] <- 0}
    }
  }
  
  C <- ape::cophenetic.phylo(tree)
  tip.heights <- ape::node.depth.edgelength(tree)[1:n]
  diag(C) <- tip.heights
  
  # defining parameters
  mu1 <- theta[1]; mu2 <- theta[2]
  mu <- c(mu1, mu2)
  gamma11 <- theta[3];  gamma22 <- theta[4]; # gamma21 <- theta[4]; gamma12 <- theta[4];
  gamma <- matrix(c(gamma11, 0, 0, gamma22), nrow = 2)
  R11 <- theta[5]; R22 <- theta[6]; # R21 <- theta[7]; R12 <- theta[7]; R22
  R <- matrix(c(R11, 0, 0, R22), 2)

  VCOV <- as.matrix(kronecker((gamma + R), C))
  
  model_loglik <- -0.5*(t(Y-D%*%mu)%*%(solve(VCOV))%*%(Y-D%*%mu)) - determinant(VCOV, logarithm = TRUE)$modulus[[1]]/2  - (n*p/2)*log(2*pi) 
  
  mu_loglik <- dnorm(mu1, mean = mu1_mean, sd = sig2mu1, log = TRUE) + dnorm(mu2, mean = mu2_mean, sd = sig2mu2, log = TRUE)  #- 1/2 * mu1^2/(2*sig2mu) - 1/2* mu2^2/(2*sig2mu)
  
  # variances gamma(2,2) prior, covariance beta(2,5) prior
  gamma_loglik <-  dgamma(gamma11, shape =gshape, scale = gscale, log = TRUE) + dgamma(gamma22, shape =gshape, scale = gscale, log = TRUE) #+ 2*dbeta(gamma21, shape1 =2, shape2 = 5, log = TRUE ) #-2*exp(-(gamma11)) -2*exp(- (gamma22))+ (log(30)*2*gamma21 -8*log(1+exp(gamma21))) 
  
  # variances gamma(2,2) prior, covariance beta(2,5) prior
  R_loglik <- dgamma(R11, shape =gshape, scale = gscale, log = TRUE) + dgamma(R22, shape =gshape, scale = gscale, log = TRUE) #+ 2*dbeta(R21, shape1 =2, shape2 = 5, log = TRUE) #-2*exp(- (R11)) -2*exp(-(R11)) + ((log(30)*2*R21 -8*log(1+exp(R21)))) 
  
  return(model_loglik + mu_loglik + gamma_loglik + R_loglik)
}


dim(X) #needs to be on form species x traits
LogLik_simple_posterior(Y = X, tree = tree, theta = params)



g_LogLik_simple_posterior <- function(Y, tree, theta){
  list_side <- rep(NA, 6)
  # check limits for variances, param must be over 0
  if (theta[3]<=0.01){ list_side[3] <- +1}
  if (theta[4]<=0.01){ list_side[4] <- +1}
  if (theta[5]<=0.01){ list_side[5] <- +1}
  if (theta[6]<=0.01){ list_side[6] <- +1}
  return(numDeriv::grad(func = LogLik_simple_posterior, x = theta, method = "Richardson", method.args = list(side= list_side), tree = tree, Y = Y))
}
g_LogLik_simple_posterior(Y = X, tree =  tree, theta = params)


L <- 15
eps_vals <- c(0.01, 0.01, rep(0.002, 4))
length(eps_vals)
N <- 500
theta <- params
set.seed(143)
#})

system.time(hmc_sample_model1 <- hmclearn::hmc(N, theta.init = theta,
                                                       epsilon = eps_vals, L = L,
                                                       logPOSTERIOR = LogLik_simple_posterior,
                                                       glogPOSTERIOR = g_LogLik_simple_posterior,
                                                       varnames = c("mu1", "mu2", "gamma11", "gamma22", "R11", "R22"),
                                                       param = list(Y = Y, tree = tree), chains = 2))



summary(hmc_sample_model1, burnin = 150)
plot(hmc_sample_model1, burnin = 150)
model1_accept_hmc <- mean(hmc_sample_model1$accept/N)
model1_hmc_rhat <- summary(hmc_sample_model1, burnin = 150)[,8]
model1_hmc_rhat
diagplots(hmc_sample_model1, burnin = 150, comparison.theta = c(4,6,NA,NA,NA,NA))


# #library(profvis)
# bench <- microbenchmark(hmclearn::hmc(N, theta.init = theta,
#                                       epsilon = eps_vals, L = L,
#                                       logPOSTERIOR = LogLik_simple_posterior,
#                                       glogPOSTERIOR = g_LogLik_simple_posterior,
#                                       varnames = c("mu1", "mu2", "gamma11", "gamma22", "R11", "R22"),
#                                       param = list(Y = Y, tree = tree), chains = 1),
#                         hmclearn::hmc(N, theta.init = theta,
#                                       epsilon = eps_vals, L = L,
#                                       logPOSTERIOR = LogLik_simple_posterior,
#                                       glogPOSTERIOR = g_LogLik_simple_posterior_simple,
#                                       varnames = c("mu1", "mu2", "gamma11", "gamma22", "R11", "R22"),
#                                       param = list(Y = Y, tree = tree), chains = 1))
# 
# 
# 
# bench2 <- microbenchmark(g_LogLik_simple_posterior_simple(Y = X, tree =  tree, theta = params),
#                          g_LogLik_simple_posterior(Y = X, tree =  tree, theta = params))
# 
# system.time(fit1 <- hmclearn::hmc(N, theta.init = theta,
#                      epsilon = eps_vals, L = L,
#                      logPOSTERIOR = LogLik_simple_posterior,
#                      glogPOSTERIOR = g_LogLik_simple_posterior,
#                      varnames = c("mu1", "mu2", "gamma11", "gamma22", "R11", "R22"),
#                      param = list(Y = Y, tree = tree), chains = 2))
# 
# 
# 
# summary(fit1)
# plot(fit1, burnin = 200)
# save.image(file='model1_sim_results.RData')


# save.image(file='model1_sim_results_simpleg.RData')



set.seed(143)



# Main function
GSS_HMC <- function(Y, tree, N, K, beta, model_log_posterior, model_gradient, initial_params, L, eps_vals){
  # Defining names for the parallel computation (does not work with global function names)
  LogLik_simple_posterior <- model_log_posterior
  g_LogLik_simple_posterior <- model_gradient
  
  #making design matrixes
  n <- nrow(Y)
  p <- ncol(Y)
  Design_matrix <- matrix(rep(NA, n*p*p), byrow = FALSE, ncol = p, nrow =n*p)
  for (i in 1:(n*p)){
    for (j in 1:p){
      if (((j-1)*n < i) & (i <= j*n)){Design_matrix[i,j] <- 1}
      else {Design_matrix[i,j] <- 0}
      }
  }
    
  C_mat <- ape::cophenetic.phylo(tree)
  tip.heights <- ape::node.depth.edgelength(tree)[1:n]
  diag(C_mat) <- tip.heights
  
  # Make vectors to save elements of estimator
  r_vec <- c()
  eta_vec <- c()
  # Activate cluster for foreach library and pass libraries
  cl <- parallel::makeCluster(detectCores()-1)
  registerDoParallel(cl)
  # parallel computation of MCMC samples
  MCMC_init <- function(theta){
    foreach(times=1:3,
            .combine = rbind,
            .export = c('model_log_posterior', 'model_gradient', 'LogLik_simple_posterior','g_LogLik_simple_posterior',
                        'N', 'L', 'eps_vals', 'Y', 'tree', 'Design_matrix')) %dopar%
      hmclearn::hmc(N, theta.init = theta,
                    epsilon = eps_vals, L = L,
                    logPOSTERIOR = model_log_posterior,
                    glogPOSTERIOR = model_gradient,
                    varnames = c("mu1", "mu2", "gamma11", "gamma22", "R11", "R22"),
                    param = list(Y = Y, tree = tree), chains = 1)$theta
  } 

  hmc_samples <- MCMC_init(params)
  stopCluster(cl)
  
  samples <- rbind.data.frame(do.call(rbind, hmc_samples[[1]])[(floor(N/3) + 1): N,], # use second two thirds of chain
                              do.call(rbind, hmc_samples[[2]])[(floor(N/3) + 1): N,], # use second two thirds of chain
                              do.call(rbind, hmc_samples[[3]])[(floor(N/3) + 1): N,]) # use second two thirds of chain
  
  # find the sample mean and sds to parametrize reference distribution
  sample_means <- sapply(samples, mean)
  sample_sd <- sapply(samples, sd)
  sample_sd <- ifelse(sample_sd <= 0.001, 0.001, sample_sd)
  
  # define new reference distribution function
  ref_dist <- function(theta, sm = sample_means, sd = sample_sd){
    # both entries of mu has normal ref dist and prior
    mu_log_lik <- (dnorm(theta[1], mean = sm[1], sd = (sd[1]), log = TRUE) 
               + dnorm(theta[2], mean = sm[2], sd = (sd[2]), log = TRUE))  #- 1/2 * mu1^2/(2*sig2mu) - 1/2* mu2^2/(2*sig2mu)
    # variances gamma(2,2) prior, covariance beta(2,5) prior
    gamma_log_lik <- (dgamma(theta[3], shape = (sm[3]^2)/sd[3]^2, scale = sd[3]^2/sm[3], log = TRUE)
                    + dgamma(theta[4], shape = (sm[4]^2)/(sd[4]^2), scale = sd[4]^2/sm[4], log = TRUE))
                   # + 2*dbeta(theta[4], shape1 = - (sm[4]*(sm[4]^2 - sm[4] + sd[4]))/sd[4], shape2 = ((sm[4] - 1)*(sm[4]^2 - sm[4] + sd[4]))/sd[4], log = TRUE)) #-2*exp(-(gamma11)) -2*exp(- (gamma22))+ (log(30)*2*gamma21 -8*log(1+exp(gamma21))) 
    # variances gamma(2,2) prior, covariance beta(2,5) prior
    R_log_lik <- (dgamma(theta[5], shape = (sm[5]^2)/sd[5]^2, scale = sd[5]^2/sm[5], log = TRUE)
                + dgamma(theta[6], shape = (sm[6]^2)/sd[6]^2, scale = sd[6]^2/sm[6], log = TRUE))
                #+ 2*dbeta(theta[7], shape1 = - (sm[7]*(sm[7]^2 - sm[7] + sd[7])/sd[7]), shape2 = ((sm[7] - 1)*(sm[7]^2 - sm[7] + sd[7]))/sd[7], log = TRUE)) #-
    return(exp(mu_log_lik + gamma_log_lik + R_log_lik))            
  }
    
  # Define log_q_beta which we will sample from in GSS
  log_q_beta <- function(Y, tree, theta, beta, sm = sample_means, sd = sample_sd,
                         model_log_posterior = LogLik_simple_posterior){
    return(beta * model_log_posterior(Y, tree, theta) + (1-beta)*log(ref_dist(theta)))
  }
  
  # define gradient for log_q_beta
  g_log_q_beta <- function(Y, tree, theta, beta){
    list_side <- rep(NA, 6)
    # check limits for variances, param must be over 0
    if (theta[3]<=0.01){ list_side[3] <- +1}
    if (theta[4]<=0.01){ list_side[4] <- +1} 
    if (theta[5]<=0.01){ list_side[5] <- +1}
    if (theta[6]<=0.01){ list_side[6] <- +1} 
    return(numDeriv::grad(func =log_q_beta, x = theta, method = "simple",
                          method.args = list(side= list_side), tree = tree, Y = Y, beta = beta))
  }
  # define function that calculates posterior divided by reference
  calc_frac <- function(Y, tree, theta, log_posterior = model_log_posterior, ref = ref_dist){
    fractions <- exp(log_posterior(Y, tree, theta) - log(ref(theta)))
    return(fractions)
  }
  cl1 <- parallel::makeCluster(detectCores()-1)
  registerDoParallel(cl1)
  for (k in (2:K)){ # We sum for beta_1 in the estimator which is the same as beta[2] 
    # calculate term for estimator
    fracs <-  apply(samples, 1, function(x) calc_frac(Y, tree, as.numeric(x)) )
                    #sapply(samples, function (theta) calc_frac(Y, tree, as.numeric(theta)))
    eta <- max(fracs)
    # calculating a single r_hat and eta
    r_vec <- c(r_vec, (1/length(fracs))*sum( exp(log(fracs)-log(eta))^(beta[k]-beta[k-1])))
    eta_vec <- c(eta_vec, (beta[k] - beta[k-1]) * log(eta))
    
    print(paste0("beta_k: ", beta[k]))

    if (k == K){
      # do nothing
    }
    else if (k != K){
      # Activate cluster for foreach library and pass libraries
      beta_temp <- beta[k]

      eps_vals <- c(0.005, 0.005, rep(0.00005,4))
      MCMC_samples <- function(theta){
        foreach(times=1:3,
                .combine = rbind,
                .export = c('log_q_beta', 'g_log_q_beta', 'LogLik_simple_posterior','g_LogLik_simple_posterior',
                            'model_log_posterior', 'N', 'L', 'eps_vals', 'Y', 'tree', 'beta_temp',
                            'ref_dist', 'sample_means', 'sample_sd')) %dopar%
          hmclearn::hmc(N, theta.init = theta,
                        epsilon = eps_vals, L = L,
                        logPOSTERIOR = log_q_beta,
                        glogPOSTERIOR = g_log_q_beta,
                        param = list(Y = Y, tree = tree, beta = beta_temp), chains = 1)$theta

      } 
      hmc_samples <- MCMC_samples(as.numeric(samples[dim(samples)[1],]))
      # stopCluster(cl1)
      samples <- rbind.data.frame(do.call(rbind, hmc_samples[[1]])[(floor(N/3) + 1): N,], # use second two thirds of chain
                                  do.call(rbind, hmc_samples[[2]])[(floor(N/3) + 1): N,], # use second two thirds of chain
                                  do.call(rbind, hmc_samples[[3]])[(floor(N/3) + 1): N,]) # use second two thirds of chain
    
    }
  }
  stopCluster(cl1)
  print(paste0("eta_vec : ", eta_vec))
  print(paste0("r_vec : ", r_vec))
  log_rhat <- sum(eta_vec) + sum(log(r_vec)) 
  return(list(log_r = log_rhat, r_vector=  r_vec, eta_vector = eta_vec))
}

L <- 15
eps_vals <- c(0.02, 0.02, rep(0.002,4))
N <- 500
theta <- params
set.seed(143)
K = 10
beta <- qbeta(seq(0, 1, length.out = K), 0.3, 1)
# mu1 mu2, gamma11, gamma21, gamma22, R11, R21, R22
initial_params <- params

unregister_dopar()
system.time(log_r <- GSS_HMC(Y = X, tree = tree, N = N, K, beta = beta, 
        model_log_posterior = LogLik_simple_posterior,
        model_gradient = g_LogLik_simple_posterior,
        initial_params = params, L = L, eps_vals = eps_vals ))

log_r

save.image(file='GSS_model1_FINAL.RData')









