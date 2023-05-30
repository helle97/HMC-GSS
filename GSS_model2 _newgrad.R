library(mvtnorm)
library(doParallel)
library(phytools)
library(ape)
library(numDeriv)
library(hmclearn)
library(mvMORPH)
library(foreach)
library(MCMCpack)
library(tidyr)
library(ggplot2)
library(bayesplot)
library(mcmcplots)

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


#Y <- t(Y)
params <- c(rnorm(1,5,0.1), rnorm(1,5,0.1), 0.3, 0.2,
            0.2, 0.2, -0.2, 0.2 )
length(params)
Y <- X

{  
n <- dim(Y)[1]
p <- dim(Y)[2]
Des_mat <- matrix(rep(NA, n*p*p), byrow = FALSE, ncol = p, nrow =n*p)
for (i in 1:(n*p)){
  for (j in 1:p){
    if (((j-1)*n < i) & (i <= j*n)){ Des_mat[i,j] <- 1}
    else {Des_mat[i,j] <- 0}
  }
}
C_mat <- ape::cophenetic.phylo(tree)
tip.heights <- ape::node.depth.edgelength(tree)[1:n]
diag(C_mat) <- tip.heights
C_mat_inv <- solve(C_mat)
}


LogLik_simple_drift_posterior <- function(Y, tree, theta, D = Des_mat, C = C_mat, C_inv = C_mat_inv){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  Y <- c(Y) # converting to n*p length vector
  # set parameters for prior distributions
  
  gshape = 35; grate = 150 #informative
  #gshape = 4; grate = 2 #weakly informative
  
  #sig2mu1 <- 5;  sig2mu2 <- 5;#weakly informative
  sig2mu1 <- 1;  sig2mu2 <- 1; #informative
  mu1_mean <- 5; mu2_mean <- 5 
  
  #sig2v1 <- 5; sig2v2 <- 5 # weakly informative
  sig2v1 <- 0.5; sig2v2 <- 0.5 # informative
  v1_mean <- -0.1; v2_mean <- 0.1
  
  # defining parameters
  mu1 <- theta[1]; mu2 <- theta[2]
  mu <- c(mu1, mu2)
  gamma11 <- theta[3];  gamma22 <- theta[4]; # gamma21 <- theta[4]; gamma12 <- theta[4];
  gamma <- matrix(c(gamma11, 0, 0, gamma22), nrow = 2)
  R11 <- theta[5]; R22 <- theta[6]; # R21 <- theta[7]; R12 <- theta[7]; R22
  R <- matrix(c(R11, 0, 0, R22), 2)
  v1 <- theta[7]; v2 <- theta[8] 
  v <- c(v1, v2) #drift param
  
  
  VCOV <- as.matrix(kronecker(solve(gamma + R), C_inv))
  
  model_loglik <- -0.5*(t(Y-(D%*%mu + c(diag(C), diag(C))*(D%*%v)))%*%(VCOV)%*%(Y-(D%*%mu + c(diag(C), diag(C))*(D%*%v)))) + determinant(VCOV, logarithm = TRUE)$modulus[[1]]/2  - (n*p/2)*log(2*pi) 
  
  mu_loglik <- dnorm(mu1, mean = mu1_mean, sd = sig2mu1, log = TRUE) + dnorm(mu2, mean = mu2_mean, sd = sig2mu2, log = TRUE)  #- 1/2 * mu1^2/(2*sig2mu) - 1/2* mu2^2/(2*sig2mu)
  
  # variances gamma(2,2) prior, covariance beta(2,5) prior
  gamma_loglik <-  dgamma(gamma11,  shape=gshape, rate=grate, log = TRUE) + dgamma(gamma22,  shape=gshape, rate=grate, log = TRUE) #+ 2*dbeta(gamma21, shape1 =2, shape2 = 5, log = TRUE ) #-2*exp(-(gamma11)) -2*exp(- (gamma22))+ (log(30)*2*gamma21 -8*log(1+exp(gamma21))) 
  
  # variances gamma(2,2) prior, covariance beta(2,5) prior
  R_loglik <- dgamma(R11,  shape=gshape, rate=grate, log = TRUE) + dgamma(R22,  shape=gshape, rate=grate, log = TRUE) #+ 2*dbeta(R21, shape1 =2, shape2 = 5, log = TRUE) #-2*exp(- (R11)) -2*exp(-(R11)) + ((log(30)*2*R21 -8*log(1+exp(R21)))) 
  
  v_loglik <- dnorm(v1, mean = v1_mean, sd = sig2v1, log = TRUE) + dnorm(mu2, mean = v1_mean, sd = sig2v1, log = TRUE)  #- 1/2 * mu1^2/(2*sig2mu) - 1/2* mu2^2/(2*sig2mu)
  
  return(model_loglik + mu_loglik + gamma_loglik + R_loglik)
}



dim(X) #needs to be on form species x traits
LogLik_simple_drift_posterior(Y = X, tree = tree, theta = params)


g_LogLik_simple_drift_posterior_simple <- function(Y, tree, theta){
  list_side <- rep(NA, 8)
  #print(cat("theta: ",theta, "\n", sep= " ,"))
  # check limits for variances, param must be over 0
  if (theta[3]<=0.09){ list_side[3] <- +1}
  if (theta[4]<=0.09){ list_side[4] <- +1} 
  if (theta[5]<=0.09){ list_side[5] <- +1}
  if (theta[6]<=0.09){ list_side[6] <- +1} 
  return(numDeriv::grad(func = LogLik_simple_drift_posterior, x = theta, method = "simple", method.args = list(side= list_side), tree = tree, Y = Y))
}
g_LogLik_simple_drift_posterior_simple(Y = X, tree =  tree, theta = params)



L <- 15
eps_vals <- c(0.05, 0.05, rep(0.0021, 4), 0.05, 0.05)
length(eps_vals)
N <- 500
theta <- params
set.seed(143)

system.time(model2_hmc_simpleg <- hmclearn::hmc(N, theta.init = theta,
                          epsilon = eps_vals, L = L,
                          logPOSTERIOR = LogLik_simple_drift_posterior,
                          glogPOSTERIOR = g_LogLik_simple_drift_posterior_simple,
                          varnames = c("mu1", "mu2", "gamma11", "gamma22", "R11", "R22", "v1", "v2"),
                          param = list(Y = Y, tree = tree), chains = 3))


summary(model2_hmc_simpleg, burnin = 150)
plot(model2_hmc_simpleg, burnin = 150)
model2_accept_hmc_new <- model2_hmc_simpleg$accept/N
model2_hmc_rhat <- summary(model2_hmc_simpleg, burnin = 150)[,8]
diagplots(model2_hmc_simpleg, burnin=150, comparison.theta=c(4,6,NA,NA,NA,NA,-1,1))
diagplots(model2_hmc_simpleg, burnin=150,plotfun = 1)
# save.image(file='model2_sim_results_simpleg.RData')



# Main function
GSS_HMC <- function(Y, tree, N, K, beta, model_log_posterior, model_gradient, initial_params, L, eps_vals){
  # Defining names for the parallel computation (does not work with global function names)
  LogLik_simple_drift_posterior <- model_log_posterior
  g_LogLik_simple_drift_posterior_simple <- model_gradient
  
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  Des_mat <- matrix(rep(NA, n*p*p), byrow = FALSE, ncol = p, nrow =n*p)
  for (i in 1:(n*p)){
    for (j in 1:p){
      if (((j-1)*n < i) & (i <= j*n)){ Des_mat[i,j] <- 1}
      else {Des_mat[i,j] <- 0}
    }
  }
  C_mat <- ape::cophenetic.phylo(tree)
  tip.heights <- ape::node.depth.edgelength(tree)[1:n]
  diag(C_mat) <- tip.heights
  C_mat_inv <- solve(C_mat)
  
  
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
            .export = c('model_log_posterior', 'model_gradient', 'LogLik_simple_drift_posterior','g_LogLik_simple_drift_posterior_simple',
                        'N', 'L', 'eps_vals', 'Y', 'tree', 'Des_mat', 'C_mat', 'C_mat_inv')) %dopar%
      hmclearn::hmc(N, theta.init = theta,
                    epsilon = eps_vals, L = L,
                    logPOSTERIOR = model_log_posterior,
                    glogPOSTERIOR = model_gradient,
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
    # variances gamma prior
    gamma_log_lik <- (dgamma(theta[3], shape = (sm[3]^2)/sd[3]^2, scale = sd[3]^2/sm[3], log = TRUE)
                      + dgamma(theta[4], shape = (sm[4]^2)/(sd[4]^2), scale = sd[4]^2/sm[4], log = TRUE))
     # variances gamma prior
    R_log_lik <- (dgamma(theta[5], shape = (sm[5]^2)/sd[5]^2, scale = sd[5]^2/sm[5], log = TRUE)
                  + dgamma(theta[6], shape = (sm[6]^2)/sd[6]^2, scale = sd[6]^2/sm[6], log = TRUE))
    # normal prior
    v_log_lik <- (dnorm(theta[7], mean = sm[7], sd = (sd[7]), log = TRUE) 
                   + dnorm(theta[8], mean = sm[8], sd = (sd[8]), log = TRUE))  #- 1/2 * mu1^2/(2*sig2mu) - 1/2* mu2^2/(2*sig2mu)
    
    return(exp(mu_log_lik + gamma_log_lik + R_log_lik + v_log_lik))            
  }
  
  # Define log_q_beta which we will sample from in GSS
  log_q_beta <- function(Y, tree, theta, beta, sm = sample_means, sd = sample_sd,
                         model_log_posterior = LogLik_simple_drift_posterior, C = C_mat, D = Des_mat, C_inv = C_mat_inv){
    return(beta * model_log_posterior(Y, tree, theta, D, C, C_inv) + (1-beta)*log(ref_dist(theta)))
  }
  
  # define gradient for log_q_beta
  g_log_q_beta <- function(Y, tree, theta, beta){
    list_side <- rep(NA, 8)
    # check limits for variances, param must be over 0
    if (theta[3]<=0.09){ list_side[3] <- +1}
    if (theta[4]<=0.09){ list_side[4] <- +1} 
    if (theta[5]<=0.09){ list_side[5] <- +1}
    if (theta[6]<=0.09){ list_side[6] <- +1} 
    return(numDeriv::grad(func =log_q_beta, x = theta, method = "simple",
                          method.args = list(side= list_side), tree = tree, Y = Y, beta = beta))
  }
  # define function that calculates posterior divided by reference
  calc_frac <- function(Y, tree, theta, log_posterior = model_log_posterior, ref = ref_dist, D = Des_mat, C = C_mat,  C_inv = C_mat_inv){
    fractions <- exp(log_posterior(Y, tree, theta, D, C, C_inv) - log(ref(theta)))
    return(fractions)
  }
  
  for (k in (2:K)){ # We sum for beta_1 in the estimator which is the same as beta[2] 
    # calculate term for estimator
    fracs <-  apply(samples, 1, function(x) calc_frac(Y, tree, as.numeric(x)) )
    #sapply(samples, function (theta) calc_frac(Y, tree, as.numeric(theta)))
    eta <- max(fracs)
    # calculating a single r_hat and eta
    r_vec <- c(r_vec, (1/length(fracs))*sum( exp(log(fracs)-log(eta))^(beta[k]-beta[k-1])))
    eta_vec <- c(eta_vec, (beta[k] - beta[k-1]) * log(eta))
    
    print(paste0("beta_k: ", beta[k]))
    # simulate more sample parameters
    cl1 <- parallel::makeCluster(detectCores()-1)
    registerDoParallel(cl1)
    if (k == K){
      # do nothing
    }
    else if (k != K){
      # Activate cluster for foreach library and pass libraries
      beta_temp <- beta[k]
      MCMC_samples <- function(theta){
        foreach(times=1:3,
                        .combine = rbind,
                        .export = c('log_q_beta', 'g_log_q_beta', 'LogLik_simple_drift_posterior','g_LogLik_simple_drift_posterior_simple',
                                    'model_log_posterior', 'N', 'L', 'eps_vals', 'Y', 'tree', 'beta_temp',
                                    'ref_dist', 'sample_means', 'sample_sd', 'Des_mat', 'C_mat', 'C_mat_inv')) %dopar%
          hmclearn::hmc(N, theta.init = theta,
                        epsilon = eps_vals, L = L,
                        logPOSTERIOR = log_q_beta,
                        glogPOSTERIOR = g_log_q_beta,
                        param = list(Y = Y, tree = tree, beta = beta_temp), chains = 1)$theta
        
      } 
      hmc_samples <- MCMC_samples(as.numeric(samples[dim(samples)[1],]))
      #stopCluster(cl1)
      samples <- rbind.data.frame(do.call(rbind, hmc_samples[[1]])[(floor(N/3) + 1): N,], # use second two thirds of chain
                                  do.call(rbind, hmc_samples[[2]])[(floor(N/3) + 1): N,], # use second two thirds of chain
                                  do.call(rbind, hmc_samples[[3]])[(floor(N/3) + 1): N,]) # use second two thirds of chain
      
    }
  }
  stopCluster(cl1)
  print(paste0("eta_vec : ", eta_vec))
  print(paste0("r_vec : ", r_vec))
  log_rhat <- sum(eta_vec) + sum(log(r_vec)) 
  return(list(estimate = log_rhat, r = r_vec, eta = eta_vec) )
}

params <- c(rnorm(1,5,0.1), rnorm(1,5, 0.1), 0.2, 0.2,
            0.2, 0.2, -0.2, 0.2 )

set.seed(143)
L <- 15
eps_vals <- c(0.02, 0.02, rep(0.002, 4), 0.02, 0.02)
#eps_vals <- c(0.02, 0.02, rep(0.0001, 4), 0.02, 0.02) # with weakly informative priors
length(eps_vals)
N <- 500
K = 10
beta <- qbeta(seq(0, 1, length.out = K), 0.3, 1)
# mu1 mu2, gamma11, gamma21, gamma22, R11, R21, R22
initial_params <- params

unregister_dopar()
system.time(log_r <- GSS_HMC(Y = X, tree = tree, N = N, K, beta = beta, 
        model_log_posterior = LogLik_simple_drift_posterior,
        model_gradient = g_LogLik_simple_drift_posterior_simple,
        initial_params = params, L = L, eps_vals = eps_vals ))

log_r

save.image(file='GSS_model2_new.RData')





#### METROPOLIS sampling

# Main function
GSS_MH <- function(Y, tree, N, K, beta, model_log_posterior, initial_params){

  #making design matricees
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
  C_mat_inv<- solve(C_mat)
  
  # Make vectors to save elements of estimator
  r_vec <- c()
  eta_vec <- c()
  
  
  cl <- parallel::makeCluster(detectCores()-1, outfile = "model2.txt")
  registerDoParallel(cl)
  # parallel computation of MCMC samples
    MCMC_init <- function(theta){
    foreach(times=1:2,
            .combine = rbind,
            .export = c('model_log_posterior', 'LogLik_simple_drift_posterior',
                        'N', 'Y', 'tree', 'Design_matrix', 'C_mat', 'C_mat_inv')) %dopar%
      MCMCpack::MCMCmetrop1R(fun = model_log_posterior,
                   Y = Y, tree = tree, D = Design_matrix, C = C_mat, C_inv = C_mat_inv,
                   theta.init = theta,
                   seed = sample(1000,1),
                   burnin = 5000,
                   mcmc = N*10,
                   thin = 10,
                   #optim.lower = c(NA,NA, 0.01, 0.01, 0.01, 0.01, NA,NA),
                   verbose = 0, logfun = TRUE)}

  samples <- as.data.frame(MCMC_init(initial_params))
  stopCluster(cl)

  # find the sample mean and sds to parametrize reference distribution
  sample_means <- sapply(samples, mean)
  sample_sd <- sapply(samples, sd)
  sample_sd <- ifelse(sample_sd <= 0.001, 0.001, sample_sd)
  
  # define new reference distribution function
  ref_dist <- function(theta, sm = sample_means, sd = sample_sd){
    # both entries of mu has normal ref dist and prior
    mu_log_lik <- (dnorm(theta[1], mean = sm[1], sd = (sd[1]), log = TRUE) 
                   + dnorm(theta[2], mean = sm[2], sd = (sd[2]), log = TRUE))  
    # variances gamma(2,2) prior, covariance beta(2,5) prior
    gamma_log_lik <- (dgamma(theta[3], shape = (sm[3]^2)/sd[3]^2, scale = sd[3]^2/sm[3], log = TRUE)
                      + dgamma(theta[4], shape = (sm[4]^2)/(sd[4]^2), scale = sd[4]^2/sm[4], log = TRUE))
   # variances gamma(2,2) prior, covariance beta(2,5) prior
    R_log_lik <- (dgamma(theta[5], shape = (sm[5]^2)/sd[5]^2, scale = sd[5]^2/sm[5], log = TRUE)
                  + dgamma(theta[6], shape = (sm[6]^2)/sd[6]^2, scale = sd[6]^2/sm[6], log = TRUE))
    
    v_log_lik <- (dnorm(theta[7], mean = sm[7], sd = (sd[7]), log = TRUE) 
                  + dnorm(theta[8], mean = sm[8], sd = (sd[8]), log = TRUE)) 
    
    return(exp(mu_log_lik + gamma_log_lik + R_log_lik))            
  }
  
  # Define log_q_beta which we will sample from in GSS
  log_q_beta <- function(Y, tree, theta, beta, sm = sample_means, sd = sample_sd,
                         log_posterior = model_log_posterior, C=C_mat, D = Design_matrix, C_inv = C_mat_inv){
    return(beta * log_posterior(Y, tree, theta, D, C, C_inv) + (1-beta)*log(ref_dist(theta)))
  }
  # define function that calculates posterior divided by reference
  calc_frac <- function(Y, tree, theta, log_posterior = model_log_posterior, ref = ref_dist,  C = C_mat, D = Design_matrix, C_inv = C_mat_inv){
    fractions <- exp(log_posterior(Y, tree, theta, D,C, C_inv) - log(ref(theta)))
    return(fractions)
  }
  cl1 <- parallel::makeCluster(detectCores()-1, outfile = "model2.txt")
  registerDoParallel(cl1)
  
  for (k in (2:K)){ # We sum for beta_1 in the estimator which is the same as beta[2] 
    # calculate term for estimator
    fracs <-  apply(samples, 1, function(x) calc_frac(Y, tree, as.numeric(x)) )
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

      MCMC_samples <- function(theta){
        foreach(times=1:2,
                .combine = rbind,
                .export = c('log_q_beta',  'LogLik_simple_drift_posterior',
                            'model_log_posterior', 'N', 'Y', 'tree', 'beta_temp',
                            'ref_dist', 'sample_means', 'sample_sd', 'Design_matrix', 'C_mat', 'C_mat_inv')) %dopar%
          MCMCpack::MCMCmetrop1R(fun = log_q_beta,
                       Y = Y, tree = tree, D = Design_matrix, C = C_mat, C_inv = C_mat_inv,
                       beta = beta_temp, sm = sample_means, sd = sample_sd,
                       log_posterior = model_log_posterior,
                       theta.init = theta,
                       seed = sample(1000,1),
                       burnin = 000,
                       mcmc = N*10,
                       optim.lower = c(NA,NA, 0.01, 0.01, 0.01, 0.01, NA,NA),
                       thin = 10,
                       verbose = 0, logfun = TRUE)
        }
      samples <- as.data.frame(MCMC_samples(initial_params))
    }
  }
  stopCluster(cl1)
  print(paste0("eta_vec : ", eta_vec))
  print(paste0("r_vec : ", r_vec))
  log_rhat <- sum(eta_vec) + sum(log(r_vec)) 
  return(log_rhat)
}


N <- 500
theta <- params
set.seed(143)
K = 10
beta <- qbeta(seq(0, 1, length.out = K), 0.3, 1)
initial_params <- params

system.time(log_rhat_mh2 <- GSS_MH(Y = Y, 
                            tree = tree,
                            N = N, K = K, beta = beta,
                            model_log_posterior = LogLik_simple_drift_posterior,
                            initial_params = initial_params))



log_rhat_mh2

accept <- c(0.203, 0.1978, 0.19840,0.19740,0.19560,0.19840,0.20820,
            0.19620, 0.20700, 0.20200, 0.19760, 0.20460, 0.20620, 0.20240,
            0.20740, 0.18300)
mean(accept)

## Trace plot and marginal posterior plots
# first make the samples
N <- 500


cl <- parallel::makeCluster(detectCores())
registerDoParallel(cl)
# parallel computation of MCMC samples
MCMC_init <- function(theta){
  foreach(times=1:2,
          .combine = rbind,
          .export = c('LogLik_simple_drift_posterior',
                      'N', 'Y', 'tree', 'Design_matrix', 'C_mat', 'C_mat_inv')) %dopar%
    MCMCpack::MCMCmetrop1R(fun = LogLik_simple_drift_posterior,
                           Y = Y, tree = tree, D = Design_matrix, C = C_mat, C_inv = C_mat_inv,
                           theta.init = theta,
                           seed = sample(1000,1),
                           burnin = 5000,
                           mcmc = N*10,
                           thin = 10,
                           optim.lower = c(NA,NA, 0.01, 0.01, 0.01, 0.01, NA, NA),
                           verbose = 0, logfun = TRUE)}

samples <- as.data.frame(MCMC_init(initial_params))
stopCluster(cl)
dim(samples)

mcmc.3chains <- mcmc.list(as.mcmc(samples[1:N, ]),as.mcmc(samples[(N+1):(2*N), ]))#, as.mcmc(samples[((2*N)+1):(3*N), ]))
varnames(mcmc.3chains) <- c( "mu1", "mu2", "gamma1", "gamma2", "R11", "R22", "v1", "v2")

# trace plot
#plot(mcmc.3chains, density = FALSE, smooth = FALSE, lty = 1, col = c("lightskyblue", "lightskyblue4", "lightskyblue2"))
bayesplot::color_scheme_set("blue")
bayesplot::mcmc_trace(mcmc.3chains)

# marginal histograms
#samples <- rbind.data.frame(samples1, samples2, samples3)
#samples <- as.data.frame(samples)
colnames(samples) <- c( "mu1", "mu2", "gamma11", "gamma22", "R11", "R22", "v1", "v2")

samples_long <- samples %>%                          # Apply pivot_longer function
  pivot_longer(colnames(samples)) %>% 
  as.data.frame()
dim((samples_long))

dummy <- data.frame(name = c( "mu1", "mu2", "gamma11", "gamma22", "R11", "R22", "v1", "v2"),
                     Z = c(4,6,NA,NA, NA,NA, -0.2, 0.2))
#dummy$names <- factor(dummy$names)
ggp1 <- ggplot(samples_long, aes(x = value, fill = name)) +    # Draw each column as histogram
  geom_histogram() + 
  facet_wrap(~ name, scales = "free") +
  theme_bw() +
  geom_vline(data = dummy, aes(xintercept = Z), color = "red")

ggp1



