# Use loo package on simulated AR(1) data 

# the underlying model is essentially an AR(1). However, the slope and variance 
# change according to months. See below for details.

library(jagsUI)
library(bayesplot)
library(loo)

theta.1 <- 1.0005  # increasing months
sigma.1 <- 10  # SD 1
mu <- 400       # X[1]
r <- 30         # Negative Binomial "size" parameter 

m <- rep(c(1:12), 21)   # months of each year

# simulate data:
set.seed(123)
X <- p <- y <- vector(mode = "numeric", length = length(m))
X[1] <- mu
p[1] <- r/(r + X[1])
y[1] <- rnbinom(1, size = r, prob = p[1])

for (k in 2:length(m)){

  X.mean <- theta.1 * X[k-1]
  
  # Could have used log-normal to avoid negative values ... 
  X[k] <- round(abs(rnorm(1, mean = X.mean, sd = sigma.1)))

  # compute the "probability" for the negative binomial
  #p[k] <- r/(r + X[k])  # this is the inverse... no?
  p[k] <- X[k]/(r + X[k])
  
  # create "Observed" value.
  y[k] <- rnbinom(1, size = r, prob = p[k])
}

# Run jags on the recorded data (y)

MCMC.n.chains <- 5
MCMC.n.samples <- 50000
MCMC.n.burnin <- 35000
MCMC.n.thin <- 1

jags.data <- list(y = y,
                  m = m,
                  T = length(m))

jags.params <- c('theta.1', 
                 'sigma.pro', "r",
                 'mu', 'X', 'deviance', 'loglik')

jags.model <- cat("model{
	for (t in 2:T){
    # process: 
    predX[t] <- theta.1 * X[t-1]
                  
    # one variance 
    X[t] ~ dnorm(predX[t], tau.pro)T(0,)
                  
    # observation
    y[t] ~ dnegbin(p[t], r)
    p[t] <- X[t]/(r + X[t])
                  
    loglik[t] <- logdensity.negbin(y[t], p[t], r)
                  
 }
                  
 X[1] <- mu
 y[1] ~ dnegbin(p[1], r)
 p[1] <- X[1]/(r + X[1])
                  
 mu ~ dunif(350, 450)
                  
 sigma.pro ~ dgamma(0.1, 0.01)
 tau.pro <- 1/(sigma.pro * sigma.pro)
                  
 r ~ dgamma(0.1, 0.01)
                  
 theta.1 ~ dnorm(0, 0.01)
}", file = "jags_model.txt")

jm <- jags(jags.data,
           inits = NULL,
           parameters.to.save= jags.params,
           model.file = 'jags_model.txt',
           n.chains = MCMC.n.chains,
           n.burnin = MCMC.n.burnin,
           n.thin = MCMC.n.thin,
           n.iter = MCMC.n.samples,
           DIC = T, parallel=T)

n.per.chain <- (MCMC.n.samples - MCMC.n.burnin)/MCMC.n.thin

loglik.obs <- jm$sims.list$loglik[, !is.na(jags.data$y)]
loglik.obs <- loglik.obs[, colSums(is.na(loglik.obs)) == 0]

Reff <- relative_eff(exp(loglik.obs), 
                     chain_id = rep(1:MCMC.n.chains, 
                                    each = n.per.chain))

loo.out <- loo(loglik.obs, r_eff = Reff)

save(list = ls(), file = "RData/loo_test_norm_negbin.RData")

