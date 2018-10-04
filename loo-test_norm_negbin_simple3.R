# Use loo package on simulated AR(1) data 

# the underlying model is essentially an AR(1). However, the slope and variance 
# change according to months. See below for details.

library(jagsUI)
library(bayesplot)
library(loo)

theta.1 <- 1.004  # increasing months
theta.2 <- 0.84  # decreasing months
sigma.1 <- 123  # SD 1
sigma.2 <- 30   # SD 2
sigma.obs <- 20
mu <- 400       # X[1]
r <- 30         # Negative Binomial "size" parameter 

m <- rep(c(1:12), 21)   # months of each year

# simulate data:
set.seed(123)
X <- p <- y <- theta <- vector(mode = "numeric", length = length(m))
X[1] <- mu
p[1] <- r/(r + X[1])
y[1] <- rnbinom(1, size = r, prob = p[1])

for (k in 2:length(m)){
  # slope is theta.1 when month is Jan through July, after July, it is theta.2
  theta[k] <- ifelse(m[k] < 8, theta.1, theta.2)
  
  # SD is sigma.1 between May and September, it is sigma.2 otherwise
  #sigma <- ifelse(m[k] < 10 & m[k] > 4, sigma.1, sigma.2)
  
  X.mean <- theta[k] * X[k-1]
  
  # Could have used log-normal to avoid negative values ... 
  X[k] <- round(abs(rnorm(1, mean = X.mean, sd = sigma.2)))
  
  # compute the "probability" for the negative binomial
  #p[k] <- r/(r + X[k])
  
  # create "Observed" value.
  #y[k] <- rnbinom(1, size = r, prob = p[k])
  y[k] <- round(abs(rnorm(1, mean = X[k], sd = sigma.obs)))
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
    # process: increasing before August, decreasing after July
    #theta[t] <- ifelse(m[t] < 8, theta.1, theta.2)
    predX[t] <- theta.1 * X[t-1]
                  
    # one variance
    X[t] ~ dnorm(predX[t], tau.pro)T(0,)
                  
    # observation
    y[t] ~ dnegbin(p[t], r)
    p[t] <- r/(r + exp(X[t]))
                  
    loglik[t] <- logdensity.negbin(y[t], p[t], r)
                  
 }
                  
 X[1] <- mu
 y[1] ~ dnegbin(p[1], r)
 p[1] <- r/(r + X[1])
                  
 mu ~ dpois(400)
                  
 sigma.pro ~ dgamma(0.1, 0.01)
 tau.pro <- 1/(sigma.pro * sigma.pro)
                  
 r ~ dgamma(0.1, 0.01)
                  
 theta.1 ~ dnorm(0, 0.01)
 #theta.2 ~ dnorm(0, 0.01)
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
                     chain_id = rep(1:MCMC.params$n.chains, 
                                    each = n.per.chain))

loo.out <- loo(loglik.obs, r_eff = Reff)

