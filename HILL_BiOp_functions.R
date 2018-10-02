#Dc_Indonesia_nesting_fcns

# these are functions that are used in this project.
# TomosFunctions.R is needed also.

#sysInfo <- Sys.info()
ifelse(Sys.info()[1] == 'Linux',
       source('~/Documents/R/tools/TomosFunctions.R'),
       source('~/R/tools/TomosFunctions.R'))

# library(ggplot2)
# library(tidyverse)
# library(lubridate)
# library(loo)

# extract.samples in TomosFunctions.R
sum.posterior <- function(yr, months = c(1:12), Xs.stats, zm) {
  Xs.stats %>%
    mutate(var.name = rownames(Xs.stats)) %>%
    filter(year == yr) %>%
    filter(month %in% months) %>%
    select(var.name) -> Xs.name
  zm.yr <- apply(Xs.name,
                 MARGIN = 1,
                 FUN = extract.samples,
                 zm)

  return(list(samples = zm.yr, var.names = Xs.name))
}

pareto.k.diag <- function(jm, MCMC.params, jags.data){
  
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
<<<<<<< HEAD
  loglik.obs <- jm$sims.list$loglik[, !is.na(jags.data$y)]
  # get rid of NA columns - even if data existed (for example the first value) - no likelihood
  # for the first data point
  loglik.obs <- loglik.obs[, colSums(is.na(loglik.obs)) == 0]
=======
  loglik.obs <- jm$sims.list$loglik[, 2:jags.data$T]
  
>>>>>>> 6bd6ecf3e1a4447aa2dd88f258ea4d68be8d5fa8
  Reff <- relative_eff(exp(loglik.obs), 
                       chain_id = rep(1:MCMC.params$n.chains, each = n.per.chain))
                       chain_id = rep(1:MCMC.params$n.chains, 
                                      each = n.per.chain))
  
  loo.out <- loo(loglik.obs, r_eff = Reff)
  return(list(loglik.obs = loglik.obs,
              Reff = Reff,
              loo.out = loo.out))
}
