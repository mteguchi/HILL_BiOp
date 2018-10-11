#time series analysis


rm(list=ls())
tic <- Sys.time()
Sys <- Sys.info()
source('HILL_BiOp_functions.R')
library(jagsUI)
library(coda)
library(ggplot2)
library(loo)


save.RData <- T
save.fig <- T
plot.fig <- T

MCMC.n.chains <- 5
MCMC.n.samples <- 500000
MCMC.n.burnin <- 350000
MCMC.n.thin <- 50

MCMC.params <- list(n.chains = MCMC.n.chains,
                    n.samples = MCMC.n.samples,
                    n.burnin = MCMC.n.burnin,
                    n.thin = MCMC.n.thin)

# get JM data first:
data.0.JM <- read.csv('data/JM_nests.csv')
# create regularly spaced time series:
data.2.JM <- data.frame(Year = rep(min(data.0.JM$Year_begin,
                                       na.rm = T):max(data.0.JM$Year_begin,
                                                      na.rm = T),
                                   each = 12),
                        Month_begin = rep(1:12,
                                          max(data.0.JM$Year_begin,
                                              na.rm = T) -
                                            min(data.0.JM$Year_begin,
                                                na.rm = T) + 1)) %>%
  mutate(begin_date = as.Date(paste(Year,
                                    Month_begin,
                                    '01', sep = "-"),
                              format = "%Y-%m-%d"),
         Frac.Year = Year + (Month_begin-0.5)/12) %>%
  select(Year, Month_begin, begin_date, Frac.Year)

data.0.JM %>% mutate(begin_date = as.Date(paste(Year_begin,
                                                Month_begin,
                                                '01', sep = "-"),
                                          format = "%Y-%m-%d")) %>%
  mutate(Year = Year_begin,
         Month = Month_begin,
         f_month = as.factor(Month),
         f_year = as.factor(Year),
         Frac.Year = Year + (Month_begin-0.5)/12,
         Nests = JM.1) %>%
  select(Year, Month, Frac.Year, begin_date, Nests) %>%
  na.omit() %>%
  right_join(.,data.2.JM, by = "begin_date") %>%
  transmute(Year = Year.y,
            Month = Month_begin,
            Frac.Year = Frac.Year.y,
            Nests = Nests) %>%
  reshape::sort_df(.,vars = "Frac.Year") %>%
  filter(Year > 1998) -> data.1.JM
#data.1.JM.2005 <- filter(data.1.JM, YEAR > 2004)

nests <- reshape2::acast(data.1.JM, Year ~ Month, value.var = "Nests")
months <- matrix(seq(1, 12), nrow = nrow(nests), ncol = 12, byrow = T)

jags.data <- list(y = log(nests),
                  m = months,
                  T = 12,
                  T0 = nrow(nests))

#load.module('dic')
jags.params <- c("r", 'sigma.Z', "sigma.X", "sigma.obs",
                 "p.pro1", "p.pro2", "df", "y", "X", 
                 "deviance", "loglik")

jm <- jags(jags.data,
           inits = NULL,
           parameters.to.save= jags.params,
           model.file = 'models/model_SS_trend_logY_norm_norm_t.txt',
           n.chains = MCMC.n.chains,
           n.burnin = MCMC.n.burnin,
           n.thin = MCMC.n.thin,
           n.iter = MCMC.n.samples,
           DIC = T, parallel=T)

#g.diag1 <- gelman.diag(jm$samples)
Rhat <- jm$Rhat
n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
# the first year is not used in the loop; t starts from 2 to T, so we need to 
# remove the first year from the data to extract loglikelihoods

data.used <- jags.data$y[2:nrow(jags.data$y), ]
loglik.obs <- matrix(NA, nrow = n.per.chain * MCMC.params$n.chains, 
                     ncol = sum(!is.na(data.used)))

# there has to be a better way to do this but can't figure out now... 10/11/2018
for (k in 1:nrow(loglik.obs)){
  loglik.1 <- t(jm$sims.list$loglik[k, 2:nrow(jm$sims.list$loglik[k,,]),])
  loglik.2 <- loglik.1[t(!is.na(data.used))]
  loglik.obs[k,] <- loglik.2
}

Reff <- relative_eff(exp(loglik.obs), 
                     chain_id = rep(1:MCMC.params$n.chains, 
                                    each = n.per.chain),
                     cores = 1)

loo.out <- loo(loglik.obs, r_eff = Reff, cores = 1)

# extract ys
ys.stats <- data.frame(low_y = as.vector(t(jm$q2.5$y)),
                       median_y = as.vector(t(jm$q50$y)),
                       high_y = as.vector(t(jm$q97.5$y)),
                       time = data.1.JM$Frac.Year,
                       obsY = data.1.JM$Nests,
                       month = data.1.JM$Month,
                       year = data.1.JM$Year)


# extract Xs - the state model
Xs.stats <- data.frame(low_X = as.vector(t(jm$q2.5$X)),
                       median_X = as.vector(t(jm$q50$X)),
                       high_X = as.vector(t(jm$q97.5$X)),
                       time = data.1.JM$Frac.Year,
                       obsY = data.1.JM$Nests,
                       month = data.1.JM$Month,
                       year = data.1.JM$Year)

Xs.year <- group_by(Xs.stats, year) %>% summarize(median = sum(median_X),
                                                  low = sum(low_X),
                                                  high = sum(high_X))

toc <- Sys.time()
dif.time <- toc - tic

p.1 <- ggplot() +
  #geom_point(data = ys.stats,
  #           aes(x = time, y = mode_y), color = "blue") +
  #geom_line(data = Xs.stats,
  #          aes(x = time, y = mode_X), color = 'blue') +
  geom_line(data = Xs.stats,
            aes(x = time, y = exp(high_X)), color = "red",
            linetype = 2) +
  geom_point(data = Xs.stats,
             aes(x = time, y = exp(median_X)), color = "red",
             alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = exp(median_X)), color = "red",
            alpha = 0.5) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5)


results <- list(data.1 = data.1.JM,
                jags.data = jags.data,
                Xs.stats = Xs.stats,
                Xs.year = Xs.year,
                ys.stats = ys.stats,
                tic = tic,
                toc = toc,
                dif.time = dif.time,
                Sys = Sys,
                MCMC.params = MCMC.params,
                Rhat = Rhat,
                jm = jm,
                loo.out = loo.out)
if (save.fig)
  ggsave(plot = p.1,
         filename = 'figures/predicted_counts_JM_logY_norm_norm_t_1999.png',
         dpi = 600)

if (save.RData)
  saveRDS(results,
          file = paste0('RData/SSAR1_logY_norm_norm_t_JM_1999_', Sys.Date(), '.rds'))

if (plot.fig){
  base_theme <- ggplot2::theme_get()
  library(bayesplot)

  # set back to the base theme:
  ggplot2::theme_set(base_theme)
  mcmc_trace(jm$samples, c("r", 'sigma.Z', "sigma.X", "sigma.obs",
                           "p.pro1", "p.pro2", "df"))
  mcmc_dens(jm$samples, c("r", 'sigma.Z', "sigma.X", "sigma.obs",
                           "p.pro1", "p.pro2", "df"))

}
