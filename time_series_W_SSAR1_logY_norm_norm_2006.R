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
data.0 <- read.csv('data/W_nests.csv')
# create regularly spaced time series:
data.2 <- data.frame(Year = rep(min(data.0$Year_begin,
                                       na.rm = T):max(data.0$Year_begin,
                                                      na.rm = T),
                                   each = 12),
                        Month_begin = rep(1:12,
                                          max(data.0$Year_begin,
                                              na.rm = T) -
                                            min(data.0$Year_begin,
                                                na.rm = T) + 1)) %>%
  mutate(begin_date = as.Date(paste(Year,
                                    Month_begin,
                                    '01', sep = "-"),
                              format = "%Y-%m-%d"),
         Frac.Year = Year + (Month_begin-0.5)/12) %>%
  select(Year, Month_begin, begin_date, Frac.Year)

data.0 %>% mutate(begin_date = as.Date(paste(Year_begin,
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
  right_join(.,data.2, by = "begin_date") %>%
  transmute(Year = Year.y,
            Month = Month_begin,
            Frac.Year = Frac.Year.y,
            Nests = Nests) %>%
  reshape::sort_df(.,vars = "Frac.Year") %>%
  filter(Year > 2005) -> data.1
#data.1.2005 <- filter(data.1, YEAR > 2004)

data.1 <- data.1[5:nrow(data.1),]

jags.data <- list(y = log(data.1$Nests),
                  m = data.1$Month,
                  T = nrow(data.1))

#load.module('dic')
jags.params <- c('theta.1', 'sigma.pro1', "sigma.obs",
                 'mu', 'y', 'X', 'deviance', 'loglik')

jm <- jags(jags.data,
           inits = NULL,
           parameters.to.save= jags.params,
           model.file = 'models/model_SSAR1_logY_norm_norm.txt',
           n.chains = MCMC.n.chains,
           n.burnin = MCMC.n.burnin,
           n.thin = MCMC.n.thin,
           n.iter = MCMC.n.samples,
           DIC = T, parallel=T)

#g.diag1 <- gelman.diag(jm$samples)
Rhat <- jm$Rhat

# extract ys
ys.stats <- data.frame(low_y = jm$q2.5$y,
                       median_y = jm$q50$y,
                       high_y = jm$q97.5$y,
                       time = data.1$Frac.Year,
                       obsY = data.1$Nests,
                       month = data.1$Month,
                       year = data.1$Year)


# extract Xs - the state model
Xs.stats <- data.frame(low_X = jm$q2.5$X,
                       median_X = jm$q50$X,
                       high_X = jm$q97.5$X,
                       time = data.1$Frac.Year,
                       obsY = data.1$Nests,
                       month = data.1$Month,
                       year = data.1$Year)

Xs.year <- group_by(Xs.stats, year) %>% summarize(median = sum(median_X),
                                                  low = sum(low_X),
                                                  high = sum(high_X))

loo.out <- pareto.k.diag(jm, MCMC.params, jags.data)
  
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


results <- list(data.1 = data.1,
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
         filename = 'figures/predicted_counts_W_logY_norm_norm_2006.png',
         dpi = 600)

if (save.RData)
  saveRDS(results,
          file = paste0('RData/SSAR1_logY_norm_norm_W_2006_', Sys.Date(), '.rds'))

if (plot.fig){
  base_theme <- ggplot2::theme_get()
  library(bayesplot)

  # set back to the base theme:
  ggplot2::theme_set(base_theme)
  mcmc_trace(jm$samples, c('theta.1', "mu", "sigma.obs",
                           "sigma.pro1"))
  mcmc_dens(jm$samples, c('theta.1', "mu", "sigma.obs",
                          "sigma.pro1"))

}
