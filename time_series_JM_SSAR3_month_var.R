#time series analysis


rm(list=ls())

tic <- Sys.time()
Sys <- Sys.info()
source('Dc_Indonesia_nesting_fcns.R')
library(rjags)

save.RData <- T
save.fig <- F

MCMC.params <- list(n.chains = 3,
                    n.iter = 50000,
                    n.adapt = 25000)

# get JM data first:
data.0.JM <- read.csv('data/NestCounts_JM_09Feb2018.csv')

#data.0 <- read.csv("data/NestCount_Warmon_27March2018.csv")


# create time-duration filed (in yrs)
# define dates with begin and end dates:
data.0.JM %>% reshape2::melt(id.vars = "YEAR",
                             variable.name = "month",
                             value.name = "count") -> data.1.JM
data.1.JM$MONTH <- unlist(lapply(data.1.JM$month, FUN = mmm2month))

data.1.JM %>% mutate(., f.month = as.factor(MONTH),
                     f.year = as.factor(YEAR))%>%
  filter(YEAR > 1999) %>%
  mutate(Frac.Year = YEAR + (MONTH-0.5)/12) %>%
  reshape::sort_df(.,vars = "Frac.Year") -> data.2

#data.1.JM.2005 <- filter(data.1.JM, YEAR > 2004)

bugs.data <- list(y = data.2$count,
                  m = data.2$MONTH,
                  T = nrow(data.2))

# inits.function <- function(){
#   mu <- rnorm(1, 0, 10)
#   theta <- rnorm(1, 0, 1)
#   phi <- rnorm(1, 0, 1)
#   #sigma.pro <- runif(1, 0, 50)
#   #sigma.obs <- runif(1, 0, 50)
#   A <- list(mu = mu, theta = theta)
#   #          sigma.pro = sigma.pro, sigma.obs = sigma.obs)
#   return(A)
# }

load.module('dic')
params <- c('theta1', 'theta2', 'theta3',
            'sigma.pro1', 'sigma.pro2', 'sigma.obs', 'mu')

jm <- jags.model(file = 'models/model_SSAR3_month_var.txt',
                 data = bugs.data,
                 #inits = inits.function,
                 n.chains = MCMC.params$n.chains,
                 n.adapt = MCMC.params$n.adapt)

# check for convergence first. This doesn't work with Xs and ys
# for Gelman diagnostics - because I think existing data points
# don't have any varibility in posterior samples so the matrix
# cannot be inverted to compute variances
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)
g.diag <- gelman.diag(zm)

# then sample y
params <- c(params, 'y', 'X', 'deviance')
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)

summary.zm <- summary(zm)

# extract ys
ys.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'y[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(ys.stats) <- c('low_y', 'mode_y', 'high_y')
ys.stats$time <- data.2$Frac.Year
ys.stats$obsY <- data.2$count
ys.stats$month <- data.2$MONTH
ys.stats$year <- data.2$YEAR

# extract Xs - the state model
Xs.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'X[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(Xs.stats) <- c('low_X', 'mode_X', 'high_X')
Xs.stats$time <- data.2$Frac.Year
Xs.stats$obsY <- data.2$count
Xs.stats$month <- data.2$MONTH
Xs.stats$year <- data.2$YEAR

Xs.year <- group_by(Xs.stats, year) %>% summarize(mode = sum(mode_X),
                                                  low = sum(low_X),
                                                  high = sum(high_X))

p.1 <- ggplot() +
  #geom_point(data = ys.stats,
  #           aes(x = time, y = mode_y), color = "blue") +
  #geom_line(data = Xs.stats,
  #          aes(x = time, y = mode_X), color = 'blue') +
  geom_line(data = Xs.stats,
            aes(x = time, y = high_X), color = "red",
            linetype = 2) +
  geom_point(data = Xs.stats,
             aes(x = time, y = mode_X), color = "red",
             alpha = 0.5, size = 2) +
  geom_line(data = Xs.stats,
            aes(x = time, y = mode_X), color = "red",
            alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = low_X), color = "red",
            linetype = 2) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5) +
  labs(title = '', x = '', y = "Nest counts")

toc <- Sys.time()
dif.time <- toc - tic

results.JM_SSAR1_month <- list(data.1 = data.2,
                               bugs.data = bugs.data,
                               summary.zm = summary.zm,
                               Xs.stats = Xs.stats,
                               Xs.year = Xs.year,
                               ys.stats = ys.stats,
                               zm = zm,
                               tic = tic,
                               toc = toc,
                               dif.time = dif.time,
                               Sys = Sys,
                               MCMC.params = MCMC.params,
                               g.diag = g.diag,
                               jm = jm)
if (save.fig)
  ggsave(plot = p.1,
         filename = 'figures/predicted_counts_SSAR3_JM_month_var_2000.png',
         dpi = 600)

if (save.RData)
  save(results.JM_SSAR1_month,
       file = paste0('RData/SSAR3_month_JM_var_2000_',
                     Sys.Date(), '.RData'))


# plot posterior densities using bayesplot functions:
# get the ggplot2 base theme:

base_theme <- ggplot2::theme_get()
library(bayesplot)

# set back to the base theme:
ggplot2::theme_set(base_theme)
#
mcmc_dens(zm, c('theta1', 'theta2'))
mcmc_dens(zm, c('sigma.pro1', 'sigma.pro2'))

#mcmc_trace(zm, 'theta')
#mcmc_dens(zm, 'phi1')
#mcmc_dens(zm, 'phi2')
#mcmc_trace(zm, 'phi2')
#mcmc_dens(zm, 'sigma.pro2')
#mcmc_trace(zm, 'sigma')
#mcmc_dens(zm, 'sigma.obs')
