#time series analysis


rm(list=ls())

tic <- Sys.time()
Sys <- Sys.info()
source('Dc_Indonesia_nesting_fcns.R')
library(rjags)
library(jagsUI)

save.RData <- T
save.fig <- F

MCMC.params <- list(n.chains = 3,
                    n.iter = 50000,
                    n.adapt = 100000)

data.0 <- read.csv("data/NestCounts_Warmon_27March2018.csv")

# create time-duration filed (in yrs)
# define dates with begin and end dates:
data.0 %>% reshape2::melt(id.vars = "YEAR",
                          variable.name = "month",
                          value.name = "count") -> data.1
data.1$MONTH <- unlist(lapply(data.1$month, FUN = mmm2month))

data.2 <- mutate(data.1, f.month = as.factor(MONTH),
                 f.year = as.factor(YEAR))%>%
  filter(YEAR < 2014 & YEAR > 2005) %>%
  #filter(YEAR > 2005) %>%
  #filter(YEAR < 2014) %>%
  mutate(Frac.Year = YEAR + (MONTH-0.5)/12) %>%
  reshape::sort_df(.,vars = "Frac.Year")

#data.1.JM.2005 <- filter(data.1.JM, YEAR > 2004)

bugs.data <- list(y = data.2$count,
                  m = data.2$MONTH,
                  T = nrow(data.2))

params <- c('theta.1', 'theta.2',
            'sigma.pro1', 'sigma.pro2',
            'sigma.obs', 'mu')
load.module("dic")
load.module("glm")
# autojags from jagsUI seems to work pretty good - can't sample missing data
# points because of the same reason as before - convergence cannot be determined.
# parallel runs don't seem to work... 
# autojags.out <- autojags(data = bugs.data,
#                          inits = NULL,
#                          parameters.to.save = params,
#                          model.file = 'models/model_SSAR1_month_var_theta_Wermon.txt',
#                          n.chains = MCMC.params$n.chains,
#                          modules = c("glm", "dic"),
#                          #parallel = TRUE,
#                          #n.cores = 2,
#                          Rhat.limit = 1.1)

jm <- jags.model(file = 'models/model_SSAR1_month_var_theta_Wermon.txt',
                 data = bugs.data,
                 #inits = inits.function,
                 n.chains = MCMC.params$n.chains,
                 n.adapt = MCMC.params$n.adapt)

# check for convergence first.
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)
g.diag <- gelman.diag(zm)

# then sample y and X
params <- c(params, 'y', 'X', 'deviance')
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)

# params <- c(params, 'y', 'X')
# jags.out <- jags(data = bugs.data,
#                  inits = NULL,
#                  parameters.to.save = params,
#                  model.file = 'models/model_SSAR1_month_var_theta_Wermon.txt',
#                  n.chains = MCMC.params$n.chains,
#                  n.adapt = MCMC.params$n.adapt,
#                  n.iter = autojags.out$mcmc.info$n.chains,
#                  modules = c("glm", "dic"),
#                  parallel = TRUE)

summary.zm <- summary(zm)

# Computing WAIC. Code is from
# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/8211df61/
WAIC.s <- jags.samples(jm, c("deviance", "WAIC"), 
                       type="mean", 
                       n.iter=10000, thin=10)
WAIC.s <- lapply(WAIC.s, unclass)
sapply(WAIC.s, sum)

# According to the website, "the monitor called WAIC does 
# indeed return the WAIC penalty, so it is currently misnamed."  So that means
# we need to convert what it comes back to the actual WAIC.  Alternatively, 
# according to Matt Denwood: "In JAGS - I have been playing with the DIC module 
# (locally but I will push to the repo). I think the easiest approach is to set 
# monitors for elpd_waic and p_waic (named for consistency with the loo package) 
# and calculate waic from this in R. To be honest WAIC was the main driver behind 
# me adding the running variance monitor to JAGS last year but I never quite got 
# around to the next step before JAGS 4.3.0..."  But the folllowing code does
# not workj - elpd_waic and p_waic cannot be found...  
# waic <- jags.samples(jm, c("elpd_waic", "p_waic"), 
#                      type = "mean",
#                      n.iter=10000, thin=10)


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

results.W_SSAR1_month_var_theta <- list(data.1 = data.2,
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
                                        jm = jm,
                                        WAIC.penalty = WAIC.s)
if (save.fig)
  ggsave(plot = p.1,
         filename = 'figures/predicted_counts_W_month_var_theta_2006To2013.png',
         height = 6, width = 8, units = "in", dpi = 600)

if (save.RData)
  save(results.W_SSAR1_month_var_theta,
       file = paste0('RData/SSAR1_month_W_var_theta_2006To2013_',
                     Sys.Date(), '.RData'))


# plot posterior densities using bayesplot functions:
# get the ggplot2 base theme:

# base_theme <- ggplot2::theme_get()
# library(bayesplot)
#
# # set back to the base theme:
# ggplot2::theme_set(base_theme)
#
# mcmc_dens(zm, c('theta.1', 'theta.2'))
# #mcmc_trace(zm, 'theta')
# #mcmc_dens(zm, 'phi1')
# #mcmc_dens(zm, 'phi2')
# #mcmc_trace(zm, 'phi2')
# mcmc_dens(zm, c('sigma.pro1', 'sigma.pro2'))
# #mcmc_dens(zm, 'sigma.pro2')
# #mcmc_trace(zm, 'sigma')
# #mcmc_dens(zm, 'sigma.obs')
