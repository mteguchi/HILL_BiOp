#time series analysis


rm(list=ls())

tic <- Sys.time()
Sys <- Sys.info()
source('Dc_Indonesia_nesting_fcns.R')
library(rjags)
library(bayesplot)

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

data.1 <- mutate(data.1, f.month = as.factor(MONTH),
                    f.year = as.factor(YEAR))%>%
  #filter(YEAR < 2014 & YEAR > 2005) %>%
  #filter(YEAR > 2005) %>%
  filter(YEAR < 2014) %>%
  mutate(Frac.Year = YEAR + (MONTH-0.5)/12) %>%
  reshape::sort_df(.,vars = "Frac.Year")

bugs.data <- list(y = data.1$count,
                  m = data.1$MONTH,
                  T = nrow(data.1))

# bugs.data <- list(y = data.1$count,
#                   T = 168)

inits.function <- function(){
  mu <- rnorm(1, 0, 10)
  theta <- rnorm(1, 0, 1)
  #phi <- rnorm(1, 0, 1)
  #sigma.pro <- runif(1, 0, 50)
  #sigma.obs <- runif(1, 0, 50)
  A <- list(mu = mu, theta = theta)
  #          sigma.pro = sigma.pro, sigma.obs = sigma.obs)
  return(A)
}

load.module('dic')
params <- c('theta', 'sigma.pro1', 'sigma.pro2',
            'sigma.obs')

jm <- jags.model(file = 'models/model_SSAR1_month_Warmon.txt',
                 data = bugs.data,
                 #inits = inits.function,
                 n.chains = MCMC.params$n.chains,
                 n.adapt = MCMC.params$n.adapt)

# check for convergence first.
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)
g.diag <- gelman.diag(zm)

# plot posterior densities using bayesplot functions:
# mcmc_dens(zm, 'theta')
# mcmc_dens(zm, 'sigma.pro')
# mcmc_dens(zm, 'sigma.obs')

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
ys.stats$time <- data.1$Frac.Year
ys.stats$obsY <- data.1$count
ys.stats$month <- data.1$MONTH
ys.stats$year <- data.1$YEAR

# extract Xs - the state model
Xs.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'X[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(Xs.stats) <- c('low_X', 'mode_X', 'high_X')
Xs.stats$time <- data.1$Frac.Year
Xs.stats$obsY <- data.1$count
Xs.stats$month <- data.1$MONTH
Xs.stats$year <- data.1$YEAR

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
  geom_line(data = Xs.stats,
            aes(x = time, y = low_X), color = "red",
            linetype = 2) +
  geom_point(data = Xs.stats,
             aes(x = time, y = mode_X), color = "red",
             alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = mode_X), color = "red",
            alpha = 0.5) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5)+
  geom_line(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5) +
  labs(x = '', y = '# nests')  +
  theme(axis.text = element_text(size = 12),
        text = element_text(size = 12))

toc <- Sys.time()
dif.time <- toc - tic

results.Warmon_SSAR1_month_To2013 <- list(data.1 = data.1,
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
         filename = 'figures/predicted_counts_SSAR1_month_Warmon_To2013.png',
         dpi = 600)

if (save.RData)
  save(results.Warmon_SSAR1_month_To2013,
       file = paste0('RData/SSAR1_month_Warmon_', Sys.Date(), '_To2013.RData'))
