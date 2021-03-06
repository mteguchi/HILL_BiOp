#time series analysis


rm(list=ls())

tic <- Sys.time()
Sys <- Sys.info()
source('HILL_BiOp_functions.R')
library(rjags)
library(bayesplot)

save.RData <- T
save.fig <- F

MCMC.params <- list(n.chains = 3,
                    n.iter = 50000,
                    n.adapt = 100000)

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

bugs.data <- list(y = data.1.JM$Nests,
                  m = data.1.JM$Month,
                  T = nrow(data.1.JM))

load.module('dic')
params <- c('theta', 'sigma.pro1', 'sigma.pro2', 'sigma.obs', 'mu')

jm <- jags.model(file = 'models/model_SSAR1_month.txt',
                 data = bugs.data,
                 #inits = inits.function,
                 n.chains = MCMC.params$n.chains,
                 n.adapt = MCMC.params$n.adapt)

# check for convergence first.
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
colnames(ys.stats) <- c('low_y', 'median_y', 'high_y')
ys.stats$time <- data.1.JM$Frac.Year
ys.stats$obsY <- data.1.JM$Nests
ys.stats$month <- data.1.JM$Month
ys.stats$year <- data.1.JM$Year

# extract Xs - the state model
Xs.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'X[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(Xs.stats) <- c('low_X', 'median_X', 'high_X')
Xs.stats$time <- data.1.JM$Frac.Year
Xs.stats$obsY <- data.1.JM$Nests
Xs.stats$month <- data.1.JM$Month
Xs.stats$year <- data.1.JM$Year

Xs.year <- group_by(Xs.stats, year) %>% summarize(mode = sum(median_X),
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
             aes(x = time, y = median_X), color = "red",
             alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = median_X), color = "red",
            alpha = 0.5) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5)

toc <- Sys.time()
dif.time <- toc - tic

results.JM_SSAR1_month <- list(data.1 = data.1.JM,
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
         filename = 'figures/predicted_counts_JM_month_1999.png',
         dpi = 600)

if (save.RData)
  save(results.JM_SSAR1_month,
       file = paste0('RData/SSAR1_month_JM_1999_', Sys.Date(), '.RData'))
