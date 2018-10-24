#time series analysis

rm(list=ls())
tic <- Sys.time()
Sys <- Sys.info()
source('HILL_BiOp_functions.R')

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

year.begin <- 1999
year.end <- 2017
loc <- "JM"
data.jags <- data.extract(location = loc, 
                          year.begin = year.begin, 
                          year.end = year.end)

jags.params <- c('theta.1', 'sigma.pro1', 'sigma.pro2', "sigma.obs",
                 'mu', 'y', 'X', 'deviance', 'loglik')

model.file = 'models/model_SSAR1_logY_norm_norm_var.txt'

jags.out <- run.jagsUI(data.jags$jags.data, 
                       jags.params, 
                       model.file = model.file, 
                       MCMC.params)

Xs.stats <- jags.out$Xs.stats

Xs.stats$time <- data.jags$data.1$Frac.Year
Xs.stats$obsY <- data.jags$data.1$Nests
Xs.stats$month <- data.jags$data.1$Month
Xs.stats$year <- data.jags$data.1$Year

ys.stats <- jags.out$ys.stats
ys.stats$time <- data.jags$data.1$Frac.Year
ys.stats$obsY <- data.jags$data.1$Nests
ys.stats$month <- data.jags$data.1$Month
ys.stats$year <- data.jags$data.1$Year

results <- list(data.1 = data.jags$data.1,
                jags.out = jags.out,
                Xs.stats = Xs.stats,
                ys.stats = ys.stats,
                MCMC.params = MCMC.params)

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
  geom_line(data = Xs.stats,
            aes(x = time, y = exp(low_X)), color = "red",
            linetype = 2) +
  
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5)

filename.root <- strsplit(strsplit(jags.out$jm$modfile, 
                                   'models/model_')[[1]][2], '.txt')[[1]][1]

if (length(strsplit(filename.root, paste0("SSAR1_", loc, "_"))[[1]]) == 1){
  
  filename.out <- paste0(filename.root, "_", loc, "_", year.begin, "_", year.end)
  
} else {
  filename.root <- strsplit(filename.root, paste0("SSAR1_", loc, "_"))[[1]][2]
  filename.out <- paste0("SSAR1_", filename.root, "_", loc, "_", year.begin, "_", year.end)
}

if (save.fig)
  ggsave(plot = p.1,
         filename = paste0("figures/", "predicted_counts_", filename.out, ".png"),
         dpi = 600)

if (save.RData)
  saveRDS(results,
          file = paste0('RData/', "jagsout_", 
                        filename.out, "_", Sys.Date(), '.rds'))


if (plot.fig){
  base_theme <- ggplot2::theme_get()
  library(bayesplot)

  # set back to the base theme:
  ggplot2::theme_set(base_theme)
  mcmc_trace(jm$samples, c('theta.1', "mu", "sigma.obs",
                           "sigma.pro1", "sigma.pro2"))
  mcmc_dens(jm$samples, c('theta.1', "mu", "sigma.obs",
                          "sigma.pro1", "sigma.pro2"))

}
