library(rstan)
library(loo)
set.seed(1)
y = rnorm(20)
N = length(y)
m = sample(1:12,size=length(y), replace=T)
y[c(3,13,14)] =NA
notna = which(!is.na(y))
y = y[notna]# drop NAs in y
n_notna = length(notna)
z = ifelse(m < 10 & m > 4, 1, 2)


data_list = list(y = y, N = N, notna = notna, n_notna = n_notna,
z = z)
mod = stan("tomo.stan", data=data_list, chains=1, iter=1000)

log_lik1 <- extract_log_lik(mod, merge_chains = FALSE)
rel_n_eff <- relative_eff(exp(log_lik1))
loo(log_lik1, r_eff = rel_n_eff, cores = 2)