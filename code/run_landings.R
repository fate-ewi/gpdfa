library(bayesdfa)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)

d = read.csv("data/port_landings_table2.csv", stringsAsFactors = FALSE)
d = dplyr::select(d, -Year)
for(i in 1:ncol(d)) {
  d[,i] = log(as.numeric(d[,i]))
}

mcmc_iter = 3000
mcmc_chains = 3
m = list()
m[[1]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 1)
m[[2]] = fit_dfa(y = t(d), iter = 1000, chains = 1, num_trends = 2)

# fit models with varying numbers of knots in a b-spline
m[[3]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=6)
m[[4]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=12)
m[[5]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=18)
m[[6]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=24)
m[[7]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=30)

# fit models with varying numbers of knots in a b-spline
m[[8]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=6)
m[[9]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=12)
m[[10]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=18)
m[[11]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=24)
m[[12]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=30)

fit = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=18)
saveRDS(fit,"landings_spline_18.rds")
