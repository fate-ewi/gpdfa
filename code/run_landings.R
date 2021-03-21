library(loo)
library(bayesdfa)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
options(mc.cores = parallel::detectCores())

d <- read.csv("data/port_landings_table2.csv", stringsAsFactors = FALSE)
d <- select(d, -Year)
for (i in 1:ncol(d)) {
  d[, i] <- log(as.numeric(d[, i]))
}

mcmc_iter <- 4000
mcmc_chains <- 3
thin <- 10
m <- list()
m[[1]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2,par_list="all")
loos[[1]] = loo::loo(m[[1]]$model,moment_match=TRUE)
saveRDS(loos, "output/landings_loos.rds")

# fit models with varying numbers of knots in a b-spline
m[[2]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "spline", n_knots = 6,par_list="all")
loos[[2]] = loo::loo(m[[2]]$model,moment_match=TRUE)
saveRDS(loos, "output/landings_loos.rds")

m[[3]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "spline", n_knots = 12,par_list="all")
loos[[3]] = loo::loo(m[[3]]$model,moment_match=TRUE)
saveRDS(loos, "output/landings_loos.rds")

m[[4]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "spline", n_knots = 18,par_list="all")
loos[[4]] = loo::loo(m[[4]]$model,moment_match=TRUE)
saveRDS(loos, "output/landings_loos.rds")

m[[5]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "spline", n_knots = 24,par_list="all")
loos[[5]] = loo::loo(m[[5]]$model,moment_match=TRUE)
saveRDS(loos, "output/landings_loos.rds")

m[[6]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "spline", n_knots = 30,par_list="all")
loos[[6]] = loo::loo(m[[6]]$model,moment_match=TRUE)
saveRDS(loos, "output/landings_loos.rds")

# fit models with varying numbers of knots in a b-spline
m[[7]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "gp", n_knots = 6,par_list="all")
loos[[7]] = loo::loo(m[[7]]$model,moment_match=TRUE)
saveRDS(m,"output/landings_models.rds")
saveRDS(loos, "output/landings_loos.rds")

m[[8]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "gp", n_knots = 12,par_list="all")
loos[[8]] = loo::loo(m[[8]]$model,moment_match=TRUE)
saveRDS(m,"output/landings_models.rds")
saveRDS(loos, "output/landings_loos.rds")

m[[9]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "gp", n_knots = 18,par_list="all")
loos[[9]] = loo::loo(m[[9]]$model,moment_match=TRUE)
saveRDS(m,"output/landings_models.rds")
saveRDS(loos, "output/landings_loos.rds")

m[[10]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "gp", n_knots = 24,par_list="all")
loos[[10]] = loo::loo(m[[10]]$model,moment_match=TRUE)
saveRDS(loos, "output/landings_loos.rds")
saveRDS(m,"output/landings_models.rds")

m[[11]] <- fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, thin = thin, num_trends = 2, trend_model = "gp", n_knots = 30,par_list="all")
loos[[11]] = loo::loo(m[[11]]$model,moment_match=TRUE)
saveRDS(loos, "output/landings_loos.rds")
saveRDS(m,"output/landings_models.rds")

