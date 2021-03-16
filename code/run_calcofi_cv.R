library(bayesdfa)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(rerddap)
library(loo)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

x <- readRDS("output/calcofi_data.rds")
# x$ts = as.numeric(as.factor(x$ts))
# x$time = x$time - min(x$time) + 1
mcmc_iter <- 4000
mcmc_chains <- 3
m <- list()
fold_ids <- x$time - min(x$time) + 1

# fit basic random walk model
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, sample = FALSE)
m[[1]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = 1, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")

# fit models with varying numbers of knots in a b-spline
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline", n_knots = 6, sample = FALSE)
m[[2]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline", n_knots = 12, sample = FALSE)
m[[3]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline", n_knots = 18, sample = FALSE)
m[[4]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline", n_knots = 24, sample = FALSE)
m[[5]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline", n_knots = 30, sample = FALSE)
m[[6]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "calcofi_models_cv.rds")

# fit models with varying numbers of knots in a gp
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 6, sample = FALSE)
m[[7]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "calcofi_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 12, sample = FALSE)
m[[8]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 18, sample = FALSE)
m[[9]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 24, sample = FALSE)
m[[10]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 30, sample = FALSE)
m[[11]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/calcofi_models_cv.rds")
