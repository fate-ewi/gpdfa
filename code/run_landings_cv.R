library(bayesdfa)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(rerddap)
library(loo)
library(rstan)

options(mc.cores = parallel::detectCores())

d <- read.csv("data/port_landings_table2.csv", stringsAsFactors = FALSE)
d <- dplyr::select(d, -Year)
for (i in 1:ncol(d)) {
  d[, i] <- log(as.numeric(d[, i]))
}

x <- data.frame(
  "obs" = c(as.matrix(d)), "time" = rep(1:nrow(d), ncol(d)),
  "ts" = sort(rep(1:ncol(d), nrow(d)))
)

# x$ts = as.numeric(as.factor(x$ts))
# x$time = x$time - min(x$time) + 1
mcmc_iter <- 4000
mcmc_chains <- 3
m <- list()
num_trends <- 2
fold_ids <- x$time - min(x$time) + 1

fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, sample = FALSE)
m[[1]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = 1, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
# fit models with varying numbers of knots in a b-spline
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "spline", n_knots = 6, sample = FALSE)
m[[2]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "spline", n_knots = 12, sample = FALSE)
m[[3]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "spline", n_knots = 18, sample = FALSE)
m[[4]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "spline", n_knots = 24, sample = FALSE)
m[[5]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "spline", n_knots = 30, sample = FALSE)
m[[6]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
# fit models with varying numbers of knots in a gp
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "gp", n_knots = 6, sample = FALSE)
m[[7]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "gp", n_knots = 12, sample = FALSE)
m[[8]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "gp", n_knots = 18, sample = FALSE)
m[[9]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "gp", n_knots = 24, sample = FALSE)
m[[10]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
fit <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = mcmc_chains, num_trends = num_trends, trend_model = "gp", n_knots = 30, sample = FALSE)
m[[11]] <- dfa_cv(fit, cv_method = "loocv", iter = mcmc_iter, chains = mcmc_chains, fold_ids = fold_ids)
saveRDS(m, "output/landings_models_cv.rds")
