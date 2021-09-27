library(bayesdfa)
library(dplyr)
library(ggplot2)
library(future)
plan(multisession)
source("code/sim_dfa2.R")

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

get_log_dens <- function(sim_dat, model) {
  e_rw <- rstan::extract(model$model)
  pred_rw <- predicted(model)
  log_density <- 0
  for (i in 1:ncol(sim_dat)) {
    for (k in 1:nrow(sim_dat)) {
      log_density <- log_density +
        log_sum_exp(dnorm(
          sim_dat[k, i],
          mean = pred_rw[, 1L, i, k],
          sd = e_rw$sigma,
          log = TRUE
        ))
    }
  }
  log_density
}

sim_and_fit <- function(sigma_obs = 0.3, sigma_loadings = 0.2,
  type_sim = c("RW", "ARMA")) {

  type_sim <- match.arg(type_sim)
  sim <- sim_dfa2(
    num_trends = 1,
    num_years = 30,
    num_ts = 6,
    type = type_sim,
    sigma = sigma_obs,
    loadings_matrix = matrix(nrow = 6, ncol = 1, rnorm(6 * 1, 1, sigma_loadings)))

  m_rw <- fit_dfa(y = sim$y_sim, iter = 300, chains = 1, scale = "none")
  m_bs <- fit_dfa(
    y = sim$y_sim, iter = 300, chains = 1,
    trend_model = "bs", n_knots = 12, scale = "none"
  )
  # actual ELPD:
  elpd_rw <- get_log_dens(sim$y_sim, m_rw)
  elpd_bs <- get_log_dens(sim$y_sim, m_bs)
  data.frame(epld_rw = elpd_rw, elpd_bs = elpd_bs)
}

# sim_and_fit(type_sim = "ARMA")

x <- furrr::future_map_dfr(1:8, ~ sim_and_fit(), .id = "iter", type = "ARMA")
