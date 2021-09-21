# this one leaves out data points and tests on those

library(bayesdfa)
library(dplyr)
library(ggplot2)
library(future)
plan(multisession, workers = 8L)
source("code/sim_dfa2.R")

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

get_log_dens <- function(sim_dat, model, out_of_sample_i) {
  e <- rstan::extract(model$model)
  pred <- predicted(model)
  log_density <- 0
  for (j in 1:nrow(out_of_sample_i)) {
    for (i in 1:col(out_of_sample_i)) {
        log_density <- log_density +
          log_sum_exp(dnorm(
            sim_dat[j, out_of_sample_i[j, i]],
            mean = pred[, 1L, out_of_sample_i[j, i], j],
            sd = e$sigma,
            log = TRUE
          ))
    }
  }
  log_density
}

sim_and_fit <- function(sigma_obs = 0.3, sigma_loadings = 0.1,
  type_sim = c("RW", "ARMA")) {

  type_sim <- match.arg(type_sim)
  n_t <- 30
  sim <- sim_dfa2(
    num_trends = 1,
    num_years = n_t,
    num_ts = 3,
    type = type_sim,
    sigma = sigma_obs,
    loadings_matrix = matrix(nrow = 3, ncol = 1,
      rnorm(3 * 1, 1, sigma_loadings)))

  left_out_i <- matrix(sample(seq_len(n_t), round(n_t * 0.1)), nrow = 1)
  for (i in seq(2, 3)) {
    left_out_i <- rbind(left_out_i, sample(seq_len(n_t), round(n_t * 0.1)))
  }
  sim_orig <- sim
  for (j in 1:nrow(left_out_i)) {
    sim$y_sim[j, left_out_i] <- NA
  }

  m_rw <- fit_dfa(y = sim$y_sim, iter = 400, chains = 1, scale = "none")
  m_bs <- fit_dfa(
    y = sim$y_sim, iter = 400, chains = 1,
    trend_model = "bs", n_knots = 13, scale = "none"
  )
  elpd_rw <- get_log_dens(sim_orig$y_sim, m_rw, left_out_i)
  elpd_bs <- get_log_dens(sim_orig$y_sim, m_bs, left_out_i)
  data.frame(epld_rw = elpd_rw, elpd_bs = elpd_bs)
}

# sim_and_fit(type_sim = "RW")

x_rw <- furrr::future_map_dfr(1:8, ~ sim_and_fit(type_sim = "RW"), .id = "iter")
plot(x_rw$epld_rw - x_rw$elpd_bs);abline(h = 0, lty = 2)

x_arma <- furrr::future_map_dfr(1:8, ~ sim_and_fit(type_sim = "ARMA"), .id = "iter")
plot(x_arma$epld_rw - x_arma$elpd_bs);abline(h = 0, lty = 2)
