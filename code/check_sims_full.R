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

sim_and_fit <- function(sigma_obs = 0.5, sigma_loadings = 0.02,
  mean_loadings = 1, scale = "none",
  type_sim = c("RW", "ARMA"), num_ts = 3, n_t = 20, stan_iter = 300,
  n_knots = 10, seed = sample.int(.Machine$integer.max, 1L), ...) {
  type_sim <- match.arg(type_sim)
  set.seed(seed)
  sim <- sim_dfa2(
    num_trends = 1,
    num_years = n_t,
    num_ts = num_ts,
    type = type_sim,
    sigma = sigma_obs,
    loadings_matrix = matrix(nrow = num_ts, ncol = 1,
      rnorm(num_ts * 1, mean_loadings, sigma_loadings)))
  left_out_i <- matrix(sample(seq_len(n_t), round(n_t * 0.1)), nrow = 1)
  for (i in seq(2, num_ts)) {
    left_out_i <- rbind(left_out_i, sample(seq_len(n_t), round(n_t * 0.1)))
  }
  sim_orig <- sim
  for (j in 1:nrow(left_out_i)) {
    sim$y_sim[j, left_out_i] <- NA
  }

  m_rw <- fit_dfa(y = sim$y_sim, iter = stan_iter, chains = 1, scale = scale)
  m_bs <- fit_dfa(
    y = sim$y_sim, iter = stan_iter, chains = 1,
    trend_model = "bs", n_knots = n_knots, scale = scale
  )
  elpd_rw <- get_log_dens(sim_orig$y_sim, m_rw, left_out_i)
  elpd_bs <- get_log_dens(sim_orig$y_sim, m_bs, left_out_i)
  data.frame(epld_rw = elpd_rw, elpd_bs = elpd_bs)
}

# sim_and_fit(type_sim = "RW")

pars <- tidyr::expand_grid(
  num_ts = c(3, 5),
  mean_loadings = c(1), sigma_loadings = c(0.001, 1),
  scale = c("none"), stan_iter = 200,
  sigma_obs = c(0.2, 0.5, 0.75),
  n_knots = c(13), iter = seq_len(20))
set.seed(1)
pars$seed <- sample.int(1e6, nrow(pars))
nrow(pars)

tictoc::tic()
out <- furrr::future_pmap_dfr(pars, sim_and_fit, .progress = TRUE)
tictoc::toc()

saveRDS(out, file = "output/check-sims.rds")

x <- bind_cols(out, pars) %>% filter(!(mean_loadings == 0 & sigma_loadings == 0.001))
x <- filter(x, !(mean_loadings == 0 & sigma_loadings == 1))

make_plot <- function(dat) {
  ggplot(dat,
    aes(sigma_obs, elpd_bs - epld_rw)) +
    geom_point(position = position_jitter(height = 0, width = 0.01)) +
    facet_grid(num_ts~mean_loadings) +
    ggtitle(paste0("Loadings = N(", unique(dat$mean_loadings), ", ",
      unique(dat$sigma_loadings), ")")) +
    geom_hline(yintercept = 0, lty = 2)
}

g <- x %>%
  # filter(sigma_obs > 0.21) %>%
  group_by(mean_loadings, sigma_loadings) %>%
  group_split() %>%
  purrr::map(make_plot)
cowplot::plot_grid(plotlist = g, ncol = 2)

# x_rw <- furrr::future_map_dfr(1:(1*8), ~ sim_and_fit(type_sim = "RW"), .id = "iter")
# plot(x_rw$epld_rw - x_rw$elpd_bs);abline(h = 0, lty = 2)
#
# x_arma <- furrr::future_map_dfr(1:8, ~ sim_and_fit(type_sim = "ARMA"), .id = "iter")
# plot(x_arma$epld_rw - x_arma$elpd_bs, ylim = c(-1, 20));abline(h = 0, lty = 2)
