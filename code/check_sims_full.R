# this one leaves out data points and tests on those
remotes::install_github("fate-ewi/bayesdfa","spline-sim")
library(bayesdfa)
library(dplyr)
library(ggplot2)
library(future)
library(viridis)
plan(multisession, workers = parallel::detectCores() / 2)

log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

get_log_dens <- function(sim_dat, model, out_of_sample_i) {
  e <- rstan::extract(model$model)
  pred <- predicted(model)
  log_density <- 0
  for (j in 1:nrow(out_of_sample_i)) {
    for (i in 1:ncol(out_of_sample_i)) {
        log_density <- log_density +
          log_sum_exp(dnorm(
            sim_dat[j, out_of_sample_i[j, i]],
            mean = pred[, 1L, out_of_sample_i[j, i], j],
            sd = e$sigma,
            log = TRUE
          ))
    }
  }
  n_samples <- nrow(e$sigma)
  log_density - log(n_samples * length(out_of_sample_i))
}

sim_and_fit <- function(sigma_obs = 0.25, sigma_loadings = 0.1,
  mean_loadings = 1, scale = "none",
  type_sim = c("rw", "bs"), num_ts = 3, n_t = 20, stan_iter = 300,
  n_knots_fit = 7, n_knots_sim = 7, est_trend_model = "bs",
  seed = sample.int(.Machine$integer.max, 1L), ...) {
  type_sim <- match.arg(type_sim)
  num_trends <- 1
  set.seed(seed)
  sim <- sim_dfa(
    num_trends = num_trends,
    num_years = n_t,
    num_ts = num_ts,
    trend_model = type_sim,
    spline_weights = matrix(ncol = n_knots_sim, nrow = num_trends,
      data = rnorm(n_knots_sim * num_trends, 0, 1)),
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
    trend_model = est_trend_model, n_knots = n_knots_fit, scale = scale
  )

  max_rhat_rw <- get_rhat(m_rw$model)
  max_rhat_bs <- get_rhat(m_bs$model)

  elpd_rw <- get_log_dens(sim_orig$y_sim, m_rw, left_out_i)
  elpd_bs <- get_log_dens(sim_orig$y_sim, m_bs, left_out_i)

  data.frame(
    epld_rw = elpd_rw, elpd_bs = elpd_bs,
    max_rhat_rw = max_rhat_rw, max_rhat_bs = max_rhat_bs
  )
}

get_rhat <- function(x) {
  s <- rstan::summary(x)$summary
  par_omit <- grep("log_lik", row.names(s))
  par_omit <- union(par_omit, grep("lp__", row.names(s)))
  max(s[-par_omit,"Rhat"])
}

pars <- tidyr::expand_grid(
  num_ts = c(3),
  mean_loadings = 1,
  sigma_loadings = 0.1,
  scale = c("none"),
  stan_iter = 1000,
  sigma_obs = c(0.25, 0.5, 1),
  type_sim = c("rw"),
  est_trend_model = "bs",
  n_knots_sim = 7,
  n_knots_fit = c(7,13,20),
  iter = seq_len(100)
)
pars2 <- pars
pars2$type_sim <- "bs"
pars <- bind_rows(pars, pars2)

set.seed(1)
pars$seed <- sample.int(1e6, nrow(pars))
nrow(pars)

tictoc::tic()
# out <- purrr::pmap_dfr(pars[1,], sim_and_fit)
out <- furrr::future_pmap_dfr(pars, sim_and_fit, .progress = TRUE)
tictoc::toc()

saveRDS(out, file = "output/check-sims-bs.rds")


pars <- tidyr::expand_grid(
  num_ts = c(3),
  mean_loadings = 1,
  sigma_loadings = 0.1,
  scale = c("none"),
  stan_iter = 1000,
  sigma_obs = c(0.25, 0.5, 1),
  type_sim = c("rw"),
  est_trend_model = "ps",
  n_knots_sim = 7,
  n_knots_fit = c(13),
  iter = seq_len(100)
)
pars2 <- pars
pars2$type_sim <- "bs"
pars <- bind_rows(pars, pars2)

set.seed(1)
pars$seed <- sample.int(1e6, nrow(pars))
nrow(pars)

tictoc::tic()
# out <- purrr::pmap_dfr(pars[1,], sim_and_fit)
out <- furrr::future_pmap_dfr(pars, sim_and_fit, .progress = TRUE)
tictoc::toc()

x <- bind_cols(out, pars)
saveRDS(x, file = "output/check-sims-ps.rds")

pars <- tidyr::expand_grid(
  num_ts = c(3),
  mean_loadings = 1,
  sigma_loadings = 0.1,
  scale = c("none"),
  stan_iter = 1000,
  sigma_obs = c(0.25, 0.5, 1),
  type_sim = c("rw"),
  est_trend_model = "gp",
  n_knots_sim = 7,
  n_knots_fit = c(20),
  iter = seq_len(100)
)
pars2 <- pars
pars2$type_sim <- "bs"
pars <- bind_rows(pars, pars2)

set.seed(1)
pars$seed <- sample.int(1e6, nrow(pars))
nrow(pars)

tictoc::tic()
# out <- purrr::pmap_dfr(pars[1,], sim_and_fit)
out <- furrr::future_pmap_dfr(pars, sim_and_fit, .progress = TRUE)
tictoc::toc()

x <- bind_cols(out, pars)
saveRDS(x, file = "output/check-sims-gp.rds")

# make_plot <- function(dat) {
#   ggplot(dat,
#     aes(sigma_obs, elpd_bs - epld_rw, group = sigma_obs)) +
#     geom_boxplot() +
#     geom_point(position = position_jitter(height = 0, width = 0.01), alpha=0.5) +
#     facet_grid(num_ts~mean_loadings) +
#     ggtitle(paste0("Loadings = N(", unique(dat$mean_loadings), ", ",
#       unique(dat$sigma_loadings), ")")) +
#     geom_hline(yintercept = 0, lty = 2)
# }
#
# g <- x %>%
#   # filter(sigma_obs > 0.21) %>%
#   # filter(sigma_loadings > 0.1) %>%
#   group_by(mean_loadings, sigma_loadings) %>%
#   group_split() %>%
#   purrr::map(make_plot)
# cowplot::plot_grid(plotlist = g, ncol = 2)

# x_rw <- furrr::future_map_dfr(1:(1*8), ~ sim_and_fit(type_sim = "RW"), .id = "iter")
# plot(x_rw$epld_rw - x_rw$elpd_bs);abline(h = 0, lty = 2)
#
# x_arma <- furrr::future_map_dfr(1:8, ~ sim_and_fit(type_sim = "ARMA"), .id = "iter")
# plot(x_arma$epld_rw - x_arma$elpd_bs, ylim = c(-1, 20));abline(h = 0, lty = 2)

# 0.25, 0.5, 1
# 3 time series
# N(0, 0.1)
# GP + 7 + 14
#
