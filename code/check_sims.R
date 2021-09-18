library(bayesdfa)
library(dplyr)
set.seed(123)
s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 4)
m <- fit_dfa(y = s$y_sim, iter = 400, chains = 1)
matplot(t(s$y_sim), type = "l")

pred <- predicted(m)
p <- rstan::extract(m$model)
matplot(t(post_sim1), type = "l")

post_sim1 <- post_sim1 + rnorm(length(as.vector(post_sim1)), 0, p$sigma[1])
matplot(t(post_sim1), type = "l")

m_rw <- fit_dfa(y = post_sim1, iter = 300, chains = 1)
m_bs <- fit_dfa(
  y = post_sim1, iter = 300, chains = 1,
  trend_model = "bs", n_knots = 10
)

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

get_log_dens(post_sim1, m_rw)
get_log_dens(post_sim1, m_bs)

# repeat!

sim_and_fit <- function(fitted_model, draw_i) {
  pred <- predicted(fitted_model)
  p <- rstan::extract(fitted_model$model)
  post_sim <- t(pred[draw_i, 1L, , ])
  set.seed(draw_i)
  post_sim <- post_sim +
    rnorm(length(as.vector(post_sim)), 0, p$sigma[draw_i])
  m_rw <- fit_dfa(y = post_sim1, iter = 300, chains = 1)
  m_bs <- fit_dfa(
    y = post_sim1, iter = 300, chains = 1,
    trend_model = "bs", n_knots = 10
  )
  data.frame(
    i = draw_i, elpd_rw = get_log_dens(post_sim, m_rw),
    elpd_bs = get_log_dens(post_sim, m_bs)
  )
}

ret <- sim_and_fit(m, 1)
ret <- sim_and_fit(m, 2)

library(future)
plan(multisession)
out <- furrr::future_map_dfr(
  seq_len(20),
  function(x) sim_and_fit(fitted_model = m, draw_i = x)
)

library(ggplot2)
ggplot(out, aes(x = elpd_bs - elpd_rw)) +
  geom_histogram()

# a version with pure simulated data:

sim_and_fit2 <- function(draw_i) {
  cat(draw_i, "\n")
  set.seed(draw_i * 1929)
  pred <- sim_dfa(num_trends = 1, num_years = 30, num_ts = 5, sigma = 0.1)
  m_rw <- fit_dfa(y = pred$y_sim, iter = 200, chains = 1)
  m_bs <- fit_dfa(
    y = pred$y_sim, iter = 200, chains = 1,
    trend_model = "bs", n_knots = 10
  )
  data.frame(
    i = draw_i, elpd_rw = get_log_dens(post_sim, m_rw),
    elpd_bs = get_log_dens(post_sim, m_bs)
  )
}

out2 <- purrr::map_dfr(
  seq_len(8),
  function(x) sim_and_fit2(draw_i = x)
)
ggplot(out2, aes(x = elpd_bs - elpd_rw)) +
  geom_histogram()

# !? look at some?
library(ggplot2)
library(dplyr)
set.seed(5 * 1929)
sim <- sim_dfa(num_trends = 1, num_years = 30, num_ts = 6, sigma = 0.5)
m_rw <- fit_dfa(y = sim$y_sim, iter = 300, chains = 1, scale = "none")
m_bs <- fit_dfa(
  y = sim$y_sim, iter = 300, chains = 1,
  trend_model = "bs", n_knots = 12, scale = "none"
)

r_rw <- rotate_trends(m_rw, invert = TRUE)
plot_trends(r_rw) +
  geom_line(data = data.frame(sim = sim$x[1, ], Time = 1:30), aes(Time, sim))

r_bs <- rotate_trends(m_bs, invert = TRUE)
plot_trends(r_bs) +
  geom_line(data = data.frame(sim = sim$x[1, ], Time = 1:30), aes(Time, sim))

p_rw <- predicted(m_rw)
p_rw <- p_rw[, 1, , ]
p_bs <- predicted(m_bs)
p_bs <- p_bs[, 1, , ]
p_df_rw <- apply(p_rw, c(2, 3), mean) %>%
  reshape2::melt() %>%
  rename(time = Var1, time_series = Var2)
p_df_bs <- apply(p_bs, c(2, 3), mean) %>%
  reshape2::melt() %>%
  rename(time = Var1, time_series = Var2)

# true values:
sim_df <- sim$pred %>%
  reshape2::melt() %>%
  rename(time = Var2, time_series = Var1)

# obs values:
sim_df_obs <- sim$y_sim %>%
  reshape2::melt() %>%
  rename(time = Var2, time_series = Var1)

mutate(p_df_rw, type = "RW fitted") %>%
  bind_rows(mutate(p_df_bs, type = "BS fitted")) %>%
  bind_rows(mutate(sim_df, type = "RW sim")) %>%
  ggplot(aes(time, value, colour = type)) +
  geom_line() +
  facet_wrap(~time_series, scale = "free_y")

# true:
log_sum_exp(dnorm(sim_df$value, p_df_rw$value,
  sd = mean(rstan::extract(m_rw$model)$sigma), log = TRUE
))
log_sum_exp(dnorm(sim_df$value, p_df_bs$value,
  sd = mean(rstan::extract(m_bs$model)$sigma), log = TRUE
))

# obs:
log_sum_exp(dnorm(sim_df_obs$value, p_df_rw$value,
  sd = mean(rstan::extract(m_rw$model)$sigma), log = TRUE
))
log_sum_exp(dnorm(sim_df_obs$value, p_df_bs$value,
  sd = mean(rstan::extract(m_bs$model)$sigma), log = TRUE
))

get_log_dens(sim$y_sim, m_rw)
get_log_dens(sim$y_sim, m_bs)
