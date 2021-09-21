library(bayesdfa)
library(dplyr)
library(ggplot2)
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

set.seed(123)

# Pick one:
# TYPE <- "ARMA"
TYPE <- "RW"

sim <- sim_dfa2(num_trends = 1, num_years = 30, num_ts = 6, sigma = 0.3,
  type = TYPE,
  loadings_matrix = matrix(nrow = 6, ncol = 1,
      rnorm(6 * 1, 1, 0.2)))
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

all <- left_join(
  rename(p_df_rw, rw_fitted = value),
  rename(p_df_bs, bs_fitted = value)) %>%
  left_join(rename(sim_df, true_value = value)) %>%
  left_join(rename(sim_df_obs, obs_value = value))

plot(all$bs_fitted, all$true_value)
plot(all$rw_fitted, all$true_value)

# true; hack the sigma:
log_sum_exp(dnorm(all$rw_fitted, all$true_value,
  sd = mean(rstan::extract(m_rw$model)$sigma), log = TRUE
))
log_sum_exp(dnorm(all$bs_fitted, all$true_value,
  sd = mean(rstan::extract(m_bs$model)$sigma), log = TRUE
))

# obs; hacked sigma:
log_sum_exp(dnorm(all$rw_fitted, all$obs_value,
  sd = mean(rstan::extract(m_rw$model)$sigma), log = TRUE
))
log_sum_exp(dnorm(all$bs_fitted, all$obs_value,
  sd = mean(rstan::extract(m_bs$model)$sigma), log = TRUE
))

# actual ELPD:
get_log_dens(sim$y_sim, m_rw)
get_log_dens(sim$y_sim, m_bs)

# check what get_log_dens() is doing:
e <- rstan::extract(m_rw$model)
pred <- predicted(m_rw)
plot(sim$y_sim[1,], apply(pred[, 1L, , 1L], 2, mean))
plot(sim$y_sim[2,], apply(pred[, 1L, , 2L], 2, mean))
plot(sim$y_sim[3,], apply(pred[, 1L, , 3L], 2, mean))
# looks good
