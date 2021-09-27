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

get_log_dens_left_out <- function(sim_dat, model, out_of_sample_i) {
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
  log_density
}

set.seed(1)

# Pick one:
# TYPE <- "ARMA"
TYPE <- "RW"

n_t <- 25
num_ts <- 4
sim <- sim_dfa2(num_trends = 1, num_years = n_t, num_ts = num_ts, sigma = 0.6,
  type = TYPE,
  loadings_matrix = matrix(nrow = num_ts, ncol = 1,
      rnorm(num_ts * 1, 1, 0.3)))
matplot(t(sim$y_sim), type = "l", lty = 2)

left_out_i <- matrix(sample(seq_len(n_t), round(n_t * 0.1)), nrow = 1)
for (i in seq(2, num_ts)) {
  left_out_i <- rbind(left_out_i, sample(seq_len(n_t), round(n_t * 0.1)))
}
sim_orig <- sim
for (j in 1:nrow(left_out_i)) {
  sim$y_sim[j, left_out_i[j,]] <- NA
}

m_rw <- fit_dfa(y = sim$y_sim, iter = 300, chains = 1, scale = "none")
m_bs <- fit_dfa(
  y = sim$y_sim, iter = 300, chains = 1,
  trend_model = "bs", n_knots = 12, scale = "none"
)

r_rw <- rotate_trends(m_rw, invert = TRUE)
plot_trends(r_rw) +
  geom_line(data = data.frame(sim = sim$x[1, ], Time = 1:n_t), aes(Time, sim))

r_bs <- rotate_trends(m_bs, invert = TRUE)
plot_trends(r_bs) +
  geom_line(data = data.frame(sim = sim$x[1, ], Time = 1:n_t), aes(Time, sim))

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

sim_df_orig <- sim_orig$y_sim %>%
  reshape2::melt() %>%
  rename(time = Var2, time_series = Var1)

missing <- sim_df_obs %>% filter(is.na(value)) %>%
  select(-value) %>%
  left_join(sim_df_orig)

mutate(p_df_rw, type = "RW fitted") %>%
  bind_rows(mutate(p_df_bs, type = "BS fitted")) %>%
  bind_rows(mutate(sim_df_obs, type = "Sim")) %>%
  ggplot(aes(time, value, colour = type)) +
  geom_line() +
  geom_point(data = mutate(missing, type = "Sim")) +
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
log_sum_exp(dnorm(all$rw_fitted, all$true_value,
  sd = mean(rstan::extract(m_rw$model)$sigma), log = TRUE
))
log_sum_exp(dnorm(all$bs_fitted, all$true_value,
  sd = mean(rstan::extract(m_bs$model)$sigma), log = TRUE
))

# actual ELPD:
get_log_dens(sim_orig$y_sim, m_rw)
get_log_dens(sim_orig$y_sim, m_bs)

# check what get_log_dens() is doing:
# e <- rstan::extract(m_rw$model)
# pred <- predicted(m_rw)
# plot(sim$y_sim[1,], apply(pred[, 1L, , 1L], 2, mean))
# plot(sim$y_sim[2,], apply(pred[, 1L, , 2L], 2, mean))
# plot(sim$y_sim[3,], apply(pred[, 1L, , 3L], 2, mean))
# looks good

check <- missing %>%
  left_join(rename(p_df_rw, rw_fitted = value)) %>%
  left_join(rename(p_df_bs, bs_fitted = value))

plot(check$rw_fitted, check$value)
plot(check$bs_fitted, check$value)

mean(rstan::extract(m_rw$model)$sigma)
mean(rstan::extract(m_bs$model)$sigma)

cor(check$rw_fitted, check$value)
cor(check$bs_fitted, check$value)

get_log_dens_left_out(sim_orig$y_sim, m_rw, left_out_i)
get_log_dens_left_out(sim_orig$y_sim, m_bs, left_out_i)
