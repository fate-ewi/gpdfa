# posterior predictive checking with splines as model drawn from:

library(bayesdfa)
library(dplyr)
library(ggplot2)
source("code/sim_dfa2.R")

set.seed(0)
n_t <- 25
num_ts <- 3
sim <- sim_dfa2(num_trends = 1, num_years = n_t, num_ts = num_ts, sigma = 0.6,
  type = "RW", # same but with this arg
  loadings_matrix = matrix(nrow = num_ts, ncol = 1,
    rnorm(num_ts * 1, 1, 0.4)))
matplot(t(sim$y_sim), type = "l")

# m_rw <- fit_dfa(y = post_sim1, iter = 300, chains = 1)
suppressWarnings({
  m_bs <- fit_dfa(
    y = sim$y_sim, iter = 200, chains = 1, scale = "none",
    trend_model = "bs", n_knots = 10, par_list = "all", seed = 29439,
  )
})

# extract spline DFA posterior and re-fit to it with both models:
post <- rstan::extract(m_bs$model)
mean(post$sigma)
sd(post$sigma)
post_sim_true <- post$pred[1,,]
matplot(t(post_sim_true), type = "l", lty = 1)
post_sim_obs <- post_sim_true +
  matrix(rnorm(length(as.vector(post_sim_true)), 0,
    post$sigma[1,1]), ncol = ncol(post$pred[1,,]))
matplot(t(post_sim_obs), type = "l", lty = 1)

sim <- list()
sim$y_sim <- post_sim_obs

m_rw <- fit_dfa(y = sim$y_sim, iter = 300, chains = 1, scale = "none", seed = 1229)
m_bs <- fit_dfa(
  y = sim$y_sim, iter = 300, chains = 1, seed = 19229,
  trend_model = "bs", n_knots = 10, scale = "none"
)

mean(rstan::extract(m_rw$model)$sigma)
mean(rstan::extract(m_bs$model)$sigma)

loo::loo_compare(list(rw = loo::loo(m_rw), bs = loo::loo(m_bs)))

n_t <- ncol(sim$y_sim)

r_rw <- rotate_trends(m_rw, invert = TRUE)
plot_trends(r_rw) +
  geom_line(data = data.frame(sim = post$x[1,,], Time = 1:n_t), aes(Time, sim))

r_bs <- rotate_trends(m_bs, invert = TRUE)
plot_trends(r_bs) +
  geom_line(data = data.frame(sim = post$x[1,,], Time = 1:n_t), aes(Time, sim))

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
sim_df <- post_sim_true %>%
  reshape2::melt() %>%
  rename(time = Var2, time_series = Var1)

# obs values:
sim_df_obs <- sim$y_sim %>%
  reshape2::melt() %>%
  rename(time = Var2, time_series = Var1)

sim_df_orig <- sim$y_sim %>%
  reshape2::melt() %>%
  rename(time = Var2, time_series = Var1)

# missing <- sim_df_obs %>% filter(is.na(value)) %>%
#   select(-value) %>%
#   left_join(sim_df_orig)

mutate(p_df_rw, type = "RW fitted") %>%
  bind_rows(mutate(p_df_bs, type = "BS fitted")) %>%
  # bind_rows(mutate(sim_df_obs, type = "Sim")) %>%
  bind_rows(mutate(sim_df, type = "Sim true")) %>%
  ggplot(aes(time, value, colour = type)) +
  geom_line() +
  # geom_point(data = mutate(missing, type = "Sim")) +
  facet_wrap(~time_series, scale = "free_y")

all <- left_join(
  rename(p_df_rw, rw_fitted = value),
  rename(p_df_bs, bs_fitted = value)) %>%
  left_join(rename(sim_df, true_value = value)) %>%
  left_join(rename(sim_df_obs, obs_value = value))

plot(all$bs_fitted, all$true_value)
plot(all$rw_fitted, all$true_value)

cor(all$bs_fitted, all$true_value)
cor(all$rw_fitted, all$true_value)
