library(dplyr)
library(ggplot2)
library(viridis)

sim_results = dplyr::filter(d, phi!=0,theta!=0, sigma!=0.2)
sim_results$model[which(sim_results$model=="rw")] = "RW"
sim_results$model[which(sim_results$model=="gp-7")] = "GP7"
sim_results$model[which(sim_results$model=="gp-13")] = "GP13"
sim_results$model[which(sim_results$model=="gp-20")] = "GP20"
sim_results$model[which(sim_results$model=="spline-7")] = "BS7"
sim_results$model[which(sim_results$model=="spline-13")] = "BS13"
sim_results$model[which(sim_results$model=="spline-20")] = "BS20"
sim_results$model = factor(sim_results$model,
                           levels = c("RW","GP7","GP13","GP20",
                                      "BS7","BS13","BS20"))
sim_results$sigma = as.factor(sim_results$sigma)

# plot of bias
p1 = dplyr::group_by(sim_results, sim, model) %>%
  dplyr::summarize(mean_bias = mean(estimate - obs), sigma=sigma[1]) %>%
  dplyr::group_by(sim) %>%
  dplyr::mutate(rel_bias = mean_bias - mean_bias[which(model=="RW")]) %>%
  dplyr::filter(model!="RW") %>%
  ggplot(aes(x = model, y = rel_bias, fill=sigma)) +
  geom_boxplot(outlier.color = NA) +
  theme_bw() +
  xlab("Model") +
  ylab(expression(paste("Bias - ",Bias[RW]))) +
  scale_x_discrete(labels = c('GP7'=expression(GP[7]),
                              'GP13'=expression(GP[13]),
                              'GP20'=expression(GP[20]),
                              'BS7'=expression(BS[7]),
                              'BS13'=expression(BS[13]),
                              'BS20'=expression(BS[20]))) +
  scale_fill_viridis(discrete = TRUE, end=0.8)

p2 = dplyr::group_by(sim_results, sim, model) %>%
  dplyr::summarize(mean_se = mean(se), sigma=sigma[1]) %>%
  dplyr::group_by(sim) %>%
  dplyr::mutate(rel_se = mean_se / mean_se[which(model=="RW")]) %>%
  dplyr::filter(model!="RW") %>%
ggplot(aes(x = model, y = rel_se, fill=sigma)) +
  geom_boxplot(outlier.color = NA) +
  theme_bw() +
  xlab("Model") +
  ylab(expression(paste("SE / ",SE[RW]))) +
  scale_x_discrete(labels = c('GP7'=expression(GP[7]),
                              'GP13'=expression(GP[13]),
                              'GP20'=expression(GP[20]),
                              'BS7'=expression(BS[7]),
                              'BS13'=expression(BS[13]),
                              'BS20'=expression(BS[20]))) +
  scale_fill_viridis(discrete = TRUE, end=0.8)

# plot ELPD
p3 = dplyr::group_by(sim_results, sim, model) %>%
  dplyr::summarize(tot_elpd = sum(elpd), sigma=sigma[1]) %>%
  dplyr::group_by(sim) %>%
  dplyr::mutate(elpd_diff = tot_elpd - tot_elpd[which(model=="RW")]) %>%
  dplyr::filter(model!="RW") %>%
  ggplot(aes(x = model, y = elpd_diff, fill=sigma)) +
  geom_boxplot(outlier.color = NA) +
  theme_bw() +
  xlab("Model") +
  ylab(expression(paste("ELPD - ",ELPD[RW]))) +
  scale_x_discrete(labels = c('GP7'=expression(GP[7]),
                              'GP13'=expression(GP[13]),
                              'GP20'=expression(GP[20]),
                              'BS7'=expression(BS[7]),
                              'BS13'=expression(BS[13]),
                              'BS20'=expression(BS[20]))) +
  scale_fill_viridis(discrete = TRUE, end=0.8)
