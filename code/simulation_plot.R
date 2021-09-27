library(dplyr)
library(ggplot2)
library(future)
library(viridis)

bs <- readRDS("output/check-sims-bs.rds")
bs$est_trend_model = paste0("bs-",bs$n_knots_fit)
ps <- readRDS("output/check-sims-ps.rds")
gp <- readRDS("output/check-sims-gp.rds")
x = rbind(bs,ps,gp)


x$est_trend_model = as.factor(x$est_trend_model)
x$sigma_obs = as.factor(x$sigma_obs)
x$type_sim = factor(x$type_sim)
levels(x$type_sim) = c("B-spline","Random walk")
# reorder levels
x$est_trend_model = factor(x$est_trend_model,
      levels = c("bs-7","bs-13","bs-20","ps","gp"))

g = ggplot(x,
           aes(est_trend_model, elpd_bs - epld_rw, col=sigma_obs)) +
  geom_boxplot(alpha=0.4, outlier.color = NA) +
  facet_wrap(~type_sim) +
  coord_cartesian(ylim=c(-17,7)) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() +
  xlab("Knots") +
  ylab(expression(paste("ELPD - ",ELPD[RW]))) +
  scale_fill_viridis(discrete=TRUE, end=0.8) +
  scale_color_viridis(discrete=TRUE, end=0.8) +
  theme(strip.background =element_rect(fill="white")) +
  scale_x_discrete(labels = c('gp'=expression(GP[20]),
                              'ps' = expression(PS[13]),
                              'bs-7'=expression(BS[7]),
                              'bs-13'=expression(BS[13]),
                              'bs-20'=expression(BS[20]))) +
  labs(col = expression(sigma[obs]))
