library(bayesdfa)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(rerddap)
library(loo)
library(rstan)

options(mc.cores = parallel::detectCores())

d <- read.csv("data/port_landings_table2.csv", stringsAsFactors = FALSE)
d <- dplyr::select(d, -Year)
for (i in 1:ncol(d)) {
  d[, i] <- log(as.numeric(d[, i]))
}

x <- data.frame(
  "obs" = c(as.matrix(d)), "time" = rep(1:nrow(d), ncol(d)),
  "ts" = sort(rep(1:ncol(d), nrow(d)))
)

# x$ts = as.numeric(as.factor(x$ts))
# x$time = x$time - min(x$time) + 1
mcmc_iter <- 4000
mcmc_chains <- 3
m <- list()
num_trends <- 2
fold_ids <- x$time - min(x$time) + 1


bs = fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = 1,
             num_trends = num_trends, trend_model = "spline", n_knots = 6, sample = TRUE)
gp = fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, chains = 1,
             num_trends = num_trends, trend_model = "gp", n_knots = 6, sample = TRUE)


rbs = rotate_trends(bs)
rgp = rotate_trends(gp)

rotated <- rbs
df <- data.frame(
  x = c(t(rotated$trends_mean)),
  lo = c(t(rotated$trends_lower)),
  hi = c(t(rotated$trends_upper)),
  trend = paste0("Trend ", sort(rep(seq_len(2), 39))),
  time = rep(1:39, 2),
  model = "bs"
)
rotated <- rgp
df2 <- data.frame(
  x = c(t(rotated$trends_mean)),
  lo = c(t(rotated$trends_lower)),
  hi = c(t(rotated$trends_upper)),
  trend = paste0("Trend ", sort(rep(seq_len(2), 39))),
  time = rep(1:39, 2),
  model = "gp"
)
df = rbind(df,df2)

p1 = ggplot(df, aes(time, x,group=model,col=model,fill=model)) +
  geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.4) +
  geom_line() + facet_wrap(~trend) +
  xlab("Time") + ylab("Trend") + theme_bw()


# compare loadings (Z)
df <- data.frame(
  x = c(apply(rbs$Z_rot,c(2,3),mean)),
  lo = c(apply(rbs$Z_rot,c(2,3),quantile,0.025)),
  hi = c(apply(rbs$Z_rot,c(2,3),quantile,0.975)),
  trend = paste0("Trend ", sort(rep(seq_len(2), 13))),
  ts = rep(1:13, 2),
  model = "bs"
)
df2 <- data.frame(
  x = c(apply(rgp$Z_rot,c(2,3),mean)),
  lo = c(apply(rgp$Z_rot,c(2,3),quantile,0.025)),
  hi = c(apply(rgp$Z_rot,c(2,3),quantile,0.975)),
  trend = paste0("Trend ", sort(rep(seq_len(2), 13))),
  ts = rep(1:13, 2),
  model = "gp"
)
df = rbind(df,df2)
p2 = ggplot(df, aes(ts, x,group=model,col=model,fill=model)) +
  geom_point(position=position_dodge(0.5)) +
  geom_linerange(aes(ymin=lo,ymax=hi),position=position_dodge(0.5))+
  facet_wrap(~trend) + xlab("Time series") +
  ylab("Rotated loadings, with 95% CIs") +
  theme_bw()

# look at predictions from these models for a given time series
pred_bs = predicted(bs)
pred_gp = predicted(gp)

df_gp = data.frame("time" = rep(1:39,13),
                "ts"=sort(rep(1:13,39)),
                "x"=c(apply(pred_gp,c(3,4),mean)),
                "lo"=c(apply(pred_gp,c(3,4),quantile,0.025)),
                "hi"=c(apply(pred_gp,c(3,4),quantile,0.975)),
                "model"="gp")
df_bs = data.frame("time" = rep(1:39,13),
                   "ts"=sort(rep(1:13,39)),
                   "x"=c(apply(pred_bs,c(3,4),mean)),
                   "lo"=c(apply(pred_bs,c(3,4),quantile,0.025)),
                   "hi"=c(apply(pred_bs,c(3,4),quantile,0.975)),
                   "model"="bs")
df = rbind(df_gp,df_bs)
df = dplyr::left_join(df,x)
ggplot(df, aes(time,x,group=model,col=model,fill=model)) +
  geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.4) +
  geom_line() +
  facet_wrap(~ts,scale="free_y") + theme_bw() +
  geom_point(aes(time,(obs - mean(obs))/sd(obs)),col="black")


# drop out 2nd to last point

bs = fit_dfa(y = x, data_shape = "long", iter = 2000, chains = 1,
             num_trends = 2, trend_model = "spline", n_knots = 6, sample = FALSE)
gp = fit_dfa(y = x, data_shape = "long", iter = 2000, chains = 1,
             num_trends = 2, trend_model = "gp", n_knots = 6, sample = FALSE)


rbs = rotate_trends(bs)
rgp = rotate_trends(gp)
