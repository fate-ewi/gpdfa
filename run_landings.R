library(bayesdfa)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)

d = read.csv("data/port_landings_table2.csv", stringsAsFactors = FALSE)
d = dplyr::select(d, -Year)
for(i in 1:ncol(d)) {
  d[,i] = log(as.numeric(d[,i]))
}

mcmc_iter = 3000
mcmc_chains = 3
m = list()
m[[1]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 1)
m[[2]] = fit_dfa(y = t(d), iter = 1000, chains = 1, num_trends = 2)

# fit models with varying numbers of knots in a b-spline
m[[3]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=6)
m[[4]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=12)
m[[5]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=18)
m[[6]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=24)
m[[7]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=30)

# fit models with varying numbers of knots in a b-spline
m[[8]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=6)
m[[9]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=12)
m[[10]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=18)
m[[11]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=24)
m[[12]] = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "gp",n_knots=30)

fit = fit_dfa(y = t(d), iter = mcmc_iter, chains = mcmc_chains, num_trends = 2, trend_model = "spline",n_knots=18)
saveRDS(fit,"landings_spline_18.rds")

loos = data.frame("loo"=rep(NA,5), "loo_se"=NA, Model="B-spline","knots"=seq(6,30,6))
for(i in 3:7) {
  loos$loo[i-2] = loo(m[[i]])$estimate["looic","Estimate"]
  loos$loo_se[i-2] = loo(m[[i]])$estimate["looic","SE"]
}

loos_gp = data.frame("loo"=rep(NA,5), "loo_se"=NA, Model="GP","knots"=seq(6,30,6))
for(i in 8:12) {
  loos_gp$loo[i-7] = loo(m[[i]])$estimate["looic","Estimate"]
  loos_gp$loo_se[i-7] = loo(m[[i]])$estimate["looic","SE"]
}
loos = rbind(loos, loos_gp)

loos_fixed = data.frame("loo"=rep(NA,1), "loo_se"=NA, Model="RW","knots"=seq(32))
loos_fixed$loo[1] = loo(m[[2]])$estimate["looic","Estimate"]
loos_fixed$loo_se[1] = loo(m[[2]])$estimate["looic","SE"]

loos = rbind(loos, loos_fixed)
loos = loos[1:11,]
#write.csv(loos,"loos_table.csv")

loos = read.csv("loos_table.csv")
pdf("plots/Figure_4.pdf")
# plot effect of knots
ggplot(loos, aes(knots, loo, group = Model, col=Model)) +
  geom_pointrange(aes(ymin=loo-loo_se, ymax=loo+loo_se),position = position_dodge(width=2), size=1,alpha=0.7) +
  theme_bw() +
  xlab("Knots") +
  ylab("Leave One Out Information Criterion (+/- SE)") +
  scale_color_viridis(discrete=TRUE,end=0.8)
dev.off()

# make plots for best model
m = readRDS("landings_spline_6.rds")
r = rotate_trends(m)

loadings = as.data.frame(cbind(apply(r$Z_rot,c(2,3),mean), apply(r$Z_rot,c(2,3),sd)))
names(loadings) = c("Trend1_mean","Trend2_mean","Trend1_sd","Trend2_sd")
names(d)[1] = "Pacific whiting"
names(d)[4] = "Pacific cod"
names(d)[5] = "Misc roundfish"
names(d)[8] = "Arrowtooth flounder"
names(d)[9:13] = c("Dover sole","English sole","Petrale sole","Misc flatfish","Misc groundfish")
loadings$species = names(d)



pdf("plots/Figure_3.pdf")
ggplot(dplyr::filter(loadings, species!="Misc roundfish",species!="Misc flatfish",species!="Misc groundfish"),
  aes(Trend1_mean,Trend2_mean,label=species)) +
  geom_point(col="darkblue") +
  geom_errorbar(aes(ymin=Trend2_mean - Trend2_sd,ymax=Trend2_mean + Trend2_sd),col="darkblue",alpha=0.3) +
  geom_errorbarh(aes(xmin=Trend1_mean - Trend1_sd,xmax=Trend1_mean + Trend1_sd),col="darkblue",alpha=0.3) +
  geom_text_repel(size=5,col="dark blue") +
  theme_bw() +
  xlab("Loading on trend 1") +
  ylab("Loading on trend 2")
dev.off()

# plot RW model, GP model, and spline model on top of eachother for predicted values
pred_rw_mean = apply(predicted(m[[2]]), c(3,4), mean)
pred_rw_lower = apply(predicted(m[[2]]), c(3,4), quantile,0.025)
pred_rw_upper = apply(predicted(m[[2]]), c(3,4), quantile,0.975)
pred_sp6_mean = apply(predicted(m[[3]]), c(3,4), mean)
pred_sp6_lower = apply(predicted(m[[3]]), c(3,4), quantile,0.025)
pred_sp6_upper = apply(predicted(m[[3]]), c(3,4), quantile,0.975)
pred_sp18_mean = apply(predicted(m[[5]]), c(3,4), mean)
pred_sp18_lower = apply(predicted(m[[5]]), c(3,4), quantile,0.025)
pred_sp18_upper = apply(predicted(m[[5]]), c(3,4), quantile,0.975)

# grab predictions for rockfish and lingcod
df_1 = data.frame("Species"="Lingcod",
  year = rep(1981:2019,3), Model = c(rep("Random walk",39), rep("B-spline, 6 knots",39), rep("B-spline, 18 knots",39)),
  pred = c(pred_rw_mean[,3], pred_sp6_mean[,3], pred_sp18_mean[,3]),
  lower = c(pred_rw_lower[,3], pred_sp6_lower[,3], pred_sp18_lower[,3]),
  upper = c(pred_rw_upper[,3], pred_sp6_upper[,3], pred_sp18_upper[,3]))
df_2 = data.frame("Species"="Pacific whiting",
  year = rep(1981:2019,3), Model = c(rep("Random walk",39), rep("B-spline, 6 knots",39), rep("B-spline, 18 knots",39)),
  pred = c(pred_rw_mean[,1], pred_sp6_mean[,1], pred_sp18_mean[,1]),
  lower = c(pred_rw_lower[,1], pred_sp6_lower[,1], pred_sp18_lower[,1]),
  upper = c(pred_rw_upper[,1], pred_sp6_upper[,1], pred_sp18_upper[,1]))

#
saveRDS(rbind(df_1,df_2), "predicted_lingcod_whiting.rds")
# bring in raw observations
d = read.csv("data/port_landings_table2.csv", stringsAsFactors = FALSE)
obs = data.frame(year = rep(1981:2019,2),
  y = c(d[["Lingcod"]],d[["P..Whiting"]]),
  "Species"=c(rep("Lingcod",39), rep("Pacific whiting",39)))
obs = dplyr::group_by(obs, Species) %>%
  dplyr::mutate(y = (y-mean(y))/sd(y),
    Model="RW")

pdf("plots/Figure_2.pdf")
ggplot(rbind(df_1,df_2), aes(year, pred, group=Model, col=Model,fill=Model)) +
  geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,col=NA) +
  geom_line() +
  theme_bw() +
  xlab("") +
  ylab("Standardized landings (mt)") +
  scale_fill_viridis(end=0.8, discrete = TRUE) +
  scale_color_viridis(end=0.8, discrete = TRUE) +
  facet_wrap(~ Species) +
  theme(strip.background =element_rect(fill="white")) +
  geom_point(data=obs, aes(year, y),size=2, col="grey30",alpha=0.5)
dev.off()

# plot trend estimates
r = rotate_trends(m[[3]])
p1 = plot_trends(r, years = 1981:2019) + theme_bw() +
  theme(strip.background =element_rect(fill="white")) + ylab("Estimate, 6 knots")
r = rotate_trends(m[[7]])
p2 = plot_trends(r, years = 1981:2019) + theme_bw() +
  theme(strip.background =element_rect(fill="white")) + ylab("Estimate, 30 knots")

# plot raw rotated trends from landings example
pdf("plots/Figure_S2.pdf")
gridExtra::grid.arrange(p1, p2, nrow=2)
dev.off()
