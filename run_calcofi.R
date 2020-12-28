library(bayesdfa)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(rerddap)
library(loo)
library(rstan)
#http://calcofi.org/publications/calcofireports/v41/Vol_41_Moser_etal.pdf
#Species included in the study were bocaccio
#(S. paucispinis), cowcod (S. levis), shortbelly rockfish
#(S. jordani), aurora rockfish (S. uuvoya), and splitnose rockfish (S. diplopyou)

calcofi_erddap <- "ERDDAP	Species
erdCalCOFIlrvcntSDtoSI	Sebastes paucispinis
erdCalCOFIlrvcntSDtoSI	Sebastes jordani
erdCalCOFIlrvcntSDtoSI	Sebastes aurora"

calcofi = read.table(textConnection(calcofi_erddap), header=TRUE, sep="\t")

for(i in 1:nrow(calcofi)) {

  out = info(as.character(calcofi$ERDDAP[i]))
  station1 = tabledap(out, fields = c("station", "line", "latitude",
    "longitude","time","scientific_name", "larvae_10m2"))
  station1$date = lubridate::as_date(station1$time)
  station1$year = lubridate::year(station1$date)
  station1$month = lubridate::month(station1$date)
  station1$yday = lubridate::yday(station1$date)

  # filter out recent years with consistent sampling
  station1 = dplyr::filter(station1, year >= 1985)

  # time windows https://calcofi.org/cruises.html
  station1$quarter = NA
  station1$quarter[which(station1$yday %in% 2:52)] = 1
  station1$quarter[which(station1$yday %in% 90:141)] = 2
  station1$quarter[which(station1$yday %in% 181:233)] = 3
  station1$quarter[which(station1$yday %in% 273:325)] = 4
  station1 = dplyr::filter(station1, !is.na(quarter))

  #https://calcofi.org/field-work/station-positions.html
  stations = read.csv("data/CalCOFIStationOrder.csv")
  stations = dplyr::rename(stations, line = Line, station = Station)
  station1 = dplyr::left_join(station1, stations[,c("line","station","StaType")])
  station1 = dplyr::filter(station1, StaType=="ROS")
  station1$larvae_10m2 = as.numeric(station1$larvae_10m2)
  dat = dplyr::filter(station1, scientific_name %in% calcofi$Species[i])

  if(i==1) {
    df_all = dat
  } else {
    df_all = rbind(df_all, dat)
  }

}

x = group_by(df_all, scientific_name, year, quarter) %>%
  dplyr::summarize(mu = mean(larvae_10m2, na.rm=T)) %>%
  dplyr::filter(year >= 1985, quarter==2)

x = dplyr::rename(x, ts = scientific_name, time = year, obs = mu) %>%
  dplyr::filter(ts %in% c("Sebastes aurora","Sebastes jordani", "Sebastes paucispinis"))
x = dplyr::select(x, -quarter)
x$ts = as.numeric(as.factor((x$ts)))

x$obs[which(x$obs==0)]=NA
x$obs = log(x$obs)

saveRDS(x, "calcofi_data.rds")
x = readRDS("calcofi_data.rds")
#x$ts = as.numeric(as.factor(x$ts))
#x$time = x$time - min(x$time) + 1
mcmc_iter = 3000
mcmc_chains = 3
m = list()
m[[1]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1)

# fit models with varying numbers of knots in a b-spline
m[[2]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline",n_knots=6)
m[[3]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline",n_knots=12)
m[[4]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline",n_knots=18)
m[[5]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline",n_knots=24)
m[[6]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "spline",n_knots=30)

# fit models with varying numbers of knots in a gp
m[[7]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp",n_knots=6)
m[[8]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp",n_knots=12)
m[[9]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp",n_knots=18)
m[[10]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp",n_knots=24)
m[[11]] = fit_dfa(y = x, data_shape="long", iter = mcmc_iter, chains = mcmc_chains, num_trends = 1, trend_model = "gp",n_knots=30)

saveRDS(m,"calcofi_models.rds")
m = readRDS("calcofi_models.rds")

# calculate summary table
loos = data.frame("loo"=rep(NA,5), "loo_se"=NA, Model="B-spline","knots"=seq(6,30,6))
for(i in 2:6) {
  loos$loo[i-1] = loo(m[[i]])$estimate["looic","Estimate"]
  loos$loo_se[i-1] = loo(m[[i]])$estimate["looic","SE"]
}

loos_gp = data.frame("loo"=rep(NA,5), "loo_se"=NA, Model="GP","knots"=seq(6,30,6))
for(i in 7:11) {
  loos_gp$loo[i-6] = loo(m[[i]])$estimate["looic","Estimate"]
  loos_gp$loo_se[i-6] = loo(m[[i]])$estimate["looic","SE"]
}
loos = rbind(loos, loos_gp)

loos_fixed = data.frame("loo"=rep(NA,1), "loo_se"=NA, Model="RW","knots"=seq(32))
loos_fixed$loo[1] = loo(m[[1]])$estimate["looic","Estimate"]
loos_fixed$loo_se[1] = loo(m[[1]])$estimate["looic","SE"]
loos = rbind(loos, loos_fixed)
loos = loos[1:11,]

write.csv(loos,"loos_table_calcofi.csv")

# use best model -- GP with ~ 12 knots
fit = m[[5]]
r = bayesdfa::rotate_trends(fit)

trend_df = data.frame("time"= seq(min(x$time),max(x$time)),
  "low" = c(r$trends_lower),
  "mean"=c(r$trends_mean),
  "hi"=c(r$trends_upper))

scaled_x = group_by(x, ts) %>%
  dplyr::mutate(obs = (obs-mean(obs))/sd(obs))
scaled_x$Species = c("aurora","shortbelly","bocaccio")[scaled_x$ts]


ggplot(trend_df, aes(time,mean)) +
  geom_ribbon(aes(ymin=low,ymax=hi),alpha=0.2) +
  geom_line() +
  theme_bw() +
  geom_point(data=scaled_x, aes(x=time,y=obs,col=Species),alpha=0.7,size=3) +
  scale_color_manual("Species",
    values = c("aurora" = "#440154FF", "shortbelly" = "#2A788EFF", "bocaccio" = "#7AD151FF")) +
  ylab("Trend estimate") +
  xlab("Year")

loadings = as.data.frame(cbind(apply(r$Z_rot,c(2,3),mean), apply(r$Z_rot,c(2,3),sd)))
names(loadings) = c("Trend_mean","Trend_sd")
names(d)[1] = "Pacific whiting"
names(d)[4] = "Pacific cod"
names(d)[5] = "Misc roundfish"
names(d)[9:13] = c("Dover sole","English sole","Petrale sole","Misc flatfish","Misc groundfish")
loadings$species = names(d)

