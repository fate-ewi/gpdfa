library(bayesdfa)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(rerddap)
library(loo)
library(rstan)
options(mc.cores = parallel::detectCores())
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

saveRDS(x, "output/calcofi_data.rds")
x = readRDS("output/calcofi_data.rds")
#x$ts = as.numeric(as.factor(x$ts))
#x$time = x$time - min(x$time) + 1
mcmc_iter = 2000
mcmc_chains = 4
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

#saveRDS(m,"calcofi_models.rds")
m = readRDS("output/calcofi_models.rds")

