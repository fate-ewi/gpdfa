library(ggplot2)
library(dplyr)
library(viridis)
library(ggrepel)
library(rerddap)
library(loo)
library(rstan)
library(bayesdfa)
options(mc.cores = parallel::detectCores())
# http://calcofi.org/publications/calcofireports/v41/Vol_41_Moser_etal.pdf
# Species included in the study were bocaccio
# (S. paucispinis), cowcod (S. levis), shortbelly rockfish
# (S. jordani), aurora rockfish (S. aurora), and splitnose rockfish (S. diplopyou)

calcofi_erddap <- "ERDDAP	Species
erdCalCOFIlrvcntSDtoSI	Sebastes paucispinis
erdCalCOFIlrvcntSDtoSI	Sebastes jordani
erdCalCOFIlrvcntSDtoSI	Sebastes aurora
erdCalCOFIlrvcntSDtoSI	Sebastes diploproa
erdCalCOFIlrvcntSDtoSI	Sebastes levis
"

calcofi <- read.table(textConnection(calcofi_erddap), header = TRUE, sep = "\t")

for (i in 1:nrow(calcofi)) {
  out <- info(as.character(calcofi$ERDDAP[i]))
  station1 <- tabledap(out, fields = c(
    "station", "line", "latitude",
    "longitude", "time", "scientific_name", "larvae_10m2"
  ))
  station1$date <- lubridate::as_date(station1$time)
  station1$year <- lubridate::year(station1$date)
  station1$month <- lubridate::month(station1$date)
  station1$yday <- lubridate::yday(station1$date)

  # filter out recent years with consistent sampling
  station1 <- dplyr::filter(station1, year >= 1985)

  # time windows https://calcofi.org/cruises.html
  station1$quarter <- NA
  station1$quarter[which(station1$yday %in% 2:52)] <- 1
  station1$quarter[which(station1$yday %in% 90:141)] <- 2
  station1$quarter[which(station1$yday %in% 181:233)] <- 3
  station1$quarter[which(station1$yday %in% 273:325)] <- 4
  station1 <- dplyr::filter(station1, !is.na(quarter))

  # https://calcofi.org/field-work/station-positions.html
  stations <- read.csv("data/CalCOFIStationOrder.csv")
  stations <- dplyr::rename(stations, line = Line, station = Station)
  station1 <- dplyr::left_join(station1, stations[, c("line", "station", "StaType")])
  station1 <- dplyr::filter(station1, StaType == "ROS")
  station1$larvae_10m2 <- as.numeric(station1$larvae_10m2)
  dat <- dplyr::filter(station1, scientific_name %in% calcofi$Species[i])

  if (i == 1) {
    df_all <- dat
  } else {
    df_all <- rbind(df_all, dat)
  }
}

x <- group_by(df_all, scientific_name, year, quarter) %>%
  dplyr::summarize(mu = mean(larvae_10m2, na.rm = T)) %>%
  dplyr::filter(year >= 1985, quarter == 2)

x <- dplyr::rename(x, ts = scientific_name, time = year, obs = mu)
x <- dplyr::select(x, -quarter)
x = dplyr::filter(x, ts%in%c("Sebastes paucispinis","Sebastes jordani","Sebastes aurora"))
x$ts <- as.numeric(as.factor((x$ts)))

x$obs[which(x$obs == 0)] <- NA
x$obs <- log(x$obs)

saveRDS(x, "output/calcofi_data.rds")
x <- readRDS("output/calcofi_data.rds")
# x$ts = as.numeric(as.factor(x$ts))
x$time = x$time - min(x$time) + 1
mcmc_iter <- 4000
mcmc_chains <- 3
mcmc_thin <- 10
m <- list()
loos = list()

m[[1]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin, chains = mcmc_chains, num_trends = 1,par_list="all")
loos[[1]] = loo::loo(m[[1]]$model,moment_match=TRUE)

# fit models with varying numbers of knots in a b-spline
m[[2]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "bs", n_knots = 6,par_list="all")
loos[[2]] = loo::loo(m[[2]]$model,moment_match=TRUE)

m[[3]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "bs", n_knots = 12,par_list="all")
loos[[3]] = loo::loo(m[[3]]$model,moment_match=TRUE)

m[[4]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "bs", n_knots = 18,par_list="all")
loos[[4]] = loo::loo(m[[4]]$model,moment_match=TRUE)

m[[5]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "bs", n_knots = 24,par_list="all")
loos[[5]] = loo::loo(m[[5]]$model,moment_match=TRUE)

m[[6]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "bs", n_knots = 30,par_list="all")
loos[[6]] = loo::loo(m[[6]]$model,moment_match=TRUE)

# fit models with varying numbers of knots in a gp
m[[7]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 6, gp_theta_prior = c(5,5),par_list="all")
loos[[7]] = loo::loo(m[[7]]$model,moment_match=TRUE)

m[[8]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 12, gp_theta_prior = c(5,5),par_list="all")
loos[[8]] = loo::loo(m[[8]]$model,moment_match=TRUE)

m[[9]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 18, gp_theta_prior = c(5,5),par_list="all")
loos[[9]] = loo::loo(m[[9]]$model,moment_match=TRUE)

m[[10]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 24, gp_theta_prior = c(5,5),par_list="all")
loos[[10]] = loo::loo(m[[10]]$model,moment_match=TRUE)

m[[11]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = 30, gp_theta_prior = c(5,5),par_list="all")
loos[[11]] = loo::loo(m[[11]]$model,moment_match=TRUE)

m[[12]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "gp", n_knots = max(x$time), gp_theta_prior = c(5,5),par_list="all")
loos[[12]] = loo::loo(m[[12]]$model,moment_match=TRUE)

# fit models with varying numbers of knots in a gp
m[[13]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "ps", n_knots = 6, gp_theta_prior = c(5,5),par_list="all")
loos[[13]] = loo::loo(m[[13]]$model,moment_match=TRUE)

m[[14]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "ps", n_knots = 12, gp_theta_prior = c(5,5),par_list="all")
loos[[14]] = loo::loo(m[[14]]$model,moment_match=TRUE)

m[[15]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "ps", n_knots = 18, gp_theta_prior = c(5,5),par_list="all")
loos[[15]] = loo::loo(m[[15]]$model,moment_match=TRUE)

m[[16]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "ps", n_knots = 24, gp_theta_prior = c(5,5),par_list="all")
loos[[16]] = loo::loo(m[[16]]$model,moment_match=TRUE)

m[[17]] <- fit_dfa(y = x, data_shape = "long", iter = mcmc_iter, thin = mcmc_thin,chains = mcmc_chains, num_trends = 1, trend_model = "ps", n_knots = 30, gp_theta_prior = c(5,5),par_list="all")
loos[[17]] = loo::loo(m[[17]]$model,moment_match=TRUE)

saveRDS(m, "output/calcofi_models.rds")
saveRDS(loos, "output/calcofi_loos.rds")

#m <- readRDS("output/calcofi_models.rds")
