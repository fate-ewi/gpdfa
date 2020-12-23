library(bayesdfa)
library(ggplot2)
library(dplyr)

df = read.csv("data/upwelling.csv")
df$station = substr(df$POSITION,1,8)
df$year = as.numeric(substr(df$POSITION,10,13))

#df = dplyr::filter(df, station=="45N 125W")

ggplot(df, aes(year, MAR, group=station, col=station)) + geom_line()

df = df %>%
  dplyr::filter(year < 2020, year >= 2000) %>%
  dplyr::rename(obs = MAR, ts = station, time = year) %>%
  dplyr::filter(!is.na(obs)) %>%
  dplyr::select(ts, time, obs) %>%
  dplyr::filter(ts %in% c("21N 107W","24N 113W","60N 146W","60N 149W")==FALSE)

df = dplyr::filter(df, ts%in%c("Lower Umpqua River","Middle Umpqua River",
  "Siltcoos Lake","South Umpqua River","Tahkenitch Lake","Tenmile Lake")==FALSE)

df$ts = as.numeric(as.factor(df$ts))
m1 = fit_dfa(y = df, data_shape = "long", iter = 1000, chains = 1, num_trends = 1)
m2 = fit_dfa(y = df, data_shape = "long", iter = 1000, chains = 1, num_trends = 2)
m3 = fit_dfa(y = df, data_shape = "long", iter = 1000, chains = 1, num_trends = 3)

m3 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="gp", n_knots = 9, num_trends = 1)
m4 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="gp", n_knots = 9, num_trends = 2)

m5 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="spline", n_knots = 20, num_trends = 3)
m6 = fit_dfa(y = df, data_shape = "long", iter = 1000, chains = 1, trend_model="spline", n_knots = 20, num_trends = 2)

loos = NA
for(i in 10:20) {
  loos[i] = loo(fit_dfa(y = df, data_shape = "long", iter = 1000, chains = 1, trend_model="spline", n_knots = i, num_trends = 2))
}

fit = fit_dfa(y = df, data_shape = "long", iter = 1000, chains = 1, trend_model="gp", n_knots = 17, num_trends = 3,
  estimate_process_sigma = TRUE)

df = df %>% gather(month, upwelling, JAN:DEC) %>%
  dplyr::select(-POSITION, -station)
df$monthn = NA
df$monthn[which(df$month=="JAN")] = 1
df$monthn[which(df$month=="FEB")] = 2
df$monthn[which(df$month=="MAR")] = 3
df$monthn[which(df$month=="APR")] = 4
df$monthn[which(df$month=="MAY")] = 5
df$monthn[which(df$month=="JUN")] = 6
df$monthn[which(df$month=="JUL")] = 7
df$monthn[which(df$month=="AUG")] = 8
df$monthn[which(df$month=="SEP")] = 9
df$monthn[which(df$month=="OCT")] = 10
df$monthn[which(df$month=="NOV")] = 11
df$monthn[which(df$month=="DEC")] = 12

df = dplyr::filter(df, year >= 2000, year<2020)

ggplot(df, aes(year, Lake.Temp.Summer.Satellite, col=lake)) + geom_line()

library(bayesdfa)
# rename

df = dplyr::rename(df, obs = upwelling, ts = year, time = monthn)


m1 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, num_trends = 1)
m2 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, num_trends = 2)

m3 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="gp", n_knots = 10, num_trends = 1)
m4 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="gp", n_knots = 10, num_trends = 2)

m5 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="spline", n_knots = 10, num_trends = 1)
m6 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="spline", n_knots = 30, num_trends = 2)


df = read.csv(file.choose())
df$station = substr(df$POSITION,1,8)
df$year = as.numeric(substr(df$POSITION,10,13))


df = read.csv(file.choose())
df$station = substr(df$POSITION,1,8)
df$year = as.numeric(substr(df$POSITION,10,13))
df = dplyr::filter(df, station%in%c("42N 125W",
  "48N 125W","54N 134W"))
df = dplyr::filter(df, year >= 2000, year < 2020)

df = df %>% pivot_longer(names_to="month", values_to="upwelling", JAN:DEC) %>%
  dplyr::select(-POSITION)
df$monthn = NA
df$monthn[which(df$month=="JAN")] = 1
df$monthn[which(df$month=="FEB")] = 2
df$monthn[which(df$month=="MAR")] = 3
df$monthn[which(df$month=="APR")] = 4
df$monthn[which(df$month=="MAY")] = 5
df$monthn[which(df$month=="JUN")] = 6
df$monthn[which(df$month=="JUL")] = 7
df$monthn[which(df$month=="AUG")] = 8
df$monthn[which(df$month=="SEP")] = 9
df$monthn[which(df$month=="OCT")] = 10
df$monthn[which(df$month=="NOV")] = 11
df$monthn[which(df$month=="DEC")] = 12

df$station = as.numeric(as.factor(df$station))
df$year = df$year - min(df$year) + 1
df = dplyr::arrange(df, station, year, month)
df = dplyr::rename(df, obs = upwelling, ts = year, time = monthn)

