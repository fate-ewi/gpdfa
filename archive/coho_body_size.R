d = readRDS("data/coho_joined_releases_recoveries_locations.rds")

d$age = as.numeric(substr(d$recovery_date,1,4)) - d$brood_year
d = dplyr::filter(d, fishery==10,
  recovery_rmis_region %in% c("CAGN","GRAY","NOOR","NWC","SOOR","WAGN","WILP"),
  age==3,
  !is.na(length))

d$recovery_year = as.numeric(substr(d$recovery_date,1,4))

stocks = dplyr::group_by(d, brood_year, release_location_rmis_basin) %>%
  dplyr::summarize(n = n()) %>%
  dplyr::group_by(release_location_rmis_basin) %>%
  dplyr::summarize(mean_recov = mean(n), n_thresh = length(which(n > 10)), u = length(unique(brood_year))) %>%
  dplyr::filter(n_thresh > 20)

sub = dplyr::filter(d, release_location_rmis_basin%in% stocks$release_location_rmis_basin) %>%
  dplyr::group_by(brood_year, release_location_rmis_basin) %>%
  dplyr::summarize(n = n(),mean_length = mean(length))
df = expand.grid(brood_year = unique(sub$brood_year),
  release_location_rmis_basin = unique(sub$release_location_rmis_basin))
df = left_join(df, sub)

df = dplyr::rename(df, ts = release_location_rmis_basin, time = brood_year,
  obs = mean_length) %>%
df = dplyr::filter(df, !is.na(obs))

ggplot(df, aes(time,obs,group=ts,col=ts)) + geom_line()

m1 = fit_dfa(y = df, data_shape="long", iter = 3000, chains = 1, num_trends = 1)
m2 = fit_dfa(y = df, data_shape="long", iter = 3000, chains = 1, num_trends = 2)

m3 = fit_dfa(y = df, data_shape="long", iter = 3000, chains = 1, trend_model="gp", n_knots = 30, num_trends = 1)
m4 = fit_dfa(y = df, data_shape="long", iter = 3000, chains = 1, trend_model="gp", n_knots = 30, num_trends = 2)

m5 = fit_dfa(y = df, data_shape="long", iter = 3000, chains = 1, trend_model="spline", n_knots = 30, num_trends = 1)
m6 = fit_dfa(y = df, data_shape="long", iter = 3000, chains = 1, trend_model="spline", n_knots = 30, num_trends = 2)

loos = NA
for(i in seq(31,40)) {
m = fit_dfa(y = df, data_shape="long", iter = 3000, chains = 1, trend_model="spline", n_knots = i, num_trends = 2)
loos[i] = loo(m)$estimates["looic","Estimate"]
}



d = readRDS("data/coho_joined_releases_recoveries_locations.rds")
# throw out lengths above or below 4 sds
d$age = as.numeric(substr(d$recovery_date,1,4)) - d$brood_year
d = dplyr::filter(d, fishery==50,
  release_location_state =="WA",
  age==3,
  sex=="F",
  !is.na(length),
  length > 304.718,
  length < 867.77)

d$recovery_year = as.numeric(substr(d$recovery_date,1,4))

sub = d %>%
  dplyr::group_by(brood_year, release_location_rmis_basin) %>%
  dplyr::summarize(n = n(),mean_length = mean(length,na.rm=T),
    region = release_location_rmis_region[1])
#df = expand.grid(brood_year = unique(sub$brood_year),
#  release_location_rmis_basin = unique(sub$release_location_rmis_basin))
df = dplyr::filter(sub, n >= 100)
df = dplyr::rename(df, ts = release_location_rmis_basin, time = brood_year, obs = mean_length) %>%
  dplyr::filter(!is.na(obs)) %>%
  dplyr::filter(ts %in% c("DUWA","EKPS","LUDA","NISQ","NOOK","PUYA","SKDO","SNOH"))

df = dplyr::filter(df, time >= 1980)

ggplot(df, aes(time, obs, group=ts,col=region))+geom_line() + xlab("Brood year") +
  ylab("Log lengths")

df$obs = log(df$obs)

m1 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, num_trends = 1)
m2 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, num_trends = 2)

m3 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="gp", n_knots = 9, num_trends = 1)
m4 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="gp", n_knots = 9, num_trends = 2)

m5 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="spline", n_knots = 30, num_trends = 1)
m6 = fit_dfa(y = df, data_shape = "long", iter = 3000, chains = 1, trend_model="spline", n_knots = 43, num_trends = 2)

