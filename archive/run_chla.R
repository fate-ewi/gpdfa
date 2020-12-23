require("ggplot2")
require("mapdata")
require("rerddap")
chlaInfo <- info('erdVHNchla3day')

df = data.frame("lat"=seq(33,48,3), lon = c(119,119,122,125,125,125))

for(i in 1:nrow(df)) {

viirsCHLA <- griddap(chlaInfo, latitude = c(df$lat[i], df$lat[i]), longitude = c(-df$lon[i], -df$lon[i]),
  time = c('2019-01-01','2019-12-31'), fields = 'chla')
if(i == 1) {
  chla = viirsCHLA$data
} else {
  chla = rbind(chla, viirsCHLA$data)
}
}

chla = dplyr::filter(chla, !is.na(chla))
chla$date = lubridate::as_date(chla$time)
chla$yday = lubridate::yday(chla$date)

ggplot(dplyr::filter(chla, lat==48.00375), aes(time, log(chla), group=lat, col=lat)) + geom_point()
