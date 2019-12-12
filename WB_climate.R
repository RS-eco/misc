####################
#install.packages('rWBclimate')
library('rWBclimate')

#Tutorial
#https://ropensci.org/tutorials/rwbclimate_tutorial.html

get_model_temp("USA","mavg",2080,2100)  ## Get model temperature data
get_model_precip()

#Historical climate data of Africa
country.list <- Africa$ISO3
country.precip <- get_historical_precip(country.list,"year")
country.temp <- get_historical_temp(country.list,"year")

ggplot(country.precip,aes(x=year,y=data,group=locator)) +
  geom_point() +
  geom_path() +
  ylab("Annual precipitation of Africa") +
  theme_bw() +
  xlab("Year") +
  stat_smooth(se=F,colour="black") +
  facet_wrap(~locator,scales="free")

ggplot(country.temp,aes(x=year,y=data,group=locator)) +
  geom_point() +
  geom_path() +
  ylab("Average annual temperature of Africa") +
  theme_bw() +
  xlab("Year") +
  stat_smooth(se=F,colour="black") +
  facet_wrap(~locator,scales="free")

# Climate anomaly of Africa

# Set the kmlpath option
options(kmlpath = "~/Protected areas/Data/KML")

# Here we use a list basins for Africa
af_basin <- create_map_df(Africa_basin)

# Map of Africa
ggplot(af_basin, aes(x=long, y=lat,group=group)) +
  geom_polygon() +
  theme_bw()

# Download ensemble temperature data
af_basin_dat <- get_ensemble_temp(Africa_basin,"annualanom",2080,2100)

#  Subset data to just one scenario, and one percentile
af_basin_dat <- subset(af_basin_dat,af_basin_dat$scenario == "a2")
af_basin_dat <- subset(af_basin_dat,af_basin_dat$percentile == 50)

# Plot map of temperature anomalies of Africa
af_map <- climate_map(af_basin,af_basin_dat,return_map = T)
af_map +
  scale_fill_continuous("Temperature \n anomaly",low="yellow",high = "red") +
  theme_bw()

####################
