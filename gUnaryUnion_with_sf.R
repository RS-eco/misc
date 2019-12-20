install.packages("rworldmap")

data(countriesLow, package="rworldmap")

library(rgeos)
continent <- gUnaryUnion(countriesLow, id=countriesLow$continent)
plot(continent)

gUnaryUnion(countriesLow, as.character(countriesLow$continent))

summarise(group_by(countriesLow_sf, continent))

library(sf)
countriesLow_sf <- st_as_sf(countriesLow)

plot(countriesLow)

library(dplyr)
continent <- countriesLow %>% select(continent) %>% st_union(by_feature=T)

plot(continent)

#' This isn't an SO question, it's to do with sf being tidier than is good for it, 
#' and way slower than sp/rgeos - correct me if I'm wrong. 
#' The wiki says use sf::aggregate.sf(), so I'll try below - I use aggregate() in teaching. 
#' Tim's example omits dplyr (needed for group_by anyway). 
#' This is anyway a lot slower than rgeos::gUnaryUnion():

system.time(gUU <- gUnaryUnion(countriesLow, as.character(countriesLow$continent)))

#   user  system elapsed 
#  1.882   0.000   1.886 
system.time(tidy <- summarise(group_by(countriesLow_sf, continent)))
#   user  system elapsed 
#  8.942   0.000   8.956 

#' It also returns a tibble, which should always be avoided, as it pretends to be a data.frame but can created unadvertised havoc in modelling. 
#' If someone wants a tibble, they should coerce explicitly - the input and output objects should be the same class (as in aggregate()).

#' This seems to work but returns only a data frame, and possibly discards four of the six output features although fastest:

system.time(continent1 <- aggregate(st_geometry(countriesLow_sf), list(countriesLow_sf$continent), head, 1))
#    user  system elapsed 
#  0.014   0.000   0.014 

#' This is also time-consuming, but fortunately does not coerce to tibble, and would be my preference despite taking a lot longer:

system.time(continent2 <- aggregate(countriesLow_sf[,"geometry"], list(countriesLow_sf$continent), head, 1))
#   user  system elapsed 
#  8.896   0.000   8.909 

#' This feature-only unary union is also slow and does not identify the features (gUnaryUnion inserts IDs in ID slots, accessed by row.names()):

system.time(continent3 <-st_as_sfc(tapply(st_geometry(countriesLow_sf), list(countriesLow_sf$continent), st_union)))
#   user  system elapsed 
#  8.874   0.000   8.887

#' OK, thanks:

system.time(continent1 <- st_union(select(countriesLow_sf, continent), by_feature=TRUE))
#   user  system elapsed 
#  0.661   0.000   0.667

#' This fixes propagation of precision through sf:::aggregate.sf

a <- aggregate(st_set_precision(countriesLow, 1e6), list(countriesLow$continent), function(f) f[1])

#' Timing is now also similar:

system.time(a <- aggregate(st_set_precision(countriesLow, 1e6), list(countriesLow$continent), function(f) f[1]))
#   user  system elapsed 
#  1.959   0.000   1.959 

data(countriesLow, package="rworldmap")
system.time(continent <- gUnaryUnion(countriesLow, id=countriesLow$continent))
#   user  system elapsed 
#  1.986   0.000   1.986 

#' The dplyr pipeline with st_set_precision(1e8) inserted did not have this problem, and has identical timing.

countriesLow = st_as_sf(countriesLow)
system.time(continent <- countriesLow %>% st_set_precision(1e8) %>% select(continent) %>% 
              group_by(continent) %>% summarize())

#   user  system elapsed 
#   1.94    0.00    1.94 

#' Interestingly, it gives a seventh continent:

continent

# Simple feature collection with 7 features and 1 field
# geometry type:  MULTIPOLYGON
# dimension:      XY
# bbox:           xmin: -180 ymin: -89.9989 xmax: 180 ymax: 83.5996
# epsg (SRID):    4326
# proj4string:    +proj=longlat +datum=WGS84 +no_defs
# precision:      1e+08 

# A tibble: 7 x 2
#  continent                                                            geometry
#  <fct>                                                      <MULTIPOLYGON [°]>
#1 Africa        (((37.85695 -46.94427, 37.81395 -46.96287, 37.61179 -46.94644,…
#2 Antarctica    (((-147.5883 -76.64984, -147.5788 -76.66276, -147.7297 -76.653…
#3 Australia     (((158.8788 -54.70974, 158.8452 -54.74927, 158.8359 -54.704, 1…
#4 Eurasia       (((-26.26413 -58.43514, -26.25984 -58.4923, -26.41533 -58.4398…
#5 North America (((-155.5813 19.01202, -155.6256 18.96391, -155.6808 18.96768,…
#6 South America (((-67.5752 -55.88966, -67.61142 -55.89173, -67.69953 -55.8731…
#7 NA            (((73.70747 -53.13707, 73.58799 -53.18461, 73.46511 -53.1842, …

#' It seems it is a factor vs character issue:

library(sf)
library(dplyr)

data(countriesLow, package = "rworldmap")
countriesLow = st_as_sf(countriesLow) %>% 
select(continent)
#> Loading required package: sp

system.time({
a = aggregate(st_set_precision(countriesLow, 1e8),
list(countriesLow$continent),
function(f) f[1])
})
#>    user  system elapsed 
#>   1.290   0.004   1.305

system.time({
b = countriesLow %>%
st_set_precision(1e8) %>%
group_by(continent) %>% 
summarize()
})
#>    user  system elapsed 
#>   1.172   0.001   1.183

system.time({
c = aggregate(st_set_precision(countriesLow, 1e8),
list(as.character(countriesLow$continent)),
function(f) f[1])
})
#>    user  system elapsed 
#>   1.216   0.002   1.227

dim(a)
#> [1] 6 3
dim(b)
#> [1] 7 2
dim(c)
#> [1] 6 3

#' ... and having nothing to do with sf:

aggregate(a~f, data.frame(a = 1:3, f = factor(c("a", "b", NA))), mean)
#  f a
#1 a 1
#2 b 2
library(dplyr)
data.frame(a = 1:3, f = factor(c("a", "b", NA)))  %>% group_by(f) %>% summarise(mean(a))
# A tibble: 3 x 2
#  f     `mean(a)`
#  <fct>     <dbl>
#1 a             1
#2 b             2
#3 NA            3