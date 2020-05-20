#' ---
#' title: "Pre-process European climate and Odonata data"
#' author: "Matthias F. Biber"
#' ---

#########################

## Set-up

# Set working directory
setwd("/home/matt/Documents/Sync+Share/@BayKliF/Daten/")

# Load packages
rm(list=ls()); invisible(gc())
library(dplyr); library(rnaturalearthhires); library(sf)

# Load high resolution country shapefile
data(countries10)
countries10 <- sf::st_as_sf(countries10)

#########################

# Read climate data
bioclim <- read.csv("climate_data/cordex_bioclim_eur.csv.xz", sep=",")

# Split data by scenario and timeframe to be able to use rasterFromXYZ
bioclim <- bioclim %>% tidyr::unite(col="gcm_rcm_rcp", c("gcm", "ensemble", "rcm", "rs", "rcp"), sep="_") %>% 
  group_split(gcm_rcm_rcp, time_frame, keep=F)

#' The rasterFromXYZ function requires that your data.frame is structured in a very specific way,
#' i.e. first column is x, second column is y and the remaining columns can only be numeric variables.

# Turn data into raster
r_bioclim <- lapply(bioclim, raster::rasterFromXYZ)
raster::plot(r_bioclim[[1]][[1]])
plot(sf::st_geometry(countries10), add=T)

r_bioclim[[1]]

#' Data has no CRS

# Define CRS
r_bioclim <- lapply(r_bioclim, function(x){
  raster::crs(x) <- "+proj=longlat +datum=WGS84 +no_defs"
  return(x)
})

r_bioclim[[1]]

#' Now data has correct CRS

#########################

# Read Odonata data
odonata <- readxl::read_xlsx("Kalkman et al._Eur_Dragonfly_AtlasData.xlsx")
odonata$Latitude <- as.numeric(gsub(",", ".", odonata$Latitude))
odonata$Longitude <- as.numeric(gsub(",", ".", odonata$Longitude))

# Turn Odonata data into Spatial Points
sf_odonata <- sf::st_as_sf(odonata, coords = c("Longitude", "Latitude"), crs = 4326)
plot(sf::st_geometry(sf_odonata %>% filter(`Species name`=="Aeshna affinis")))
plot(sf::st_geometry(countries10), add=T)

# Project data to UTM
source("/home/matt/Documents/LatLonToUTMEPSGCode.R")

sf_odonata_utm <- sf_odonata %>%
  do(cbind(., st_coordinates(.))) %>%
  mutate(EPSG = LatLonToUTMEPSGCode(lat = Y, lon = X)) %>%
  group_by(EPSG) %>%
  do(cbind(as.data.frame(.),
           st_coordinates(st_transform(., crs = head(.$EPSG, 1))))) %>%
  ungroup()
epsg_names <- unique(sf_odonata_utm$EPSG)
sf_odonata_utm <- sf_odonata_utm %>% group_split(EPSG)

sf_odonata_utm <- lapply(1:length(sf_odonata_utm), function(x){
  dat <- sf::st_as_sf(sf_odonata_utm[[x]], coords = c("X", "Y"), crs = epsg_names[x])
  return(dat)
})
# sf_odonata_utm[[10]]

# Now we plot the data in UTM for EPSG 32364 and species Aeshna affinis
plot(st_geometry(sf_odonata_utm[[10]] %>% filter(Species.name == "Aeshna affinis")))

# Transform and plot the country shapefile in the same EPSG
plot(st_transform(st_geometry(countries10), epsg_names[10]), add=T)

#' **Note:** Plotting the country shapefile with the same command as before, 
#' won't show you anything, as the coordinate reference system needs to match
#' 

# Turn odonata_utm data into raster files
r_odonata_utm <- lapply(1:length(sf_odonata_utm), function(x){
  dat <- as.data.frame(as(sf_odonata_utm[[x]], "Spatial")) %>% select(coords.x1, coords.x2, Species.name) %>%
    mutate(coords.x1 = round(coords.x1, 0), coords.x2 = round(coords.x2, 0)) %>%
    mutate(presence = 1) %>% group_split(Species.name)
  r_dat <- lapply(dat, function(y){
    raster::rasterFromXYZ(y %>% dplyr::select(-Species.name), res=c(50000, 50000), digits=0, 
                          crs=raster::crs(paste0("+init=epsg:", epsg_names[x])))
  })
  r_ext <- lapply(r_dat, function(k){
    ext <- raster::extent(k)
    ext <- data.frame(xmn=ext@xmin, xmx=ext@xmax,
                      ymn=ext@ymin, ymx=ext@ymax)
    return(ext)
  })
  r_ext <- bind_rows(r_ext)
  r <- raster::raster(xmn=min(r_ext$xmn), xmx=max(r_ext$xmx), 
                      ymn=min(r_ext$ymn), ymx=max(r_ext$ymx), res=c(50000, 50000),
                      crs=raster::crs(paste0("+init=epsg:", epsg_names[x])))
  r_dat <- lapply(r_dat, function(z) raster::extend(z, r))
  r_dat <- lapply(r_dat, function(z) raster::resample(z, r))
  r_dat <- raster::stack(r_dat)
  r_dat[r_dat >= 0.5] <- 1
  return(r_dat)
})
raster::plot(r_odonata_utm[[10]][[1]])

# Turn raster data into WGS84
r_odonata_wgs84 <- lapply(r_odonata_utm, function(x) raster::projectRaster(x, crs="+init=epsg:4236"))

# Merge different rasters by species, time frame, ...


#########################

# Re-project climate data to projection of Odonata data
r_bioclim <- projectRaster(r_bioclim, r_odonata)

# Aggregate climate data to resolution of Odonata data
r_bioclim <- aggregate(r_bioclim, r_odonata)

# Resample climate data to resolution of Odonata data
r_bioclim <- resample(r_bioclim, r_odonata)

# Turn bioclim data back to data.frame
