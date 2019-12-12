# Interrupted Goode homolosine

library(tidyverse)
library(cowplot)   # for theme_minimal_grid()
library(sf)        # for manipulation of simple features objects
library(rworldmap) # for getMap()

# Obtain basic world map
world_sf <- st_as_sf(getMap(resolution = "low"))

# Goode homolosine projection
crs_goode <- "+proj=igh"

# projection outline in long-lat coordinates
lats <- c(
  90:-90, # right side down
  -90:0, 0:-90, # third cut bottom
  -90:0, 0:-90, # second cut bottom
  -90:0, 0:-90, # first cut bottom
  -90:90, # left side up
  90:0, 0:90, # cut top
  90 # close
)
longs <- c(
  rep(180, 181), # right side down
  rep(c(80.01, 79.99), each = 91), # third cut bottom
  rep(c(-19.99, -20.01), each = 91), # second cut bottom
  rep(c(-99.99, -100.01), each = 91), # first cut bottom
  rep(-180, 181), # left side up
  rep(c(-40.01, -39.99), each = 91), # cut top
  180 # close
)

goode_outline <- 
  list(cbind(longs, lats)) %>%
  st_polygon() %>%
  st_sfc(
    crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  )

# now we need to work in transformed coordinates, not in long-lat coordinates
goode_outline <- st_transform(goode_outline, crs = crs_goode)

# get the bounding box in transformed coordinates and expand by 10%
xlim <- st_bbox(goode_outline)[c("xmin", "xmax")]*1.1
ylim <- st_bbox(goode_outline)[c("ymin", "ymax")]*1.1

# turn into enclosing rectangle
goode_encl_rect <- 
  list(
    cbind(
      c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]), 
      c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
    )
  ) %>%
  st_polygon() %>%
  st_sfc(crs = crs_goode)

# calculate the area outside the earth outline as the difference
# between the enclosing rectangle and the earth outline
goode_without <- st_difference(goode_encl_rect, goode_outline)

# Plot map
ggplot(world_sf) + 
  geom_sf(color = "black", size = 0.5/.pt) +
  geom_sf(data = goode_without, fill = "white", color = "NA") +
  geom_sf(data = goode_outline, fill = NA, color = "gray30", size = 0.5/.pt) +
  coord_sf(crs = crs_goode, xlim = 0.95*xlim, ylim = 0.95*ylim, expand = FALSE) +
  theme_minimal_grid()

# Plot map with raster
library(raster)
r1 <- raster(nrows=360, ncols=720)
values(r1) <- rnorm(ncell(r1), 10, 40)
r1 <- projectRaster(r1, crs=crs_goode)
df <- as.data.frame(rasterToPoints(r1))

ggplot() + geom_raster(data=df, aes(x=x,y=y, fill=layer)) +
  geom_sf(data=world_sf, fill=NA, color = "black", size = 1/.pt) +
  geom_sf(data = goode_without, fill = "white", color = "NA") +
  geom_sf(data = goode_outline, fill = NA, color = "black", size = 1/.pt) +
  coord_sf(crs = crs_goode, xlim = 0.95*xlim, ylim = 0.95*ylim, expand = FALSE, ndiscr=200) +
  labs(x="", y="")
