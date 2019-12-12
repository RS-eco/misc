# Load packages
library(plotrix)
library(tidyverse)
library(stringr)
library(threejs)
library(maps)
library(readr)
library(dplyr)

# Get the image of the globe - NASA
earth <- "http://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73909/world.topo.bathy.200412.3x5400x2700.jpg"

#Display the data on the globe
globejs(
 img = earth, lat = 35.68, long = 139.69,
 arcsOpacity = 0.20, arcsHeight = 0.5, arcsLwd = 1,
 arcsColor = "orange", atmosphere = TRUE, height = 700,
 width = 700, bg = "white", value = 4
)

data(flights)
# Approximate locations as factors
dest   <- factor(sprintf("%.2f:%.2f",flights[,3], flights[,4]))
# A table of destination frequencies
freq <- sort(table(dest), decreasing=TRUE)
# The most frequent destinations in these data, possibly hub airports?
frequent_destinations <- names(freq)[1:10]
# Subset the flight data by destination frequency
idx <- dest %in% frequent_destinations
frequent_flights <- flights[idx, ]
# Lat/long and counts of frequent flights
ll <- unique(frequent_flights[,3:4])

# Plot frequent destinations as bars, and the flights to and from
# them as arcs. Adjust arc width and color by frequency.
globejs(lat=ll[,1], long=ll[,2], arcs=frequent_flights, bodycolor="#aaaaff",
        arcsHeight=0.3, arcsLwd=2, arcsColor="#ffff00", arcsOpacity=0.15,
        atmosphere=TRUE, color="#00aaff", pointsize=0.5)

library(threejs)
library(maps)
data(world.cities, package="maps")
cities <- world.cities[order(world.cities$pop, decreasing=TRUE)[1:1000],]
value  <- 100 * cities$pop / max(cities$pop)
col <- colorRampPalette(c("cyan", "lightgreen"))(10)[floor(10 * value/100) + 1]
globejs(lat=cities$lat, long=cities$long, value=value, 
        color=col, atmosphere=TRUE)

# Plot the data on the moon:
moon <- system.file("images/moon.jpg", package="threejs")
globejs(img=moon, bodycolor="#555555", lightcolor="#aaaaaa",
        lat=cities$lat, long=cities$long,
        value=value, color=col)
