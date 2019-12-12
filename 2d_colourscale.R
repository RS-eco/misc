library(raster)

# See also this link for more information
# https://stackoverflow.com/questions/11070101/2d-color-gradient-plot-in-r

# get some climate data
clim <- getData('worldclim', var='bio', res=10, path="/media/matt/Data/Documents/Wissenschaft/Data")

# select precip of the wettest and driest months
clim <- subset(clim, c("bio13", "bio14"))

# crop to colorado
clim <- crop(clim, extent(-109, -102, 37, 41))

# convert to long-format matrix
v <- values(clim)

# map color to the climate variables
source("colormap.R")
colors <- colors2d(v, c("green", "yellow", "black", "blue"))

# plot the data in climate space
par(mfrow=c(1,2))
plot(v, col=colors, pch=16, 
     xlab="precip of wettest month (mm)", 
     ylab="precip of driest month (mm)")

# and in geographic space
col <- clim[[c(1,1,1)]]
col[] <- t(col2rgb(colors))
plotRGB(col)
