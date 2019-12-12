# Basic remote sensing tasks in R

# Data aquired from USGS, Landsat 8 (bands are different from Landsat 7 -> +1)
# Times points are February 2015 and July 2015

# aim: (1) affetion of tradewinds on greening of the northwestern side?
# (2) to prove that La Palma is in fact la Isla bonita with wonderful vegetation throughout 
#     the year (comparison of winter - summer images)
# (3) what are the ecosystem types that contribute most (area) to the greeness


# load required packages for the whole task
library("rgdal")
library("randomForest")
library("RStoolbox")
library("caret")
library("raster")
library("landsat")
library("landsat8")

# Extent: coordinates derived from QGIS analysis
# SO 237857.732,3148660.798
# NW 200326.588,3198826.188

# images La Palma, February 2015
### trying to get 2 different pictures Jan, Jul -> is in fact the green isle
# with this data set, i have to do a layer stack first, to import all band combos
Feb_1<-raster("Februar/LC82080402014010LGN00_B1.tif")
crop<-c(200326.588,237857.732,3148660.798,3198826.188) # desired extent of study area
Feb1crop<-crop(Feb_1,crop) # only region of interest cropped

Feb_2<-raster("Februar/LC82080402014010LGN00_B2.tif")
Feb2crop<-crop(Feb_2,crop)

Feb_3<-raster("Februar/LC82080402014010LGN00_B3.tif")
Feb3crop<-crop(Feb_3,crop)

Feb_4<-raster("Februar/LC82080402014010LGN00_B4.tif")
Feb4crop<-crop(Feb_4,crop)

Feb_5<-raster("Februar/LC82080402014010LGN00_B5.tif")
Feb5crop<-crop(Feb_5,crop)

Feb_6<-raster("Februar/LC82080402014010LGN00_B6.tif")
Feb6crop<-crop(Feb_6,crop)

Feb_7<-raster("Februar/LC82080402014010LGN00_B7.tif")
Feb7crop<-crop(Feb_7,crop)

# mask for clouds
# turn raster into grid for landsat package
g <- as(Feb2crop, 'SpatialGridDataFrame')
f <- as(Feb7crop, 'SpatialGridDataFrame')

# mask the clouds
cloudsFeb<-clouds(g,f)
image(cloudsFeb) # sample data in landsat package has the same result
# tried with RStoolbox, but since stacking is not possible for different resolutions, it is not
# possible to mask the clouds with cloudMask; Landsat TIR bands have 100m resolution
# so clouds will be in the classification afterwards marked as clouds

## Import Landsat example subset
data(lsat) 
## We have two tiny clouds in the east
ggRGB(lsat, stretch = "lin")

## Calculate cloud index
cldmsk    <- cloudMask(lsat, blue = 1, tir = 6)
ggR(cldmsk, 2, geom_raster = TRUE) 

## Define threshold (re-use the previously calculated index)
## Everything above the threshold is masked
## In addition we apply a region-growing around the core cloud pixels
cldmsk_final <- cloudMask(cldmsk, threshold = 0.1, buffer = 5) 
ggRGB(lsat, stretch = "lin") +
  ggR(cldmsk_final[[1]], ggLayer = TRUE, forceCat = TRUE, geom_raster = TRUE) +
  scale_fill_manual(values = "red", na.value = NA)

## Non-interactively. Pre-defined shadow displacement estimate (shiftEstimate)
shadow <- cloudShadowMask(lsat, cldmsk_final, shiftEstimate = c(-16,-6))

## Plot
csmask <- raster::merge(cldmsk_final[[1]], shadow)
ggRGB(lsat, stretch = "lin") +
  ggR(csmask, ggLayer = TRUE, forceCat = TRUE, geom_raster = TRUE) +
  scale_fill_manual(values = c("blue", "yellow"), 
                    labels = c("shadow", "cloud"), na.value = NA)

# stack them together
Feb_allbands<-stack(Feb1crop,Feb2crop,Feb3crop,Feb4crop,Feb5crop,Feb6crop,Feb7crop)
plot(Feb_allbands)
writeRaster(Feb_allbands,"Februar/Feb_allbands.tif","GTiff")
RGB<-stack(Feb4crop,Feb3crop,Feb2crop)
writeRaster(RGB,"Februar/RGB.tif","GTiff")
All_Feb<-brick("Februar/Feb_allbands.tif")
writeRaster(Feb_allbands,"Februar/Feb_allbands.grd")
# Feb<-brick("Februar/Feb_allbands.tif")
# plot(Feb)


# NDVI of La Palma
# modify to Landsat 8
ndvi<-(Feb5crop-Feb4crop)/(Feb5crop+Feb4crop) # calculating NDVI
plot(ndvi)
writeRaster(ndvi,"Februar/ndvi_Feb.tif","GTiff")
All_Feb_ndvi<-stack(All_Feb,ndvi)

#########
# 
# check for haziness of the image
metaData <- readMeta("Februar/LC82080402014010LGN00_MTL.txt")
lsat <- stackMeta(mtlFile)

## Convert DN to top of atmosphere reflectance and brightness temperature
lsat_ref <- radCor(All_Feb_ndvi, metaData = meta, method = "apref")
## Correct DN to at-surface-reflecatance with DOS (Chavez decay model)
lsat_sref <- radCor(lsat, metaData = metaData, method = "dos")
## Correct DN to at-surface-reflecatance with simple DOS
## Automatic haze estimation
hazeDN <- estimateHaze(lsat, hazeBands = 1:4, darkProp = 0.01, plot = TRUE)
lsat_sref <- radCor(lsat, metaData = metaData, method = "sdos",
                    hazeValues = hazeDN, hazeBands = 1:4)


# Use QGIS in order to identify the different land classes by using different band combinations
# watch for the same coordinate system in shapefile

# natural colour: 432 (Band 4, Band 3, Band 2)
# false colour (urban): 764
# Agriculture: 652
# atmospheric penetration: 765
# healthy vegetation: 562
# Land/water: 564
# natural with atmospheric removal: 753
# shortwave infrared: 754
# vegetation analysis: 654
# Source: ArcGIS user manual

# do the landclasses of La Palma
Landclasses<-readOGR("Shp","LandclassesFeb4")

# plot Landclasses
plot(All_Feb,1)
plot(Landclasses, add=TRUE)

# trainings and validation data
saveRDS(Landclasses,"Landclasses.rds")

train<-readRDS("Landclasses.rds")

# fit classification
SC<-superClass(All_Feb_ndvi, nSamples=1000, trainData = train, responseCol = "id",
               model = "rf", trainPartition = 0.7) # classification with random forest
SC  

# plot result
plot(SC$map)


# write raster
writeRaster(SC$map,"Februar/classification_Feb.tif","GTiff")



###########################
##########################
# July 2015

Jul_1<-raster("July/LC82080402015189LGN00_B1.tif")
crop<-c(200326.588,237857.732,3148660.798,3198826.188) # desired extent of study area
Jul1crop<-crop(Jul_1,crop) # only region of interest cropped

Jul_2<-raster("July/LC82080402015189LGN00_B2.tif")
Jul2crop<-crop(Jul_2,crop)

Jul_3<-raster("July/LC82080402015189LGN00_B3.tif")
Jul3crop<-crop(Jul_3,crop)

Jul_4<-raster("July/LC82080402015189LGN00_B4.tif")
Jul4crop<-crop(Jul_4,crop)

Jul_5<-raster("July/LC82080402015189LGN00_B5.tif")
Jul5crop<-crop(Jul_5,crop)

Jul_6<-raster("July/LC82080402015189LGN00_B6.tif")
Jul6crop<-crop(Jul_6,crop)

Jul_7<-raster("July/LC82080402015189LGN00_B7.tif")
Jul7crop<-crop(Jul_7,crop)

# ndvi July
ndvi_Jul<-(Jul5crop-Jul4crop)/(Jul5crop+Jul4crop) # calculating NDVI
plot(ndvi_Jul)
writeRaster(ndvi_Jul,"July/ndvi_Jul.tif","GTiff",overwrite=TRUE)



# stack them together
Jul_allbands<-stack(Jul1crop,Jul2crop,Jul3crop,Jul4crop,Jul5crop,Jul6crop,Jul7crop,ndvi_Jul)
plot(Jul_allbands)
writeRaster(Jul_allbands,"July/Jul_allbands.tif","GTiff",overwrite=TRUE)
RGB_Jul<-stack(Jul4crop,Jul3crop,Jul2crop)
writeRaster(RGB_Jul,"July/RGB_Jul.tif","GTiff",overwrite=TRUE)
All_Jul<-brick("July/Jul_allbands.tif")


# do the landclasses of La Palma
LandclJul<-readOGR("Shp","LandclassesJul")

# plot Landclasses
plot(All_Jul,1)
plot(LandclJul, add=TRUE)


# trainings and validation data
saveRDS(LandclJul,"LandclassesJul.rds")

train<-readRDS("LandclassesJul.rds")


# fit classification
SCJul<-superClass(All_Jul, nSamples=1000, trainData = train, responseCol = "id",
               model = "rf", trainPartition = 0.7) # classification with random forest
SCJul

# plot result
plot(SCJul$map)


# write raster
writeRaster(SCJul$map,"July/classification_Jul.tif","GTiff",overwrite=TRUE)

######################
# calculate differences in ndvi
ndvi_diff<-ndvi_Jul-ndvi
plot(ndvi_diff)
writeRaster(ndvi_diff,"ndvi_difference.tif","GTiff")
