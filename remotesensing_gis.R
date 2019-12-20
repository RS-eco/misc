---
title: "Remote sensing and GIS for Ecologists"
---

#' #Spatial Data and Software

#' ## Basic introduction to R

# Load required packages
library(sp)

# Import the csv with x and y coordinates
gps_in <- read.table("/path/to/mytrack.csv")

# #### data integrity check

# check first entries
head(gps_in)

# overview statistics
summary(gps_in)

# plot data
plot(gps_in)

# data analysis

# ### data conversion for ...

#' ##Turning R into GIS

install.packages("raster") # execute only once to install a package
library("raster") # execute every time you start a new R session

#' ##Notation

#+ eval=FALSE
myAgg <- raster::aggregate(myRaster, 2)

#' # Where to obtain spatial data?

#' ## First spatial data for our study area

#' Let us download the vector outline of administrative boundaries of Brazil:

library(raster)
brazil <- getData("GADM", country = "BRA", level=1)
plot(brazil)

#' *Fig. 1.* Administrative boundaries (first level subdivision) of Brazil using the GADM data.

#' Such information can be retrieved for any selected country. 
#' Simply add the country acronym from the list by querying the row of your country:

x <- getData("ISO3")
x[x[, "NAME"] == "South Africa",]

#' You can also use the `getData()` command to also download raster data, such as elevation or climate variables:

prec <- getData("worldclim", var = "prec", res = 2.5)
names(prec) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
plot(prec, 1)

#' *Fig. 2.* The global precipitation pattern for January based on the Worldclim data set.

#' ### Landsat

#' Landsat data is available as scenes that are divided into paths and rows using the WRS2.
#' In R, the `wrspathrow`package offers a helpful way to obtain the corresponding path and row.

library(rgdal)
studyarea <- readOGR("vector_data/", "study_area_ll")

library(wrspathrow)
pathrowSA <- pathrow_sum(studyarea)
pathrowSA

#' In addition, we want to know where our study area is located within the Landsat scene:

pathrowSAsp <- pathrow_sum(studyarea, as_polys=TRUE)
plot(pathrowSAsp)
plot(studyarea, add=TRUE)

#' *Fig. 3.* Location of our study area within the Landscate scene (path 224, row 63).

#' ## Downloading data from Web portals

#' Global Land Cover Facility Earth Science Data Interface (ESDI) (http://www.landcover.org);
#' USGS Global Visualization Viewer (http://glovis.usgs.gov);
#' EarthExplorer (http://earthexplorer.usgs.gov);
#' Earth Observing System Data and Information System (http://reverb.echo.nasa.gov/reverb);
#' Portal of the Brazilian Space Agency INPE (http://www.dgi.inpe.br/CDSR/);
#' ESA Sentinel Data Hub (http://scihub.esa.int/)

#' ### Visualizing Landsat data availability using R 

#' This is achieved by searching the portal and downloading the scene list for the Landsat data you are looking for. These files, which hold the Landsat scene information, need to be unzipped:

zipfiles <- list.files(path = "raster_data/EarthExplorer/", pattern=".zip", full.names=TRUE)

#Unzip all files and return the file paths in a list
eeResults <- lapply(zipfiles, unzip, exdir="raster_data/EarthExplorer")

#' Now, we can import the EarthExplorer results using the `RStoolbox::readEE()`.

eeResults <- readEE(eeResults)

#' Then we plot the available Landscat scenes throught the year with the `ggplot2` package.

ggplot(subset(eeResults, Year < 2015 & Cloud.Cover < 20)) +
  geom_tile(aes(x = Doy, y = Year, alpha = Cloud.Cover, fill = Satellite),
            width=2, size=2) +
  scale_y_continous(breaks = c(1984:2014)) + 
  scale_alpha_continous(name = "Cloud Cover (%)", range = c(1, 0.5))

#' *Fig. 4.* Availability of Landsat data for the study area (Landsat TM, ETM+ and OLI) with cloud cover lower than 20%. TM = Thematic Mapper, ETM+ = Enhanced Thematic Mapper Plus, OLI = Operationl Land Imager.

#' ## Remote sensing products

#' Land cover, globally: GLC2000, GlobCover 2005, GlobCover 2009, MODIS Land Cover

#' GlobCover can be downloaded from the ESA website (http://due.esrin.esa.int/page_globcover.php)

#' Land cover, Europe: CORINE Land Cover (http://www.eea.europa.eu/data-and-maps#tab-datasets)

#' Land cover change data: **Global Forest Change** (http://earthengingepartners.appspot.com/science-2013-global-forest)

#' ### Spatial Non-remote sensing data

#' You will find an up-to-date list of available data portals at http://remote-sensing-biodiversity.org/resources.

#' The `letsR` package provides information on the geographic distribution of species and environmental variables based on the IUCN Red List of Threatened Species database. Several packages allow you to search for and retrieve species data including the Global Biodiversity Information facility (GBIF) (`rgbif`), the eBird database (`rebird`), FishBase (`rfishbase`) or from a variety of databases (`spooc`). There are many more packages, most of which are managed by the Transforming Science through Open Science project (https://www.ropensci.org).

#' Another valuble source of data is Open Street Map (OSM, http://www.openstreetmap.org).
#' The data can be access via R:


library(osmar)
studyarea_bb <- center_bbox(-49.75, -4.72, 30000, 30000)
OSMdata <- get_osm(studyarea_bb, source = osmsource_api())
plot(OSMdata)


#' More details on the OSM data set and how you can extract just roads or bike paths can be found on this web page: http://www.rpubs.com/RobinLovelace/12696. Moreover, a variety or further R packages are providing access to spatial data, such as `ocedata` for ocean data. World Bank climate data can further be easily retrieved with R using the `rWBclimate` package.

#' # Spatial data analysis for ecologists: First steps

#' ## Spatial data in R

#' ### Importing raster data

#' Let us start by importing the 2011 Landsat 5 TM imagery we downloaded. The Landsat data are delievered in a separate GeoTIFF file for each band. To load only a single band, we would need the following command:


p224r63_2011_B1 <- raster("raster_data/LT52240632011210/LT52240632011210CUB01_B1.TIF")


#' Printing the object name will return a handy summary of frequently needed information about the raster, such as clas, dimensions or projection and a summary of values for each layer:


p224r63_2011_B1


#' The `raster()` command can also be used to access single layers from multi-band files, by specifying the `band` argument. To import all the layers we could load each band separately and `stack()` the *RasterLayers* afterwards:


tmp1 <- raster("raster_data/LT52240632011210/LT52240632011210CUB01_B1.TIF")
tmp2 <- raster("raster_data/LT52240632011210/LT52240632011210CUB01_B2.TIF")
tmp3 <- raster("raster_data/LT52240632011210/LT52240632011210CUB01_B3.TIF")
# and so on ...

p224r63_2011 <- stack(tmp1, tmp2, tmp3)


#' But since this is time-consuming, we will simply list all files matching a given pattern and load them in one go with `stack()`, which can take a list of file paths as input:


allBands <- list.file("raster_data/LT52240632011210/", pattern = "TIF", full.names=TRUE)
p224r63_2011 <- stack(allBands)


#' Since the automatically created band names are rather long and redundant, we want to shorten them. Layer names can be queried, set and changed using `names()`. We simply remove the scene ID with search and replace:


names(p224r63_2011) <- gsub(pattern="LT52240632011210CUB01_", replace= "", x = names(p224r63_2011))p224r63_2011


#' Now that we have all bands in a single `Raster Stack`, we could save it to a single multi-band GeeoTIFF file on our hard disk:


writeRaster(pp224r63_2011, filename = "raster_data/LT52240632011210.tif")


#' Finally, to import all the layers from a multi-band file, we use a processed version of the current scene that was stored in the .grid file format of the raster package:


p224r63_2011 <- brick("raster_data/final/p224r63_2011.grd")


#' Only small rasters will be read fully into memory. To check the state of a raster object, let us see whether it is in memory or just linked to a file:


inMemory(p224r63_2011)


#' In addition to the classes defined by the raster package, you are likely to encounter *SpatialPixelsDataFrames* and *SpatialGridDataFrames* defined by the *sp* package.


p224r63_2011_sp <- as(p224r63_2011, "SpatialPixelsDataFrame")
class(p224r63_2011_sp)

p224r63_2011_ras <- as(p224r63_2011_sp, "RasterBrick")
# is equivalent to:
p224r63_2011_ras <- brick(p224r63_2011_sp)
class(p224r63_2011_ras)


#' Note that after the conversion from *sp* to raster classes the "*raster*" objects will remain loaded into memory. To free some memory, you will have to remove the object with `rm("p224r63_2011_ras")` or explicitly write them to disk:


p224r63_2011_ras <- writeRaster(p224r63_2011_ras, filename="raster_data/p224r63.tif")


#' ### Export raster data

#' As introduced earlier, saving raster data to disk should always be accomplished using `writeRaster()`.

#' Currently our *RasterBricks* links to the following file:


rasterfile <- filename(p224r63_2011_ras)
rasterfile


#' Now we will delete it externally using `file.remove()` and see what happens. The raster object still looks valid and would print its summary information just fine, but if we wanted to access some actual values it would generate an error message.


file.remove(rasterfile)
p224r63_2011_ras[1]


#' ### File formats

#' Many file formats from other vendors are very similar to the raster format. In fact, they often only differ in the format of the header file. We can make use of this by writing, for example, an additional ENVI header file (.hdr) to a raster grid file.
#' To write such a header, we use `raster::hdr()`:


hdr(p224r63_2011, format="ENVI")


#' Setting the global raster option `addheader` to "*ENVI*" will prompt `writeRaster()`to write an additional .hdr file whenever you save a .grd file, making it really easy to switch between R and QGIS:


rasterOptions(addheader="ENVI")


#' ### Importing vector data


vec <- readOGR(dsn = "path/to/directory", layer="myShape")
vec <- readOGR(dsn = ".", layer="filename")


#' Note that the filename is specified without a file extension.

#' Now, we import a shapefile containing a polygon that delineates our study area and a shapefile containing the outlines of two PAs.


studyArea <- readOGR("vector_data", "study_area_UTMnorth")
pa <- readOGR("vector_data", "PAs_UNEP_WCMC_p224r63")


#' In addition to the polygon data sets, you will also find a vector of the largest roads in the example data that we want to import:


roads <- readOGR("vector_data", "roads_p224r63_UTMsouth")


#' ### GPS Data

#' The syntax needed to import GPS data (.gpx files) is a little different than for shapefiles.


gpsPoints <- readOGR(dsn = "vector_data/field_measurements.gpx", layer="waypoints")


#' ### CSV Data

#' Often you will just have a text file to begin with, exported from your speadsheet software, with X and Y coordinates and some columns with data


csv_file <- read.csv("vector_data/csv_file_locationdata.csv")
csv_file


#' We will convert this object into a *SpatialPoints* object, because we know that it is point data, as opposed to line or polygon data:


csv.sp <- SpatialPoints(coords = csv_file[, c("X", "Y")])
csv.sp


#' As you can see, the *SpatialPoints* object does not have any values associated with it; it just holds the coordinates. Therefore, we build a *SpatialPointsDataFrame* that holds the data values in columns 3-5 in an attribute table:


csv.spdf <- SpatialPointsDataFrame(csv_file[, c("X", "Y")], data = csv_file[3:5])
csv.spdf


#' A helpful shortcut to transform a *data.frame* into a *SpatialPointsDataFrame* is to use the `coordinates()` function and specify the coordinate columns directly.


class(csv_file)

coordinates(csv_file) <- c("X", "Y")
class(csv_file)


#' Now that we have the spatial data structure all that is left is defining the projection.


projection(csv.spdf) <- projection(p224r63_2011)


#' ### Other vector classes

#' Sometimes the function you need resides in another package that does not use the sp *Spatial* classes introduced earlier. For example, to convert *SpatialPolygons* to a *PolySet*, two maptools commands can be used:


PBSpoly <- SpatialPolygons2PolySet(studyArea)

# and vice versa
spdf <- PolySet2SpatialPolygons(PBSpoly)


#' ### Extent objects


ex <- extent(p224r63_2011)
ex


#' Extent objects have a special characteristic; they can be easily resized.


ex*0.5


#' This can be very useful behaviour, for example, for a quick expansion of plotting dimensions without having to fiddle with coordinates or for quick subsetting of raster data to a smaller extent.

#' ### Exporting vector data


writeOGR(csv.spdf, "vector_data/processed/", "csv_spdf_as_shp", driver = "ESRI Shapefile")


#' Note that the column names in an ESRI Shapefile are limited by design to only 10 characters. If the column names are longer, they will be shortened and `writeOGR()` will print a warning.

#' The `writeOGR()` function that exports your vector data to *KML* should be like this:
#' Note that .kml files require unprojected coordinates, that is, latitude and longitude!


writeOGR(gpsPoints, dsn = "vector_data/processed/gpsPoints_GE.kml", layer="myLayername", driver = "KML")


#' Finally, you can export vector data to use them on your GPS device. This can be useful, for example, if you planned a random sampling scheme in R and need to find the selected sites in the field. Exporting into GPS Exchange Format (GPX) comes with a few restrictions. One of them is that you need a column labelled "name" and as with KML files you data must be in unprojected coordinates.


projection(gpsPoints)

"name" %in% names(gpsPoints)

writeOGR(gpsPoints, "vector_data/places_to_go.gpx", layer= "waypoints", driver = "GPX", dataset_options = "GPX_USE_EXTENSIONS=YES")


#' ## Plotting spatial data in R


# Plot all bands separately
plot(p224r63_2011)
# Plot only the fifth layer
plot(p224r63_2011, 5)
# Change the colour palette
plot(p224r63_2011, 5, col = grey.colors(100))
# Set layer transparency to 50%
plot(p224r63_2011, 5, alpha = 0.5)
# Alternative function (less flexible than plot and no legend)
image(p224r63_2011, 5)


#' If plotting large raster files is slow or your experience a memory shortage, you can use the *maxpixels* argument to reduce the number of pixels that are used for plotting.


plot(p224r63_2011, 5, maxpixels = 2e+05)


#' If you prefer the *lattice* plotting system, you can use `raster::spplot()`for *raster* and *sp* objects, for example, to plot the first four bands with a common scale:


spplot(p224r63_2011, 1:4)


#' When using the *ggplot2* package functions, the `RStoolbox::ggR()` command can be used to plot single layers with different options:


ggR(p224r63_2011, 5)

# with legend
ggR(p224r63_2011, 5, geom_raster=TRUE)

# with a custom legend
ggR(p224r63_2011, 5, geom_raster=TRUE) + scale_fill_gradientn(colours=rainbow(100))


#' For a true colour composite of the Landsat 5 TM scene, we assigned band three to *red*, band two to *green* and band one to *blue* and then defined the maximum value (scale):


plotRGB(p224r63_2011, r = 3, g = 2, b = 1, scale = 1)


#' If you plot is very dark or even just black, it might be due to the majority of values being very low. Therefore, different approaches so-called stretch methods are commonly used to increase contrast in the image.


plotRGB(p224r63_2011, r = 3, g = 2, b = 1, stretch = "lin")


#' or


ggRGB(p224r63_2011, r = 3, g = 2, b = 1, stretch = "lin")


#' We can also plot a false colour composite with a NIR band (4) to *red*, the red band (3) to *green* and the breen band (2) to *blue*. This band combination is used frequently when highlighting differences in vegetation:


plotRGB(p224r63_2011, r = 4, g = 3, b = 2, stretch = "lin")


#' Next, we want to add vector data as additional layers to our plots.


plot(p224r63_2011, 5)
points(csv.spdf, col = "blue") # equivalent to plot(..., add=TRUE)
plot(studyArea, add=TRUE, lwd= 10) # lwd=linewidth


#' The plot command also supports plotting *extent* objects:


plot(roads)
plot(extent(roads), add=TRUE)


#' When using the *ggplot2* system, you need to turn the *Spatial* classes into non-spatial data frames first. For *SpatialPointsDataFrames*, simply do the following:


df_pts <- as.data.frame(csv.spdf)


#' For the more complex *SpatialLines* amd *SpatialPolygons*, *ggplot2* provides a function called `fortify()` which will return a properly structured data frame ready for plotting


df_poly <- fortify(studyArea)
head(df_poly, 3)


#' To use rasters in conjunction with *ggplot2*, you have three options. You can can use `ggR()/ggRGB()` or convert the whole raster to a dataframe (not recommended for large rasters) or you use `fortify()` (currently provided by *RStoolbox*), which subsamples your raster just like the *raster* plot command does:


df_ras <- as.data.frame(p224r63_2011, xy=TRUE)
df_ras <- fortify(p224r63_2011, maxpixels= 1e+05)



ggRGB(p224r63_2011, 4, 3, 2, stretch="lin") + 
  geom_point(data = df_pts, aes(x = X, y = Y), size = 5, col = "yellow") + 
  geom_path(data = df_poly, aes(x = long, y = lat, group = group), size = 2, col = "blue") + 
  coord_equal()


#' Apart from the packages and functions we have discussed so far, the `rasterVis` package provides several useful features for data visualization, including time sieres visualization. Digital elevation models, for example, can be visualized in an interactive three-dimensional plot with `plot3D`.

#' ## Basic Spatial Data Handling in R

#' ### Indexing raster data


values <- p224r63_2011_B1[1,] # First row
head(values)

values <- p224r63_2011_B1[,2] # Second column
values <- p224r63_2011_B1[,7:300] # Columns 7 to 300
values <- p224r63_2011_B1[1,2] # Single pixel: 1st row & 2nd column


#' The four corner pixels could also be queried by first extracting the number of rows and columns:


nr <- nrow(p224r63_2011_B1)
nc <- ncol(p224r63_2011_B1)
p224r63_2011_B1[c(1, nr), c(1, nc)]


#' Apart from rows and columns, each pixel has a unique cell ID that is used when you query *without* a comma.
#' Query by cell ID indexing works with:


p224r63_2011_B1[c(3141, 5926)]


#' The syntax expands naturally to select *all* pixels in the same way you could address all values in a matrix


values <- p224r63_2011_B1[ , ]
# equivalent to
values <- p224r63_2011_B1[]


#' Next we move on to addressing layers of stacks or bricks, which is similar to a list. Layers are addressed with double square brackets, which return a raster object and not just the pixel values:


p224r63_2011[[2]]


#' Querying more than one layer is also done with double square brackets (numbers or layer names can be used) or for single layers with the $ sign:


subs <- p224r63_2011[[1:5]]
subs <- p224r63_2011[[c("B1_sre, B6_bt")]]
subs <- p224r63_2011$B1_sre


#' To remove layers, you can use negative indices and double square brackets or the `dropLayer()` command:


p224r63_2011.drop15 <- p224r63_2011[[-c(1,5)]]
# equivalen to
p224r63_2011.drop15 <- dropLayer(p224r63_2011, c(1,5))
# equivalent to
p224r63_2011.drop15 <- dropLayer(p224r63_2011, c("B1_sre", "B5_sre"))


#' We can also add more layers to our raster brick. For example, we can divide one band by two and add it to the raster brick with `addLayer()`, or specify a new layer name in square brackets:


newLayer <- p224r63_2011[[1]]/2
names(newLayer) <- "new"

p224r63_2011.add <- addLayer(p224r63_2011, newLayer)
# equivalent to
p224r63_2011.add[["myNewLayer"]] <- newLayer


#' Querying pixels of a stack or brick with single square brackets will return a matrix with one column per layer and pixels in rows:


values <- p224r63_2011[3:7, 84]


#' ### Indexing vector data

#' As described earlier SpatialPointsDataFrames carry an attribute table in the *data* slot.


csv.spdf@data


#' To access the coordinates:


coordinates(csv.spdf)


#' ### Adding information to vector data

#' Another important step that is often needed is to combine the features of separate spatial objects into a single object. To combine features, they must be of the same type, for example, you can only merge points with points and polygons with polygons.


pts_random <- spsample(studyArea, 5, type="random")
pts_combined <- spRbind(csv.spdf, pts_random)


#' When combining polygons unique IDs are mandatory because each entry in the attribute table must be referenced to many coordinates defining the shape of the polygon. Therefore, we have to change the polygon IDs manually before we can merge the polygons into one object.


studyArea2 <- studyArea
newNames <- paste0("p", 1:length(studyArea))
spChFIDs(studyArea2) <- newNames
studyArea_double <- spRbind(studyArea, studyArea2)


#' Alternatively, you can use the `rbind` command defined by the *sp* package, which assigns unique IDs automatically.


studyArea_triple <- rbind(studyArea, studyArea, studyArea2, makeUniqueIDs = TRUE)


#' ### Coordinate reference systems and data projections

#' Detailed information about any projection can be found at the Spatial Reference website (http://spatialreference.org).

#' Check if two objects have the same projection string:


identicalCRS(pa_utmS, roads)


#' ### Raster-Vector Conversions

#' The conversion from vector to raster data involveds a loss of accuracy and should only be used if truly necessary.


pa_utmN_rast <- rasterize(pa_utmN, p224r63_2011, field = "wpaid")


#' Now, it is likely that several vector features fall inside a target pixel. The default is to use only the last feature that was discovered.
#' However, we could also count the number of points falling into a pixel, to estimate the frequency of species observed in a given area.


xy <- coordinates(csv.spdf)
x <- rnorm(3000, mean = xy[, 1], sd = 3000)
y <- rnorm(3000, mean = xy[, 2, sd = 3000])
obs <- data.frame(X = x, Y = y)
coordinates(obs) <- c("X", "Y")

target <- raster(studyArea, resolution = 1000)
obs_frequency <- rasterize(obs, target, fun = "count")


#' The commands `rasterToPolygons()` and `rasterToPoints()` allow you to vectorize raster data. Perhaps the most relevant function worth mentioning is `rastertoContour()`, which creates isolines, typically from a DEM, but of course also from any other raster. Let us load the DEM from the example data to calculate isolines at 50 m intervales:


dem <- raster("raster_data/dem/srtm_dem.tif")
contourLines <- rasterToContour(dem, levels = seq(50, 350, 50))


#' ### Spatial Queries: Extracting Values

#' Let us assume we want to query the Landsat pixel values corresponding to our newly created *csv.spdf SpatialPointsDataFrame*. The quickest way to achieve this is to use single square brackets, which recognize that we providing *SpatialPoints* and search for the corresponding pixels to return their values:


ptVals <- p224r63_2011[csv.spdf]


#' A more obvious way of doing the same thing is to use `raster::extract()`; it has many additional options, such as buffers or summary statistics that are not accessible via the square bracket notation.


ptVals <- extract(p224r63_2011, y = csv.spdf)

pa_values <- extract(p224r63_2011[[1:4]], pa_utmN)


#' By default, the pixel centre determines whether a pixel is considered inside or ouside a polygon. Extracting by polygons will return a list of the same order as the polygons. Each list entry consists of a matrix with pixels in rows and layers in columns (unless you extract from a single layer, then it is just a vector). To look at the values of the Parakana PAs we could do the following:


para <- which(pa_utmN[["name"]] == "Parakan")
para_values <- pa_values[[para]]
head(para_values, 3)


#' Next we might wnt to calculate some summary statistics for each PA.


pa_values_mean <- extract(p224r63_2011[[1:4]], pa_utmN, fun = mean, na.rm = TRUE)


#' We can append the information to the data slot of the *SpatialPolygonsDataFrame*:


pa_utmnN@data <- cbind(pa_utmN@data, pa_values_mean)


#' Sometimes, extracting the single pixel behind a point might not be adequate because a point may not refer to this single location but can have a larger catchment area.


pts_sd <- extract(p224r63_2011[[3]], y = csv.spdf, buffer = 500)

library(reshape2)
pts_sd_df <- melt(pts_sd)
boxplot(value ~ L1, pts_sd_df)


#' Now we move on to use vectors to extract information from other vectors. For example, we can extract the PA *wdpaid* and *name* in which each point falls.


pa_per_point <- over(x = csv.spdf, y = pa_utmN)
pa_per_point[, c("wdpaid", "name")]


#' A common task is to count how many points fall within a polygon. For this we need to query the points per PA and then simply count(`length()`) the number of positive matches:


npts_per_pa <- over(x = pa_utmN, y = csv.spdf, fn = length)
data.frame(pa_utmN$name, npts_per_pa$id)


#' ### Computing distances and buffers


pa_utmN_dist <- distance(pa_utmN_rast)

csv.spdf_dist <- distanceFromPoints(p224r63_2011, csv.spdf)


#' The last important distance function is `raster::pointDistance()`, which calculates the pairwise distance matrix between spatial points either in spherical (*lonlat = TRUE*) or planar coordinates (*lonlat = FALSE*). 

#' An operation closely related to distance analysis is *buffering*, that is, enforcing thresholds on top of a distance calculation. For buffering we use `gBuffer()`from the *rgeos* package, which provides tools for almost any vector geometry operations you might be interested in.


csv.spdf_buf <- gBuffer(csv.spdf, width = 5000, byid = TRUE)

csv.spdf:buf2 <- gBuffer(csv.spdf, width = c(csv.spdf$value * 100), byid=TRUE, id = NULL)


#' Polygons can also be buffered from their boundaries towards the inside by specifying a negative radius. This could be useful if we wanted to plan a disturbance-free core zone for a PA:


pa_coreZone <- gBuffer(pa_utmN, width = -4000, byid= TRUE)


#' ### Raster processing: Cropping and Masking

#' #### Vector processing: Intersecting and Clipping


intersec_pa_poi <- gIntersection(pa_tmN, csv.spdf, byid=TRUE)

# or

pa_per_point <- over(csv.spdf, pa_utmN)
inside_pa <- !is.na(pa_per_point$wdpaid)
points_inside_pa <- csv.spdf[inside_pa,]

buf_inside_pa <- gIntersection(pa_utmN, csv.spdf_buf, byid=TRUE)


#' ### Raster Arithmetic

#' ### Managing Spatial Resolution


p224r63_2011_agg <- aggregate(p224r63_2011, fact = 2, fun = mean)
res(p224r63_2011_agg)

p224r63_2011_agg_dis <- disaggregate(p224r63_2011_agg, fact = 2)
res(p224r63_2011_agg_dis)

modis <- raster(studyArea, resolution = 250)
p224r63_2011_250m <- resample(p224r63_2011, modis, method = "bilinear")


#' ### Raster Tricks

#' #### Writing directly to file


calcResult <- calc(p224r63_2011, fun = function(x) {3 * x + 2 * x^2}, filename = "raster_data/p224r63_2011_mycalc.tif", datatype = "FLT4S")


#' #### Setting global options


rasterOptions(tmpdir = "folder/on/another/harddrive")
# Query current options
rasterOptions()


#' #### parallel computing

#' A few raster functions, such as *projectRaster()* or *resample()* will automatically run in parallel if you run the `beginCluster()` command in advance. Other functions, such as `calc()`, `overlay()`, or `predict()` can be run in parallel if run inside the `raster::clusteR()`function. Most functions in the RStoolbox package will also run in parallel if `beginCluster()`has been called beforehand.


beginCluster(3)

sqrt_1 <- calc(p224r63_2011, fun = sqrt)

sqrt_2 <- clusterR(x = p224r63_2011, fun = calc, args = list(fun = sqrt))

endCluster()


#' ## Pre-processing Remote Sensing Data


library(RStoolbox)
meta2011 <- readMeta("raster_data/LT52240632011210/LT52240632011210CUB01_MTL.txt")
summary(meta2011)

p224r63_2011 <- stackMeta(meta2011)
p224r63_2011


#' ### Data conversion: From Digital Numbers to Meaningful Units

#' To recalculate the actual radiation measured at the sensor, we need to apply sensor- and band-specific calibration parameters to the digital numbers (DNs). Typically these parameters consist of a multiplicative term (often referred to as *agin* or *mult*) and an additive term (*offset* or *add*), which are applied to each pixel as `gain*x+offset`. The bandwise conversion parameters can be extracted thus:


dn2rad <- meta2011$CALRAD
dn2rad


#' These gain/offset information can now be applied to the single bands:


p224r63_2011_rad <- p224r63_2011 * dn2rad$gain + dn2rad$offset


#' The same can be accomplished with the `RStoolbox::radCor()` command:


p224r63_2011_rad <- radCor(p224r63_2011, metaData = meta2011, method = "rad")


#' We have now converted the Dns to **top-of-atmosphere radiance** to physically meaningful units (W*m-2*srad-1*microm-1).

If we normalize the radiance to the approximate solar irradiation at a given day of the year, we receive the **top-of-atmosphere reflectance**, also referred to as **apparaent reflectance**, a unit-less qunatity with valid values ranging from 0 to 1.


p224r63_2011_rad <- radCor(p224r63_2011, metaData = meta2011, method = "apref")


#' In addition to calculating reflectance, the command `radCor()` will also convert the thermal bands to brightness temperature (in degrees Kelvin), which is the temperature a black body would have, given the measured radiation. Brightness temperature is an approximation of land surface temperature. For true land surface temperature, the brightness temperature must be adjusted for land cover emissivity, which is different for different land cover types, for example, water and forest.

#' ### Atmospheric correction

#' Different **dark object subtraction** (DOS) approaches exist. The simplest (*method = "sdos"*) requires one haze value per band, that is, the measured radiation in pixels that should have reflected almost no radiation, estimates haze contribution and subtracts it from all pixels in the image. To estimate the haze value, we can use `RStoolbox::estimateHaze`on the DN data by specifying the proportion of dark pixels we expect in the image (`darkProk`). Note that you should stick to bands in the visible and NIR wavelength range, as the effect of haze in wavelengths > 800 nm is negligible. The haze threshold is estimated using bands 1-4:


haze <- estimateHaze(p224r63_2011, darkProp = 0.01, hazeBand = c("B1_dn", "B2_dn", "B3_dn", "B4_dn"))


#' And bands 1-4 are corrected for haze:


p224r63_2011_sdos <- radCor(p224r63_2011, metaData=meta2011, hazeValues = haze, hazeBands = c("B1_dn", "B2_dn", "B3_dn", "B4_dn"), method = "sdos")


#' The more sophisticated `method = "dos"`, relies on an atmospheric scattering model that assumes scattering is highest in the blue wavelengths and decreases exponentially towards the NIR. By estimating haze in the blue band, we can also derive haze in all other bands. Note that unlike the simple DOS approach, this method requires sensor-specific parameters, that is, it must be explicitly implemented for the sensor you wish to correct:


p224r63_2011_sdos <- radCor(p224r63_2011, metaData=meta2011, darkProp = 0.01, method = "dos")


#' The DOS-based approaches are very simple and primarily address the effect of atmospheric haze. More complex alternatives have been developed, however, most of them are commercial. A free alternative we will not cover within the scope of this chaptr, yet worth mentionig, is the `i.atcorr` commnd in `GRASS` which relies on the MODIS 6S atmospheric correction routine.

#' Alternatively, use Landsat CDR imagery:


xml_meta2011 <- readMeta("raster_data/LT52240632011210/LT52240632011210CUB01.xml")
p224r63_2011_cdr <- stackMeta(xml_meta2011, quantity = c("sre", "bt"))


#' The original reflectance values (0-1) in the Landsat CDR product were scaled by a factor 10000 to store them as integers while maintaining a reasonable level of precision.


scaleF <- getMeta(p224r63_2011_cdr, xml_meta2011, what = "SCALE_FACTOR")
scaleF

p224r63_2011_cdr <- p224r63_2011_cdr * scaleF


#' ### Topographic illumination correction

#' In order to run, `topCor()` needs the solar angles of your scene, which are stored in the meta-data file. Using the `solarAngles` argument instead you could also provide solar azimuth and solar zenith manually.


dem <- raster("raster_data/dem/srtm_dem.tif")

p224r63_2011_cdr_ilu <- topCor(p224r63_2011_cdr, dem = dem, metaData = meta2011, method = "C")



#' A note of caution: many of the topographic illumination corrections have difficulty dealing with low sun elevation angles and slopes without direction insolation, which you will find in very rugged terrain. Also, these approaches only deal with terrain-based illumination. Shadows cast by mountains will not be corrected for. If these are a serious problem in your scene, you should look into the *GRASS* function `r.sunmask` which can calculate a cast-shadow mask that may be useful for reflectance reconstruction or temporal interpolation from alternative scenes.

#' ### Geolocation correction

#' #### Automated co-registration

#' Selecting ground control points (GCPs). Selecting GCPs can be very tedious and depending on the scene characteristics, for example, a feature-poor forested area, it can be difficult to find enough suitable GCPs. Automated image co-registration routines try to find the best shift of an image relative to a reference image by maximizing the similarity of two images. As a result, the assumption is that similarity is maximized once both images overlap perfectly. A common problem with most automated co-registration approaches is that they get stuck in local optima, that is, they falsely identify an optimal shift.
#' One of the best image matching approaches to be performed is based on the *mutual information* criterion from information theory.


p224r63_2011_shifted <- shift(p224r63_2011, 30, 60)

p224r63_2011_corre <- coregisterImages(p224r63_2011, p224r63_2011_shifted, verbose=TRUE)


#' ### Masking invalid values

#' #### Clouds and cloud shadows

#' A typical property distinguishing clouds from land surface is very low thermal radiation and high reflectance in the visible region of the spectrum.


xml_meta88 <- readMeta("raster_data/LT52240631988227/LT52240631988227CUB02.xml")

p224r63_1988 <- stackMeta(xml_meta88)

ex <- extent(c(625500, 635300, -493500, -48500))

p224r63_1988_crop <- crop(p224r63_1988, ex)


#' Identify cloud pixels in a semi-automated fashion


cmask <- cloudMask(p224r63_1988_crop, treshold = 0.8, buffer = 3, blue = "B1_Sre", tir = "B6_bt")

shadow <- cloudShadowMask(p224r63_1988_crop, cmask)


#' This approach is useful for quick, interactive cloud masking, however, it is not as effective in areas with snow or rugged terrain. Automated cloud detection algorithms are usually much more complex in that they use many additional empirical tresholds and band ratios. These are computationally intensive but do no require user interaction. The Landsat CDR data ship with two such cloud products: the internal LEDAPS cloud masks and the product of the *fmask* algorithm (Zhu et al. 2015).


p224r63_1988_qa <- stackMeta(xml_meta88, category = "qa")
p224r63_1988_crop_qa <- crop(p224r63_1988_qa, ex)

# Then the cloud and cloud shadow layers are summed up
ccs_mask <- sum(p224r63_1988_crop_qa[[c("QA_cloud_shadow", "QA_cloud", "QA_adjacent_cloud")]])

# Compare the layers
plot(stack(p224r63_1988_crop_qa[[c("QA_cloud_shadow", "QA_cloud", "QA_adjacent_cloud")]], clouds))


#' The fmask is encoded into five integers: land surface (0), water (1), cloud shadow (2), snow (3) and clouds (4). To create the mask, we select all pixels of values 2 and 4 and replace them with NA.


cfmask <- p224r63_1988_crop_qa[["QA_cfmask"]]
cfmask[cfmask %in% c(2,4)] <- NA

p224r63_1988_crop <- mask(p224r63_1988_crop, ccs_mask, maskvalue=NA)


#' ### Guidelines: When to Correct for What?

#' It is imperative that your data have a precise geolocation. Whenever you combine spatial data you must ensure they spatially fit to another, and if necessary correct it. Topographic illumination correction should be applied if your scene contains pronounced terrain aind you notice illumination differences solely due to terrain orientation. For single-scene classification, it is generally acceptable just using the DNs. As soon as you plan to compare spectral values between scenes or time steps, you should at least use top-of-atmosphere reflectance correction and, where feasible, atmospherically correct to at-surface reflectance data. When combining imagery of different sensors, especially those which have different radiometric resolutions, you should, at the minimum, convert DNs to radiance values. However, as this often means that you will have multi-date imagery, you should also use reflectance data here also. The calculation of spectral indices also benefits from atmospherically corrected data.

#' ### Quality layers

#' Landsat 8 - http://landsat.usgs.gov/L8QualityAssessmentBand.php

#' Landsat L-LDOPE Toolbelt (http://landsat.usgs.gov/L-LDOPE_Toolbelt.php)

#' In R, you can use `encodeQA()` and its reverse counterpart `decodeQA` to calculate QA values by narrowing down which pixel properties you want, i.e. extract cirrus clouds.


justClouds <- encodeQA(cloud = "high", cirrus = "high", snow = "low", water = "low")
justClouds


#' The binary form is provided by:


decodeQA(justClouds)


#' The USGS provides the same QA access tools for MODIS under LDOPE Tools (http://lpdaac.usgs.gov/tools/ldope_tools)

#' ### Artefacts

#' Occassionally you will find pixels with negative radiance. This does not make sense of course, as radiance is an absolute physical quantity. Thus, it is a valid approach to set these pixels to 0.


p224r63_2011_rad_bad <- p224r63_2011_rad < 0
cellStats(p224r63_2011_rad_bad, sum)
plot(p224r63_2011_rad_bad)

p224r63_2011_rad[p224r63_2011_rad_bad] <- 0


#' ### Image matching and mosaicking


p224r63_2 <- brick("raster_data/mosaic/LT2240632011178CUB02.grd")

layers <- names(p224r63_adjacent)
ls_multiscene <- merge(p224r63_adjacent, p224r63_2011[[layers]])

ls_multiscene_m <- mosaic(p224r63_adjacent, p224r63_2011[[layers]], fun=mean)


#' #### Contrast matching


cm <- raster("raster_data/mosaic/LT52240632011178CUB02_cloudMask.grd")
p224r63_adjacent <- mask(p224r63_adjacent, cm)

p224r63_adjacent_hm <- histMatch(p224r63_adjacent, p224r63_2011[[layers]])

ls_multiscene_hm <- merge(p224r63_adjacent_hm, p224r63_2011[[layers]])


#' #### Pan-sharpening


pan <- raster("raster_data/pan/LC82240632013167LGN00_B8.tif")

p224r63_2011_sub <- crop(p224r63_2011_cdr, pan)

# We then apply the IHS method and Brovey transformation

ihsSharp <- panSharpen(p224r63_2011_sub, pan = pan, method = "ihs", r = 3, g = 2, b = 1)

broSharp <- panSharpen(p224r63_2011_sub, pan = pan, method = "brovey", r = 3, g = 2, b = 1)

res(p224r63_2011_sub)
res(ihrSharp)


#' ## Field Data for Remote Sensing Data Analysis

#' With the generic function `spsample()` we can easily sample point locations within a spatial object, which can be either a simple square area, a grid, polygons or spatial lines. different sample patterns can be set (random, regular, stratified, non-aligned, hexagonal, clustered, Fibonacci).


set.seed(1)
RandomPoints <- spsample(x = SpatialObject, n = 100, type = "random")

study_area <- readOGR("vector_data/", "study_area_UTMnorth")
ranPoi_random <- spsample(study_area, n = 100, type = "random")
ranPoi_regular  <- spsample(study_area, n = 100, type = "regular")
ranPoi_nonaligned  <- spsample(study_area, n = 100, type = "nonaligned")
ranPoi_stratified  <- spsample(study_area, n = 100, type = "stratified")


#' Sampling random points on top of a raster file is feasible with the `randomSample()` command:


library(dismo)
ranPoi_lsat <- randomPoints(p224r63_2011, 100)

library(rgdal)
ranPoi_lsat.sp <- SpatialPoints(ranPoi_lsat)

projection(ranPoi_lsat.sp) <- projection(p224r63_2011)


#' #### Constraining the sampling


road_dist <- raster("raster_data/road_dist.tif")

road_Dist[road_dist >= 2000] <- NA
sample_road <- randomPoints(road_dist, 100)


#' ## From Spectral to Ecological Information

#' ### Vegetation indices

#' NDVI

#' ### Comparison of vegetation indices


p224r63_2011m_viStack <- spectralIndices(p224r63_2011m, red = "B3_sre", nir = "B4_sre", indices = c("NDVI", "MSAVI2", "DVI"))

p224r63_2011_viStack_sd <- calc(p224r63_2011m_viStack, fun=sd)
p224r63_2011_NDVI_MSAVI_cor <- corLocal(p224r63_2011m_viStack$NDVI, p224r63_2011m_viStack$DVI, ngb=11, method="spearman")


#' ### Linear transformations: principal components


p224r63_2011_pca <- rasterPCA(p224r63_2011)
summary(p224r63_2011_pca$model)

loadings(p224r63_2011_pca$model)


#' Another prominent spectral transformation is the *tasseled-cap* transformation (TCT), which was originally developed for Landsat data. The primary result of the TCT are three components: brightness, greenness and soil wetness.


tctBands <- c("B1_sre", "B2_sre", "B3_sre", "B4_sre", "B5_sre", "B7_sre")
p224r63_TCT <- tasseledCap(p224r63_2011[[tctBands]], sat = "Landsat5TM")


#' ## Land Cover or Image Classification Approaches

#' ### Practical examples of land cover classifications

#' #### Running an Unsupervised Classification


set.seed(6)
p224r63_2011_uc <- unsuperClass(p224r63_2011, nClasses=4, nStarts=50, nSamples=10000, norm = TRUE)
plot(p224r63_2011_uc$map)


#' or manually


p224r63_2011_sub <- p224r63_2011[[c(1:4)]]

p224r63_2011.kmeans <- kmeans(p224r63_2011_sub[], centers = 5, iter.max = 100, nstart = 100)

p224r63_2011_uc <- raster(p224r63_2011)
p224r63_2011_uc[] <- p224r63_2011.kmeans$cluster


#' #### Running a supervised classification

#' - Collect training and validation data


training_2011 <- readOGR("vector_data/", "training_2011")
validation_2011 <- readOGR("vector_data/", "validation_2011")

p224r63_2011_sc <- superClass(img = p224r63_2011, model = "rf", trainData = training_2011, responseCol = "class_name")
p224r63_2011_sc


#' or manually


uniqueClasses <- unique(training_2011$class_name)
uniqueClasses

set.seed(25)
for(i in 1:length(uniqueClasses)){
  # Select only polygons of current class
  class_data <- subset(training_2011, class_name == uniqueClasses[i])
  # Get random points for these polygons
  classpts <- spsample(class_data, type = "random", n = 100)
  # Add class column to SpatialPoints object
  classpts$class <- rep(uniqueClasses[i], length(classpts))
  if (i == 1) {
    xy <- classpts
  } else{
    xy <- rbind(xy, classpts)
  }
}

plot(p224r63_2011, 4)
points(xy)

trainvals <- extract(p224r63_2011, y = xy, cellnumbers = TRUE)
trainvals <- data.frame(response = xy$class, trainvals)

head(trainvals[,1:5], 3)

any(duplicated(trainvals$cells))
trainvals <- trainvals[!duplicated(trainvals$cells), -2]

names(trainvals)

table(Trainvals$response)

model_1 <- randomForest(response ~ . , data = trainvals, na.action = na.omit, confusion = TRUE)

p224r63_2011_sc <- predict(p224r63_2011, model_1, filename="results/p224r63_2011_sc.tif", format="GTiff", datatype="INT1U", type="response", overwrite=TRUE)


#' Get area covered by classes


p224r63_2011_sc.freq <- freq(p224r63_2011_sc, useNA = "no")

resLsat <- res(p224r63_2011_sc)
area_km2 <- p224r63_2011_sc.freq[,"count"] * prod(resLsat) * 1e-06
data.frame(landcover = uniqueClasses, area_km2 = area_km2)


#' Look at the internal tree of a fitted random forest model


library(party)
ctree_model <- ctree(response ~ ., data = trainvals)
plot(ctree_model)


#' #### Accuracy assessment in practice


p224r63_2011_sc <- superClass(img = p224r63_2011, model = "rf", trainData = training_2011, responseCol = "class_name", valData = validation_2011)

p224r63_2011_sc <- superClass(img = p224r63_2011, model = "rf", trainData = training_2011, responseCol = "class_name", trainPartition = 0.7)


#' Manual validation


set.seed(7)
xy_val <- lapply(uniqueClasses, function(class_i){
  class_data <- subset(training_2011, class_name == class_i)
  classpts <- spsample(class_data, type = "random", n = 100)
  classpts$class <- rep(class_i, length(classpts))
  return(classpts)
})

xy_val <- do.call("rbind", xy.val)

pred <- extract(p224r63_2011_sc$map, xy_val, cellnumbers = TRUE)

dup <- duplicated(pred)
pred <- pred[!dup, "class_name"]
obs <- xy_val$class[!dup]

valFactor <- uniqueClasses[pred]
head(valFactor)

confusionMatrix(obs, reference = valFactor)


#' #### Classification uncertainty maps


# Class probabilities p
p <- predict(p224r63_2011, model_1, type = "vote", index = 1:3)

pLog <- log(p)
pLog[pLog == -Inf] <- 0
forestAgreement <- -1 * sum(p* pLog)

groupIndex <- rep(1:5, length.out = nrow(trainvals))
groupIndex <- sample(groupIndex)

predictions <- lapply(1:5, function(i){
  model <- randomForest(response ~., data = trainvals[groupIndex != i,])
  predict(p224r63_2011, model)
})

prediction_stack <- stack(predictions)

nClasses <- calc(prediction_stack, function(x){
  occuringClasses <- unique(x)
  length(!is.na(occuringClasses))
})


#' Compare three models


models <- c("rf", "svmRadial", "pls")
ensemble <- lapply(models, function(mod){
  set.seed(5)
  sc <- superClass(p224r63_2011, trainData = training_2011,
                   responseCol = "class_name", model = mod)
  return(sc$map)
})

prediction_stack <- stack(ensemble)
names(prediction_stack) <- models

modelEntropy <- rasterEntropy(prediction_stack)


#' ## Land cover change or change detection 

#' ### Common approaches to change detection:

#' #### Post-classification comparison

#' #### spectral change vector analysis


tc2011 <- tasseledCap(dropLayer(p224r63_2011, "B6_bt"), sat = "Landsat5TM")
tc1988 <- tasseledCap(dropLayer(p224r63_1988, "B6_bt"), sat = "Landsat5TM")

changeVec <- rasterCVA(tc2011[[1:2]], tc1988[[1:2]])


#' #### multi-date classification


p224r63_88_11 <- stack(p224r63_1988, p224r63_2011)

change_classes <- readOGR("vector_data/", "change_classes_1988_2011")

p224r63_88_11_change <- superClass(img = p224r63_88_11, nSamples = 1000, trainData = change_classes, responseCol = "class")

# Display class <- ID mapping
p224r63_88_11_change$classMapping
plot(p224r63_88_11change$map)


#' ### Landscape Statistics


p224r63_88_11_change.freq <- freq(p224r63_88_11_change$map, useNA = "no")
p224r63_88_11_change.freq[,"class"] * 30^2 * 1e-06


#' ### Alternative Online Resources


study_area_ll <- readOGR("vector_data/", "study_area_ll")

tiles <- calc_gfc_tiles(study_area_ll)
length(tiles)

plot(tiles)
plot(study_area_ll, add=TRUE, lty=2, col="blue")

dir.create("raster_data/GFC")

download_tiles(tiles, "raster_data/GFC", first_and_last=FALSE)

gfc_data <- extract_gfc(study_area_ll, data_folder = "raster_data/GFC", filename = "raster_Data/GFC/GFC_studyArea.tif")

gfc_data <- mask(gfc_data, gfc_data$datamask, maskvalue=2)
ggRGB(gfc_data, r = "loss", g = "gain", b = NULL)

lf <- freq(gfc_data$lossyear, useNA = "no")

barplot(lf[-1, "count"], names.arg = 2001:2012)


#' ## Continuous Land Cover Information

#' ### Existing Data Sets

#' - Tree Cover Continuous Fields from the Global Land Cover Facility (http://landcover.org/data/treecover)
#' - Landsat Treef Cover Continuous Field from the Global Land Cover Facility  (http.//landcover.org/data/landsatTreecover)
#' - MODIS Land Vegetation Continuous Fields (http://modis-land.gsfc.nasa.gov/vcc.html)
#' - Fraction of Green Vegetation Cover from Copernicus (http://land.copernicus.eu/global/?q=products/fCover)

#' ### Fraction of vegetation cover mapping


p224r63_2011_sc <- raster("raster_data/products/p224r63_2011_sc.tif")

# Assemble filenames for bands 1 to 7 (alternative: list.files())
filenames <- paste0("raster_data/MODIS/", "MCD43A4.MRTWEB.A2011193.005.Nadir_Reflectance_Band", 1:7, ".tif")

filenames[1:3]

# Import stack files
mcd43a4_193_2011 <- stack(filenames)

# Make pretty names
names(mcd43a4_193_2011) <- paste0("b", 1:7)

mcd43a4_193_2011

plot(mcd43a4_193_2011[[6]])
studyArea_sinu_modis <- drawPoly()

projection(studyArea_sinu_modis) <- projection(mcd43a4_193_2011)

studyArea_sinu_modis <- as(study_area_extent, "SpatialPolygonsDataFrame")
writeOGR(studyArea_sinu_modis, "vector_data/", "study_area_extent_sinu", driver = "ESRI Shapefile")

mcd43a4_193_2011 <- crop(mcd43a4_193_2011, studyArea_sinu_modis)

mcd43a4_193_2011 <- projectRaster(mcd43a4_193_2011, method = "ngb", crs = crs(p224r63_2011_sc), filename= "raster_data/MODIS/mcd43a4_193_2011_UTMwgs84_crop.tif")

fCover_result <- fCover(classImage = p224r63_2011_sc, model = "rf", predImage=mcd43a4_193_2011, classes=c(1:3))
fCover_result$map


#' or manually


commonExt <- intersect(extent(mcd43a4_193_2011), extent(p224r63_2011_sc))
modisSampleCells <- sampleRandom(mcd43a4_193_2011, ext = commonExt, size = 1000, xy = TRUE, na.rm = TRUE)

# Which classes occur?
uniqueClasses <- unique(p224r63_2011_sc)

fc_Samples <- matrix(nrow = nrow(modisSampleCells), ncol = length(uniqueClasses) + 1)
colnames(fc_samples) <- c(paste0("Class", uniqueClasses), "NAs")

# MODIS pixel size
modRes <- res(mcd43a4_193_2011)

# Loop through all sampled MODIS pixels and calculate fractional cover
for(i in 1:nrow(modisSampleCells)){
  # Get the x and y coordinate from the center of the MODIS pixel
  centerCoords <- modisSampleCells[i, c("x", "y")]
  
  # Calculate the extent of the selected MODIS pixel
  ext <- extent(centerCoords[1] - modRes[1]/2,
                centerCoords[1] + modRes[1]/2,
                centerCoords[2] - modRes[2]/2,
                centerCoords[2] + modRes[2]/2)
  
  # Retrieve all classified image pixels that fall inside the MODIS pixel
  classCellValues <- extract(p224r63_2011_sc, ext)
  
  # Calculate the frequency of each class
  counts <- table(factor(classCellValues, levels = uniqueClasses), 
                  useNA = "always")
  fc_samples[i, ] <- counts/sum(counts)
}

head(fc_samples, 3)

# Drop x and y columns from sample cells
modisSampleCells <- modisSampleCells[, -c(1,2)]

# Create a training data set
trainingData <- data.frame(response = fc_samples[, "Class1"], modisSampleCells)

# Include only samples with complete data
trainingData <- trainingData[fc_samples[, "NAs"] == 0, ]

model_class1 <- randomForest(response ~ ., data = trainingData)
model_class1

fCover_class1 <- predict(mcd43a4_193_2011, model_class1)

# Enforce [0,1] value range
fCover_class1 <- clamp(fCover_class1, lower = 0, upper = 1, filename="raster_data/results/fCover_mcd_lsat_class1.tif")


#' Define new classes or maps


# Make a copy of the first layer of the fCover() result
fCover_90 <- fCover_result$map[[1]]

# Binary mask (0 = fCover < 90%, 1 = fCover > 90%)
fCover_90_mask <- fCOver_90 > 0.9

# Alternative continuous mask
fCover_90[fCover_90 > 0.9] <- NA

# Define class intervals
from <- c(0, 0.25, 0.5, 0.9)
to <- c(0.25, 0.5, 0.9, 1)
becomes <- c(1,2,3,4)
classDef <- cbind(from, to, becomes)
classDef

newClasses <- reclassify(fCover_result$map[[1]], rcl= classDef)

ecotoneClass <- fCover_result$map[[1]] > 0.5 & fCover_result$map[[2]] > 0.2


#' ## Time Series Analysis


install.packages("MODIS", repos = "http://R-Forge.R-project.org")
library(MODIS)

MODISoptions(localArcPath = getwd(), outDirPath = getwd())

viname <- "ndvi"
product <- "MOD13Q1"
ofilename <- paste0(product, "_", viname, "brick.grd")

pth <- paste0(getwd(), "raster_data/", product)
fileout <- paste(pth, "/", ofilename, sep = "")

if(!file.exists(fileout)){
  if(!file.exists(pth)){
    print("the outfolder does not exist and will be created")
    print(pth)
    dir.create(pth)
  }
}

tileH <- 13
tileV <- 9
begin <- "2002.01.01"
end <- "2014.01.31"

modis.hdf <- getHdf(product = product, begin = begin, end = end, tileH = tileH, tileV = tileV, checkIntegrity=TRUE)

studyarea <- readOGR("vector_data", "study_area_extent_sinu")

for(i in 1:length(modis.hdf[[1]])){
  ifl <- unlist(strsplit(unlist(strsplit(modis.hdf[[1]][i],
                                         c("[/]")))[5], "[.]"))[1]
  print(ifl)
  fn <- paste("cvi_", ifl, ".tif", sep="")
  
  if(is.na(ifl) | file.exists(fn)){
    print("file exists or is not available on the server")
  }else{
   sds <- get_subdatasets(modis.hdf[[1]][i]) 
  }
}

#Extract NDVI
tmp <- rasterTmpFile()
extension(tmp) <- "tif"
gdal_translate(sds[1], dst_dataset = tmp)
ndvi <- crop(x = raster(tmp), y = studyarea)/(10000^2)

#Extract MODIS quality layer
tmp2 <- rasterTmpFile()
extension(tmp2) <- "tif"
gdal_translate(sds[12], dst_dataset = tmp)
rel <- crop(x = raster(tmp2), y = studyarea, filename = paste("rel", ifl, ".tif", sep = ""), dataType = "INT2U", formate = "GTiff")

# Function to discard NDVI time series pixels that do not fulfil a certain minimum criterion for quality
f_cleanvi <- function(vi, rel){
  i <- (rel <=1)
  res <- matrix(NA, length(i), 1)
  if(sum(i, na.rm = TRUE) > 0){
    i <- which(i)
    res[i] <- vi[i]
  }
  res
}

clvi <- overlay(ndvi, rel, fun = f_cleanvi, filename = fn, dataType = "INT2U", overwrite=TRUE)
rm(ndvi, clvi, rel, tmp, tmp2)


ofl <- list.files(pattern = glob2rx("rel*.tif"))
s <- do.call("brick", lapply(ofl, raster))
names(s) <- as.Date(modis.hdf[[1]], "A%y%j")

modis <- writeRaster(crop(s, studyarea), bandorder = "BIL", filename = "cvi_MOD13Q1_ndvi_brick.grd", options = c("COMPRESS=NONE"), dataType = "INT2U")


#' ### Change Detection


modis <- brick("raster_data/MODIS/NDVIMOD13Q1/cvi_MOD13Q1_ndvi_brick.grd")
modis.raw <- brick("raster_data/MODIS/NDVIMOD13Q1/ndvi_MOD13Q1_ndvi_brick.gri")

modis.2002 <- subset(modis, grep(paste("MOD13Q1A2002", sep=""), names(modis)))

# Median
modis2002.median <- calc(modis.2002, fun = function(x) median(x, na.rm=TRUE))

# Variance
modis2002.var <- calc(modis.2002, fun = function(x) var(x, na.rm = TRUE))

plot(modis2002.median/10000)


#' #### Summarizing temporal data


modis <- brick("raster_data/MODIS/NDVIMOD13Q1/cvi_MOD13Q1_ndvi_brick.grd")

modis
summary(modis)

# Function to convert day of year to an actual date
doy2date <- function(year, doy){
  as.Date(doy -1, origin = paste=(year, "-01-01"))
}

data.dates <- as.character(doy2date(2002, doy = as.numeric(substr(names(modis.2002), 18, 20))))

modis.stack <- stack()

# Calculate monthly MODIS composites from 16-day composite
for(i in 1:12){
  bands_i <- which(sapply(strsplit(data.dates, "-"),
                          function(x) as.numeric(x[2])) == i)
  modis.stack <- stack(modis.stack, calc(subset(modis.2002, bands_i)))
}

names(modis.stack) <- month.abb

library(rasterVis)
levelplot(modis.stack, par.settings = RdBuTheme)
bwplot(modis.stack)


#' #### Change Detection using Modis Satellite Images


library(bfast) # change detection
library(zoo) # time series handling
library(strucchange) # break detection

plot(modis, 1)
ndvi <- as.numeric(click(modis, n = 1))/10000 # now click on the map

ndvi <- ts(ndvi, start = c(2000, 4), frequency = 23)
plot(ndvi, ylab="NDVI", type="b")


#' ##### Change analysis


bfastmonitor(ndvi, formula = response ~ trend + harmon, start = c(2010, 1))
plot(bfm, ylab = "NDVI")

# Extract break point and magnitude for the single pixel time series
bfm$breakpoint
bfm$magnitude


#' Apply `bfastmonitor()` across the whole study area


# Define helper function
f_bfm <- function(x) {
  x <- ts(x, start = c(2000, 4), frequency = 23)/10000
  bfm <- bfastmonitor(data = x, start = c(2010, 1))
  return(cbind(bfm$breakpoint, bfm$magnitude))
}

rbfm <- calc(modis, fun = function(x) {t(apply(x, 1, f_bfm))})

timeofbreak <- raster(rbfm, 1)
magnitudeofbreak <- raster(rbfm, 2)
# only visualise the magnitude of the detected breaks
magnitudeofbreak[is.na(timeofbreak)] <- NA

plot(magnitudeofbreak, zlim = c(-0.2, 0.2), col = rev(diverge_hcl(7)))
ndvi <- ts(x, start = c(2000, 4), frequency = 23)
bfm <- bfastmonitor(ndvi, start = c(2010, 1))
plot(bfm, ylab = "NDVI")


#' #### Comparison of different cover change information


GFC <- brick("raster_data/GFC/GFC_studyArea.tif")
GFCsinu <- projectRaster(GFC, modis)modis.stack
GFC.bfast <- zonal(GFcsinu[[1]], magnitudeofbreak)


#' r.series in GRASS7
#' MODISTools Package
#' bfast and greenbrown package
#' spacetime-vis package
#' bfastSpatial - https://dutri001.github.io/bfastSpatial
#' bfastPlot - https://github.cmo/bendv/bfastPlot
#' timeSyncR - https://github.com/bendv/timeSyncR

#' ## Spatial Land Cover Pattern Analysis

#' ### Applying moving window techniques


p224r63_2011m.ndvi <- raster("raster_data/products/p224r63_2011m_ndvi.tif")

window <- matrix(1, nrow = 3, ncol = 3)
window

# Calculate variance
ndvi2011_3x3_var <- focal(p224r63_2011m.ndvi, w = window, fun = var)

#Calculate sd with various window sizes
focals <- lapply(c(3, 7, 11, 15), function(winDim){
  # Define window
  window <- matrix(data = 1, ncol = winDim, nrow = winDim)
  
  # Calculate standard deviation
  focal(p224r63_2011m.ndvi, w = window, fun = sd)
})

# Merge list of layers into single RasterStack
focal_sd_winSizes <- stack(focals)
names(focal_sd_winSizes) <- c("w3x3", "w7x7", "w11x11", "w15x15")

focal_diff <- focal_sd_winSizes[["w3x3"]] - focal_sd_winSizes[["w15x15"]]
plotRGB(focal_var_winSizes, r = 1, g = 2, b = 3, stretch = "lin")


#' #### Other Statistics


p224r63_2011_sc <- raster("raster_data/products/p224r63_2011_sc.tif")

#Count number of classes per 7 x 7 neighbourhood
countUnique <- function(x){
  if (anyNA(x)){
    return(NA) # we don't want to count NA as a class
  }else{
    return(length(unique(x)))
  }
}

p224r63_2011_sc_clean <- focal(p224r63_2011_sc, w = matrix(1,7,7), fun = countUnique)


#' #### Focal as a smoothing filter


window <- matrix(1/15^2, ncol=15, nrow = 15)
sum(window)

# Mean
ndvi_15x15_mean <- focal(p224r63_2011m.ndvi, w = window)

#Minimum
window <- matrix(1, ncol=15, nrow=15)
sum(window)
ndvi_15x15_min <- focal(p224r63_2011m.ndvi, w = window, fun=min)


#' #### Post-classification filtering


window(1, 3, 3)
p224r63_2011_sc_majority <- focal(p224r63_2011_sc, w = window, fun = modal)


#' #### Adding spatial weigths


fw_circle <- focalWeight(p224r63_2011m.ndvi, d = 160, type = "circle")

fw_circle[fw_circle > 0] <- 1
ndvi_circ160_var <- focal(p224r63_2011m.ndvi, w = fw_circle, fun = var)

# Used for smoothing
fw_gauss <- focalWeight(p224r63_2011m.ndvi, 20, type = "Gauss")
ndvi_sd20_gauss <- focal(p224r63_2011m.ndvi, w = fw_gauss)


#' #### Texture Based on the Grey Level Co-Occurrence Matrix


library(glcm)

ndvi_5x5_glcm <- glcm(p224r63_2011m.ndvi, window=c(5,5))


#' ### Combination with other approaches


training_2011 <- readOGR("vector_data/", "training_2011")
validation_2011 <- readOGR("vector_data/", "validation_2011")
p224r63_2011 <- brick("raster_data/final/p224r63_2011.grd")

set.seed(1)
p224r63_2011_sc <- superClass(p224r63_2011m, trainData = training_2011, valData = validation_2011, responseCol = "class_name")

ndvi_5x5_glcm <- glcm(p224r63_2011.ndvi, window=c(5,5))
p224r63_2011_plusText <- stack(p224r63_2011, ndvi_5x5_glcm$glcm_variance)

set.seed(1)
p224r63_2011_plusText_sc <- superClass(p224r63_2011_plusText, trainData = training_2011, valData = validation_2011, responseCol = "class_name")

# Compare validation
p224r63_2011_sc$validation$performance
p224r63_2011__plusText_sc$validation$performance


#' Also see r.li and r.pi modules in GRASS

#' ## Modelling Species Distributions

#' ### Data Pre-Processing

#' #### Combining Species and Environmental Information


pa_points <- readOGR("vector_data/", layer="occurrence")

env <- extract(p224r63_2011, pa_points)

pa_data <- data.frame(env, occurrence = pa_points$occurrence)


#' #### Visualization


plotRGB(p224r63_2011, stretch="lin")
pointSize <- pa_data$occurrence * 5 + 1
points(pa_points, cex = pointSize, pch=20)

boxplot(pa_data[,1:6])

boxplot(subset(data.frame(pa_data), occurrence == 0, select = c(1:6)))
boxplot(subset(data.frame(pa_data), occurrence == 1, select = c(1:6)))


#' #### Collinearity


library(car)
scatterplotMatrix(~. | occurrence, data = pa_Data, smoother = FALSE, reg.line = FALSE)

# Add correlation coefficients function
panel.cor <- function(x, y, digits = 2, cor_thresh = 0.7, col = c("black", "red"), ...){
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0,1,0,1))
  r <- abs(cor(x,y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  text(0.5, 0.5, txt, col = ifelse(r < cor_thresh, col[1], col[2]))
}


scatterplotMatrix(~. | occurrence, data = pa_data, smoother = FALSE, reg.line = FALSE, upper.panel = panel.cor)

pa_data_sel <- data.frame(pa_data[, c("occurrence", "B3_sre", "B4_sre")])


#' ### Modelling

#' #### Model fitting


library(randomForest)
library(dismo)
library(rms)
library(mgcv)

pa_data_sel$species <- factor(pa_data_sel$occurrence)

set.seed(2)
rfmodel <- randomForest(species ~ B3_sre + B4_sre, data = pa_data_sel)


#' #### Performance indices


round(val.prob(predict(rfmodel, newdata=pa_data_sel, type="prob")[,2], pa_data_sel[["occurrence"]])[c("C (ROC)", "R2", "Slope")], 2)


#' #### Cross-validation


set.seed(2)
fold <- kfold(pa_data_sel, k=5, by = pa_data_sel$species)

pa_data_sel$cv_pred <- NA
for(i in unique(fold)){
  traindata <- pa_data_sel[fold != i,]
  testdata <- pa_data_sel[fold == i,]
  cv_rfmodel <- RandomForest(species ~ B3_sre + B4_Sre, data =traindata, ntree = 500, nodesize = 1)
  pa_data_sel$cv_pred[fold == i] <- predict(cv_rfmodel, newdata = testdata, type = "prob")[,2]
}

round(val.prob(pa_data_sel[["cv_pred"]], pa_data_sel[["occurrence"]], p1 = false)[c("C (ROC)", "R2", "Slope")], 2)


#' #### Interpretation

#' use full model


rfmodel <- randomForest(species ~ B3_sre + B4_sre, data = pa_data_sel, ntree = 2000, nodesize = 5)


#' #### Spatial Patterns of Predictions


rfmap <- predict(p224r63_2011, rfmodel, type="prob", index=2)



#' Uncertainty of model predictions


fold_values <- sort(unique(fold))
set.seed(2)
models <- list()
for(i in fold_values){
  models[[i]] <- randomForest(species ~ B3_sre + B4_sre, data = pa_data_sel[fold != i,], ntree = 2000, nodesize = 5)
}
predictions <- lapply(models, function(x){
  predict(p224r63_2011, x, type = "prob", index = 2)
})

prediction_stack <- stack(predictions)
prediction_sd <- calc(prediction_stack, sd)


#' #### Predicting new environmental conditions

#' Multivariate Environmental Similarity Surfaces (MESS)


mess_Result <- mess(subset(p224r63_2011, c("B3_sre", "B4_sre")), pa_data_sel[, c("B3_sre", "B4_sre")])


#' #### Variable Importance


rfmodel <- randomForest(species ~ B3_sre + B4_sre, data = pa_data_sel, ntree = 2000, nodesize = 5, importance = TRUE)
rfImp <- importance(rfmodel)
varImpPlot(rfmodel, type = 1, scale = FALSE)


#' #### Response functions


par(mfrow = c(1,2))
partialPlot(rfmodel, traindata, B3_sre, xlab = "band 3", which.class = 1)
partialPlot(rfmodel, traindata, B4_sre, xlab = "band 4", which.class = 1)



#' #### GAM


AIC()

gamimp <- varImpBiomod()

plot(gammodel, page = 1, seWithMean = TRUE)


#' #### MaxEnt (Presence-only Data)


response(... , expand = 0)


#' ## Introduction to the Added Value of Animal Movement

#' ### Data inspection


move_csv <- read.csv("vector_data/move_data.csv", as.is=TRUE)

# We also want to load our Landsat scene
p224r63_2011m <- brick("raster_data/final/p224r63_masked.grd")

head(move_csv)

coordinates(move_sp) <- c("x", "y")
projection(move_sp) <- projection(p224r63_2011m)

library(move)
moveobj <- move(
  # set columns for coordinates
  x = move_csv$x, y = move_csv$y,
  # conversion of Date and Time from character to POSIX object
  time = as.POSIXct(move_csv$DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC"),
  # set projection
  proj = "+proj=utm +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",
  # set id
  animal = move_csv$id,
  # define data to use and attach to the spatial object
  data = move_csv)

moveobj

projection(moveobj)

plot(p224r63_2011m.ndvi)
lines(moveobj)
points(moveobj)

saveRDS(moveobj, "vector_data/movement_vector_w_time.rds")
moveobj <- readRDS("vector_data/movement_vector_w_time.rds")


#' ###First statistics


# Calls the functions distance(), time(), speed() and angle()
summary(moveobj)

# Calculates time lag between locations
timeLag(moveobj)

# Calculates the total number of locations
n.locs(moveobj)

# Summary statistics for both angle and speed
angleSummary(moveobj)
speedSummary(moveobj)

speeds <- speed(moveobj)
hist(speeds, breaks=seq(0,5,0.05))

moveobj.angle <- angle(moveobj)

library(circular)
plot(as.circular(moveobj.angle, type="directions", units="degrees", template="geographics", modulo="asis", rotation="clock", zero=0), stack=TRUE)

hist(speeds, breaks= "FD", main=NA, xlab="Speed in m/s")
hist(speeds[speeds > 0.05 & speeds < 3], breaks="FD", xlab="Speed in m/s for speed > 0.1")

# Add dummy entry (NA) at the beginning to have equal length
# i.e. each speed will be assigned to the second point of each pair
move_spatialpoint$speed <- c(NA, speeds)

#Plot log speed
spplot(move_spatialpoint, "speed", scales=list(draw=T), do.log=T)

# Reproject our moveobj to geographic coordinates
moveobj.ll <- spTransform(moveobj, CRS = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Calculate turning angles vs. speed
# Plot
plot(abs(turnAngleGc(moveobj.ll)),
     rowMeans(cbind(speed(moveobj.ll)[-length(speed(moveobj.ll))],
                    speed(moveobj.ll[-1]))))

# Clean data from outliers
move_clean$angle <- c(NA, angle(move_clean))
move_clean$speed <- c(NA, speed(move_clean))

plausibleSpeed <- move_clean$speed < 3
plausibleAngle <- !(move_clean$speed > 1 &
                      move_clean$angle < 90 |
                      move_clean$angle > 270)

# Set the first NA value to TRUE, i.e. keep the first observation
plausibleAngle[1] <- TRUE
plausibleSpeed[1] <- TRUE

while(any(!plausibleSpeed) | any(!plausibleAngle)){
  # Discard flagged observations
  move_clean <- move_clean[plausibleSpeed & plausibleAngle,]

  # Update speed and angles
  move_clean$angle <- c(NA, angle(move_clean))
  move_clean$speed <- c(NA, speed(move_clean))
  
  # Update logical vectors
  plausibleSpeed <- move_clean$speed < 2
  plausibleAngle <- !(move_clean$speed < 1 & 
                      move_clean$angle < 90 |
                      move_clean$angle > 270)
  plausibleSpeed[1] <- TRUE
  plausibleAngle[1] <- TRUE
  }


#' ### The Home Range Concept


library(adehabitatHR)

move_pts <- as(move_clean, "SpatialPoints")
move_pts.mcp <- mcp(move_pts, percent=90, unin = "m", unout = "km2")

move_pts.mcp

plotRGB(p224r63_2011m, 3,2,1, stretch="lin")
points(move_pts.mcp, col="yellow")
points(move_pts, pch=20)

move_pts.kernel <- kernelUD(move_pts, h = "href")
move_pts.kernelarea <- kernel.area(move_pts.kernel, percent = 90, unin = "km", unout="ha")

plot(move_pts, xlab=NA, ylab=NA)
plot(getverticeshr(move_pts.kernel, percent=90), lwd=2, col="grey90", add=TRUE)
points(move_pts, pch=16, cex=0.75, col = "#A020F050")

template <- raster(p224r63_2011m)
move_clean.BBMM <- brownian.bridge.dyn(move_clean, raster = template, ext = 0.1, location.error = 10, margin = 5, window.size=21, time.step=4)


#' ### Connecting movement data with remotely sensed environmental information


move_pts.fCover <- extract(fCover, move_pts)
kernelPolygon <- getverticeshr(move_pts.kernel, percent=90)
move_kern90.fCover <- extract(fCover, kernelPolygon)

boxplot(move_pts.fCover)
boxplot(move_kern90.fCover)

fCover.mask.pts <- mask(fCover, move_pts)
fCover.mask.kern90 <- mask(fCover, kernelPolygon)

library(rasterVis)
bwplot(fCover.mask.pts)
bwplot(fCover.mask.kern90)

library(dismo)

fCover.outside.kern90 <- mask(fCover, kernelPolygon, inverse=TRUE)
fCover.