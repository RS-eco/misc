# title         : worldmaps.R
# purpose       : download, import, resampling and export of publically availaible global datasets (gridded maps) + some examples of how these maps could be used for spatial modelling (species distribution modeling);
# reference     : [http://spatial-analyst.net/wiki/index.php?title=Global_datasets]
# producer      : Prepared by T. Hengl
# last update   : In Amsterdam, NL, 8 May 2010.
# inputs        : repository of world maps, various publicly available point datasets;
# outputs       : resampled grids and visualizations of thematic maps;
# remarks 1     : in order to download the maps, you need to be on-line of course;

library(maptools)
library(rgdal)
library(RSAGA)
library(gstat)
library(spatstat)
library(RCurl)
library(XML)
library(R2HTML)
library(colorspace)

# ------------------------------------------------------------
# Download of maps and export to SAGA GIS format
# ------------------------------------------------------------

# location of maps:
URL <- "http://spatial-analyst.net/worldmaps/"
# list all maps available:
items <- strsplit(getURL(URL), "\n")[[1]]
# convert to a character vector:
zip.list <- items[grep(items, pattern=".zip")]
zip.list <- unlist(lapply(strsplit(zip.list, '\"> '), function(x){x[length(x)]}))
zip.list <- unlist(lapply(strsplit(zip.list, '</a></li>'), function(x){x[length(x)]}))
zip.list

# download the zipped maps one by one:
for(i in 1:length(zip.list)) {
  download.file(paste(URL, zip.list[i], sep=""), destfile=paste(getwd(), "/", zip.list[i], sep=""))
  unzip(zipfile=zip.list[i], exdir=getwd())
#  unlink(zip.list[i])
}

# Get info about a map:
GDALinfo(set.file.extension(zip.list[1], ".tif"))
worldmaps <- readGDAL(set.file.extension(zip.list[1], ".tif")) 
names(worldmaps) <- sub(".zip", "", zip.list[1])
proj4string(worldmaps) <- CRS("+proj=longlat +datum=WGS84")

# ------------------------------------------------------------
# Working with meta-data
# ------------------------------------------------------------

# To get the real metadata use e.g.:
meta.list <- list(rep(NA, length(zip.list)))
for(j in 1:length(zip.list)){
tmp <- t(matrix((unlist(strsplit(readLines(file(set.file.extension(zip.list[j], ".rdc"))), split=": "))), ncol=2, byrow=T))
tmp.meta <- unclass(tmp[2,])
names(tmp.meta) <- as.character(tmp[1,])
meta.list[[j]] <- tmp.meta
close(file(set.file.extension(zip.list[j], ".rdc")))
}
# print all metadata:
print(attr(meta.list[[j]], "")[1:35])
# proj4string:
print(meta.list[[j]][12])

# write to table:
meta.table <- data.frame(title=rep(NA, length(zip.list)), description=rep(NA, length(zip.list)), source=rep(NA, length(zip.list)), URL=rep(NA, length(zip.list)), scale=rep(NA, length(zip.list)), ref=rep(NA, length(zip.list)), owner=rep(NA, length(zip.list)))
for(j in 1:length(zip.list)){
   meta.table[j,"title"] <- meta.list[[j]][2]
   meta.table[j,"description"] <- meta.list[[j]][11]
   meta.table[j,"source"] <- meta.list[[j]][30]
   meta.table[j,"URL"] <- meta.list[[j]][31]
   meta.table[j,"scale"] <- meta.list[[j]][32]
   meta.table[j,"ref"] <- meta.list[[j]][33]
   meta.table[j,"owner"] <- meta.list[[j]][34]
}
write.table(meta.table, file="meta_worldmaps.txt", sep=";", quote=FALSE, row.names=FALSE)
# write metadata to a HTML file:
HTML(as.data.frame(meta.list[[1]]), file=set.file.extension(zip.list[1], ".html")) 

# Write categorical maps to SAGA format (with color lookup table)
# define map name:
mapname <- "glc2000"
writeGDAL(worldmaps[mapname], paste(mapname, ".sdat", sep=""), "SAGA")
# read the metadata from the *.rdc file:
tmp <- t(matrix((unlist(strsplit(readLines(file(set.file.extension(mapname, ".rdc"))), split=" : "))), ncol=2, byrow=T))
meta <- unclass(tmp[2,])
# read the color legend:
tmp <- read.table(set.file.extension(mapname, ".PAL"), col.names=c("Class", "R", "G", "B"))
# convert to BGR codes:
BGR <- (tmp$B * 65536) + (tmp$G * 256) + tmp$R
# write a lookup table for SAGA GIS:
filename <- file(set.file.extension(mapname, ".txt"), "w", blocking=FALSE)
write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
write(paste(BGR[1], meta[36+i], paste("CL", tmp$Class[1], sep=""), (tmp$Class[1])-0.9, (tmp$Class[1])+0.1, sep="\t"), filename, append=TRUE)
for(i in 2:length(BGR)){
write(paste(BGR[i], meta[36+i], paste("CL", tmp$Class[i], sep=""), (tmp$Class[i-1])+0.1, (tmp$Class[i])+0.1, sep="\t"), filename, append=TRUE)}
close(filename)

# generate a new .rdc file:
mapname <- "PRECm"
# read the BLANK metadata sheet:
tmp <- t(matrix((unlist(strsplit(readLines(file(paste(URL, "BLANK.rdc", sep=""))), split=" : "))), ncol=2, byrow=T))
map.info <- GDALinfo(paste(mapname, ".tif", sep=""))
driver.name <- paste(gdalDrivers()[gdalDrivers()$name==attr(map.info, "driver"), "long_name"])
# write the GDLAinfo metadata into the sheet:
filename <- file(set.file.extension(mapname, ".rdc"), "w", blocking=FALSE)
for(i in 1:ncol(tmp)){
if(i==1){ write(paste(tmp[1,1], " : ", driver.name, " file format", sep=""), filename) }
if(i==2){ write(paste(tmp[1,2], " : ", strsplit(attr(map.info, "file"), ".tif")[[1]][1], sep=""), filename, append=TRUE) }
if(i==3){ write(paste(tmp[1,3], " : ", format(Sys.Date(), "%d-%b-%Y"), sep=""), filename, append=TRUE) }
if(i>3&i<6){ write(paste(tmp[,i], collapse=" : "), filename, append=TRUE) }
if(i==6){ write(paste(tmp[1,6], " : ", paste(attr(map.info, "df")$GDType), sep=""), filename, append=TRUE) }
if(i==7){ write(paste(tmp[1,7], " : ", paste(attr(map.info, "driver")), sep=""), filename, append=TRUE) }
if(i==8){ write(paste(tmp[1,8], " : ", paste(map.info[2]), sep=""), filename, append=TRUE) }
if(i==9){ write(paste(tmp[1,9], " : ", paste(map.info[1]), sep=""), filename, append=TRUE) }
if(i>9&i<12){ write(paste(tmp[,i], collapse=" : "), filename, append=TRUE) }
if(i==12){ write(paste(tmp[1,12], " : ", paste(attr(map.info, "projection")), sep=""), filename, append=TRUE) }
if(i>12&i<16){ write(paste(tmp[,i], collapse=" : "), filename, append=TRUE) }
if(i==16){ write(paste(tmp[1,16], " : ", paste(map.info[4]), sep=""), filename, append=TRUE) }
if(i==17){ write(paste(tmp[1,17], " : ", paste(map.info[4]+map.info[2]*map.info[6]), sep=""), filename, append=TRUE) }
if(i==18){ write(paste(tmp[1,18], " : ", paste(map.info[5]), sep=""), filename, append=TRUE) }
if(i==19){ write(paste(tmp[1,19], " : ", paste(map.info[5]+map.info[1]*map.info[6]), sep=""), filename, append=TRUE) }
if(i>19&i<21){ write(paste(tmp[,i], collapse=" : "), filename, append=TRUE) }
if(i==21){ write(paste(tmp[1,21], " : ", paste(map.info[6]), sep=""), filename, append=TRUE) }
if(i==22){ write(paste(tmp[1,22], " : ", paste(attr(map.info, "df")$Bmin), sep=""), filename, append=TRUE) }
if(i==23){ write(paste(tmp[1,23], " : ", paste(attr(map.info, "df")$Bmax), sep=""), filename, append=TRUE) }
if(i>23&i<28){ write(paste(tmp[,i], collapse=" : "), filename, append=TRUE) }
if(i==28){ write(paste(tmp[1,28], " : ", paste(attr(map.info, "df")$Bmax), sep=""), filename, append=TRUE) }
if(i>28){ write(paste(tmp[,i], collapse=" : "), filename, append=TRUE) }
}
close(filename)
# other meta-data you will need to edit manually!

# ------------------------------------------------------------
# Creating ground overlays
# ------------------------------------------------------------

# Example: make a ground overlay
worldmaps.kml <- GE_SpatialGrid(worldmaps, maxPixels=7200)
png(file=set.file.extension(zip.list[1], ".png"), width=worldmaps.kml$width, height=worldmaps.kml$height, bg="transparent")
par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
image(as.image.SpatialGridDataFrame(worldmaps[1]), col=rainbow(length(levels(as.factor(worldmaps@data[,1])))), xlim=worldmaps.kml$xlim, ylim=worldmaps.kml$ylim) # this can take time!!!
kmlOverlay(worldmaps.kml, kmlfile=set.file.extension(zip.list[1], ".kml"), imagefile=set.file.extension(zip.list[1], ".png"), name=set.file.extension(zip.list[1], ".tif"))
dev.off()

# ------------------------------------------------------------
# Species distribution modelling
# Spreading of the Sturnella magna over the USA
# ------------------------------------------------------------

# create an output directory for maps:
Sys.chmod(getwd(), mode="7777") # allow read/write
dir.create(path="out"); 
out.env <- rsaga.env(workspace=paste(getwd(), "/out", sep=""))
grids.env <- rsaga.env(workspace=paste(getwd(), "/grids", sep=""))
# USA Coordinate system:
AEA <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Get the occurrence records for Sturnella magna from AKN e.g. [http://www.avianknowledge.net]; GBIF [http://data.gbif.org/species/13815891/] unfortunately limits the results to 250,000 records!
# Load the AKN records to R (point data):
load(url("http://spatial-analyst.net/book/system/files/AKN_Sturnella_magna_USA.RData"))
                      
# download and import the USA borders map:
download.file("http://www.census.gov/geo/cob/bdy/st/st00shp/st99_d00_shp.zip", destfile=paste(getwd(), "/usa_borders.zip", sep=""))
load(url("http://spatial-analyst.net/book/system/files/usa_borders.RData"))
usa_borders <- spTransform(usa_borders.ll, CRS(AEA))

# plot density counts;
AKN.pntXY.2009 <- subset(AKN.pntXY, AKN.pntXY$Observation.Date>as.Date("2008-12-31")&AKN.pntXY$Observation.Date<as.Date("2009-12-31"))
bubble.plt <- bubble(AKN.pntXY.2009, "Dens", scales=list(draw=T, cex=.7), maxsize=2, pch=19, lwd=1.5, col="blue", sp.layout=list("sp.lines", usa_borders), main="Sturnella magna (AKN) in 2009")
print(bubble.plt)

# focus on some categorical maps:
map.list.f <- c("ecoflor.zip", "glwd31.zip", "globcov.zip")
# and a list of continous maps:
map.list <- c("cforest.zip", "PCEVI1.zip", "PCEVI2.zip", "PCEVI3.zip", "PCEVI4.zip", "dcoast.zip", "pcnligh1.zip", "rcrops.zip", "slope.zip", "globedem.zip", "LSTNs.zip", "PRECm.zip")

# create temp directory to put resampled maps:
dir.create(path="grids")

# resample worldmaps to local coordinate system (Continental USA):
for(j in 1:length(map.list)) {
# resample to 5 km local grid:
 unlink("grids/tmp.sdat")
 system(paste("C:\\PROGRA~2\\FWTOOL~1.7\\bin\\gdalwarp ", set.file.extension(map.list[j], ".tif"), " -t_srs \"", AEA, "\" grids\\tmp.sdat -of \"SAGA\" -r bilinear -te -2405000 265000 2295000 3225000 -tr 5000 5000", sep=""))
 rsaga.geoprocessor(lib="io_grid", module=0, param=list(GRID="tmp.sgrd", FILE=paste("m_", set.file.extension(map.list[j], ".asc"), sep=""), FORMAT=1, PREC=3), env=grids.env, show.output.on.console=FALSE) 
}

# categorical maps:
for(j in 1:length(map.list.f)) {
# resample to 5 km local grid:
 unlink("grids/tmp.sdat")
 system(paste("C:\\PROGRA~2\\FWTOOL~1.7\\bin\\gdalwarp ", set.file.extension(map.list.f[j], ".tif"), " -t_srs \"", AEA, "\" grids\\tmp.sdat -of \"SAGA\" -r near -te -2405000 265000 2295000 3225000 -tr 5000 5000", sep=""))
 rsaga.geoprocessor(lib="io_grid", module=0, param=list(GRID="tmp.sgrd", FILE=paste("m_", set.file.extension(map.list.f[j], ".asc"), sep=""), FORMAT=1, PREC=0), env=grids.env, show.output.on.console=FALSE) 
}

# read the resampled worldmaps to R:
gridmaps <- readGDAL(paste("grids/m_", set.file.extension(map.list[1], ".asc"), sep=""))
for(j in 2:length(map.list)) {
  gridmaps@data[j] <- round(readGDAL(paste("grids/m_", set.file.extension(map.list[j], ".asc"), sep=""),  silent=TRUE)$band1, 3) 
}
# categorical maps:
for(j in 1:length(map.list.f)) {
   gridmaps@data[(length(map.list)+j)] <- as.factor(readGDAL(paste("grids/m_", set.file.extension(map.list.f[j], ".asc"), sep=""),  silent=TRUE)$band1)
}
names(gridmaps) <- sub(".zip", "", c(map.list, map.list.f))
str(gridmaps@data)
proj4string(gridmaps) <- CRS(AEA)
cellsize <- gridmaps@grid@cellsize[1]
spplot(gridmaps[1], col.regions=rainbow(20), sp.layout=list("sp.lines", usa_borders), Main=map.list[1])

# fix some ASCII layers (problems with writing NA values):
fix.layers <- c("PCEVI1", "PCEVI2", "PCEVI3", "PCEVI4", "LSTNs")
for(j in 1:length(fix.layers)) { write.asciigrid(gridmaps[fix.layers[j]], fname=paste("grids/m_", fix.layers[j], ".asc", sep="")) }
unlink("grids/m_dcoast.asc")

# ------------------------------------------------------------
# Species density distribution 
# (kernel smoothing)
# ------------------------------------------------------------

count_13815891.ppp <- ppp(coordinates(AKN.pntXY.2009)[,1], coordinates(AKN.pntXY.2009)[,2], marks=AKN.pntXY.2009$Dens, window=as(gridmaps["globcov"], "owin"))  # [http://data.gbif.org/species/13815891/]
count_13815891.dens <- density(count_13815891.ppp, sigma=2*cellsize)
gridmaps$dens <- as(count_13815891.dens, "SpatialGridDataFrame")$v
gridmaps$densr <- gridmaps$dens/(max(gridmaps$dens, na.rm=T))
plt1 <- spplot(gridmaps["densr"], col.regions=grey(rev((1:60)^3/60^3)), sp.layout=list("sp.lines", usa_borders), main="Kernel density", at=seq(0,1,1/60)) 

# ------------------------------------------------------------
# Habitat Suitability analysis using MaxEnt:
# MaxEnt software v3.1 [http://www.cs.princeton.edu/~schapire/maxent/]
# ------------------------------------------------------------

# position of MaxEnt and local directories:
MaxEnt <- "C:\\MaxEnt\\maxent.jar"
MaxEnt.layers <- paste(gsub("/", "\\\\", getwd()), "\\grids", sep="")
MaxEnt.out <- paste(gsub("/", "\\\\", getwd()), "\\MEout", sep="")
MaxEnt.samples <- paste(gsub("/", "\\\\", getwd()), "\\MEsamples", sep="")
dir.create(path="MEout"); dir.create(path="MEsamples"); dir.create(path="out")

# prepare the presence only maps:
writeOGR(AKN.pntXY.2009["Dens"], "out/pr_13815891.shp", "Dens", "ESRI Shapefile")
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(GRID="pr_13815891.sgrd", INPUT="pr_13815891.shp", FIELD=1, LINE_TYPE=1, USER_CELL_SIZE=cellsize, USER_X_EXTENT_MIN=gridmaps@bbox[1,1]+cellsize/2, USER_X_EXTENT_MAX=gridmaps@bbox[1,2]-cellsize/2, USER_Y_EXTENT_MIN=gridmaps@bbox[2,1]+cellsize/2, USER_Y_EXTENT_MAX=gridmaps@bbox[2,2]-cellsize/2), env=out.env)
# read to R:
gridmaps$pr_13815891 <- readGDAL("out/pr_13815891.sdat")$band1
dens.weight <- as.im(as.image.SpatialGridDataFrame(gridmaps["pr_13815891"]))
pr.13815891 <- rpoint(as.integer(sum(AKN.pntXY.2009$Dens, na.rm=TRUE)), f=dens.weight) # how to determine how many points to simulate?
plot(pr.13815891)

# export the presence only values and run MaxEnt:
species.csv <- data.frame(x=pr.13815891$x, y=pr.13815891$y, sp=rep("Sturnella_magna", length(pr.13815891$x)))
write.csv(species.csv[,c("sp","x","y")], "MEsamples/Sturnella_magna.csv", quote=FALSE, row.names=FALSE)
# run MaxEnt:
system(command=paste("java -mx1000m -jar ", MaxEnt, " environmentallayers=", MaxEnt.layers, " samplesfile=", MaxEnt.samples, "\\Sturnella_magna.csv", " togglelayertype=m_", sub(".zip", "", map.list.f[1]), " togglelayertype=m_", sub(".zip", "", map.list.f[2]), " togglelayertype=m_", sub(".zip", "", map.list.f[3]), " outputdirectory=", MaxEnt.out, " randomtestpoints=25 maximumiterations=100 redoifexists autorun nowarnings notooltips", sep=""))

# print the results:
MEout.list <- htmlTreeParse("MEout/Sturnella_magna.html")
MEout.list[[3]][[1]][[2]][[65]]; MEout.list[[3]][[1]][[2]][[52]]; MEout.list[[3]][[1]][[2]][[56]]; MEout.list[[3]][[1]][[2]][[67]]

# import the HSI map to R:
gridmaps$HSI <- readGDAL("MEout/Sturnella_magna.asc")$band1 
plt2 <- spplot(gridmaps["HSI"], col.regions=bpy.colors(30), main="Habitat index (0-100)", sp.layout=list("sp.lines", usa_borders))

print(plt1, split=c(1,1,1,2), more=T)
print(plt2, split=c(1,2,1,2), more=F)

# end of script