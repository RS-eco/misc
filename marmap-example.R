###############################################
### Pante and Simon Bouhet - marmap package ###
###       Figures 1 and 2, PLoS ONE         ###
###############################################

### LOAD MARMAP
library(marmap)
library(lattice)

### LOAD DATA
data(nw.atlantic)
data(hawaii)
data(hawaii.sites)
as.bathy(nw.atlantic) -> atl
getNOAA.bathy(lon1=140,lon2=155,lat1=-13,lat2=0, resolution=1) -> papoue

###############################################

### Figure 1A (top left)
colorRampPalette(c("purple","lightblue","cadetblue2","cadetblue1","white")) -> blues
plot(atl, image=T, bpal=blues(100), deep=c(-6500,0), shallow=c(-50,0), step=c(500,0), 
     lwd=c(0.3,0.6), lty=c(1,1), col=c("black","black"), drawlabels=c(F,F))
scaleBathy(atl, deg=3, y=42.5,x=-74)

###############################################

### Figure 1B (middle left, 2D)
get.transect(atl,x1=-67.9, y1=41.2,x2=-55.2,y2=33.6, loc=FALSE, dis=TRUE) -> test
points(test$lon,test$lat,type="l",col=col2alpha("blue", alpha=0.5),lwd=2,lty=1)
points(min(test$lon),max(test$lat),col="blue")
points(max(test$lon),min(test$lat),col="blue")
plotProfile(test)

###############################################

### Figure 1C (middle left, 3D)
get.box(atl,x1=-67.9, y1=41.2,x2=-55.2,y2=33.6, width=2,col=2) -> out
wireframe(out, shade=T, zoom=1.1,
	aspect=c(1/4, 0.1), 
	screen = list(z = -20, x = -50),
	par.settings = list(axis.line = list(col = "transparent")),
	par.box = c(col = rgb(0,0,0,0.1)))
	
###############################################

### Figure 1D (bottom left, 3D)	
wireframe(unclass(atl), shade=T, aspect=c(1/2, 0.1),
screen = list(z = 0, x = -50),
par.settings = list(axis.line = list(col = "transparent")),
par.box = c(col = rgb(0,0,0,0.1)))

###############################################

### Figure 1E (top right)
colorRampPalette(c("red","purple","blue","cadetblue1","white")) -> blues 
plot(papoue, image=T, bpal=blues(100), #xlim=c(141,155),
	deep=c(-9000,-3000,0), shallow=c(-3000,-10, 0), step=c(1000, 1000, 0),
	col=c("lightgrey","darkgrey","black"), lwd=c(0.3,0.3,0.6), lty=c(1,1,1),
	drawlabel=c(F,F,F))

###############################################

### Figure 1F (middle right)
sites <- hawaii.sites[-c(1,4),] ; rownames(sites) <- 1:4

trans1 <- trans.mat(hawaii,min.depth=-1000)
trans2 <- trans.mat(hawaii,min.depth=-4000)

out1 <- lc.dist(trans1,sites,res="path")
out2 <- lc.dist(trans2,sites,res="path")

plot(hawaii, xlim=c(-161,-154), ylim=c(18,23), 
     deep=c(-4000,-1000,0), shallow=c(-4000,-1000,0), 
     col=c("grey80","grey40", "black"), step=c(1,1, 1), 
     lty=c(1,1,1), lwd=c(0.6,0.6,1.2), draw=c(F,F,F))
points(sites,pch=21,col="blue",bg=col2alpha("blue",.9),cex=1.2)
text(sites[,1],sites[,2],lab=rownames(sites),pos=c(3,4,1,2),col="blue")
lapply(out2,lines,col=col2alpha("blue",alpha=0.3),lwd=1,lty=1) -> dummy
lapply(out1,lines,col=col2alpha("red",alpha=0.3),lwd=1,lty=1) -> dummy

###############################################

### Figure 1G (bottom right)
plot(hawaii, lwd=0.2)

get.area(hawaii, level.inf=-4000, level.sup=-1000) -> bathyal
get.area(hawaii, level.inf=min(hawaii), level.sup=-4000) -> abyssal

image(bathyal$Lon, bathyal$Lat, bathyal$Area, col=c(rgb(0,0,0,0), rgb(.7,0,0,.3)),    add=T)
image(abyssal$Lon, abyssal$Lat, abyssal$Area, col=c(rgb(0,0,0,0), rgb(.7,.7,0.3,.3)), add=T)

round(bathyal$Square.Km, 0) -> ba
round(abyssal$Square.Km, 0) -> ab

legend(x="bottomleft",legend=c(paste("bathyal:",ba,"km2"), paste("abyssal:",ab,"km2")), 
		bty="n", fill=c(rgb(.7,0,0,.3), rgb(.7,.7,0,.3)))
		
###############################################

### Figure 2: projected Europe

bath <- getNOAA.bathy(-20,20,30,80,res=10,keep=FALSE)

r1 <- as.raster(bath)
projection <- "+proj=ortho"
r2 <- raster::projectRaster(r1,crs=projection)
as.bathy(r2) -> bath2

pal <- colorRampPalette(c("#4F503E","#3073BB","#4698CD","#C6D7D5","#93856B"))

tiff(filename="Europe_projected.tiff", 
     width=17,height=17, unit="cm", res=650, pointsize=12, 
     compression="lzw", bg="white", type=c("cairo","Xlib","quartz"))
plot.bathy(bath2,lwd=.03,n=1000,land=T,col=pal(1000),xlab="",ylab="",axes=F)
plot.bathy(bath2,lwd=.5,deep=0,shallow=0,step=0,col=pal(1), add=T)
dev.off()
