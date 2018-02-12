# Convert county shape files to raster in new projection
root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/"

dir_dem<-paste(in.root,"DEM/",sep="")
dir_counties<-paste(root,"counties/",sep="")
print(dir_counties)
dem<-raster(paste(dir_dem,"dem.tif",sep=""))
plot(dem)

# Define projections
latlong = "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"
new.crs<-"+init=epsg:3857"

# Load county shape file and crop to SW area and reproject 
hcounties.shp<- readOGR(dsn = path.expand(paste(root,"OSdata/bdline_essh_gb/Data/Supplementary_Historical",sep="")), layer = "Boundary-line-historic-counties_region")
swcounties.shp<-(crop(hcounties.shp,dem))
swcounties.shp<-spTransform(swcounties.shp, new.crs)

# First create raster template in new projection for whole SW
e<-extent(dem)
template.r<-projectExtent(dem,crs=new.crs)
res(template.r)<-c(100,100)
print(template.r)

# Create and write county raster masks in new projection 
sel <- which(swcounties.shp$Name == "Cornwall")
county<-swcounties.shp[sel,]
e<-extent(county)
countymask.r<-crop(template.r,e)
r<-raster(nrows=nrow(countymask.r), ncols=ncol(countymask.r), ext=extent(countymask.r),
          crs=crs(countymask.r), resolution=res(countymask.r), vals=NULL)
r<- rasterize(county, r)
plot(r)
writeRaster(r,file=paste(dir_counties,"cornwallmask.tif",sep=""),overwrite=TRUE)

sel <- which(swcounties.shp$Name == "Devon")
county<-swcounties.shp[sel,]
e<-extent(county)
countymask.r<-crop(template.r,e)
r<-raster(nrows=nrow(countymask.r), ncols=ncol(countymask.r), ext=extent(countymask.r),
          crs=crs(countymask.r), resolution=res(countymask.r), vals=NULL)
r<- rasterize(county, r)
plot(r)
writeRaster(r,file=paste(dir_counties,"devonmask.tif",sep=""),overwrite=TRUE)

sel <- which(swcounties.shp$Name == "Dorset")
county<-swcounties.shp[sel,]
e<-extent(county)
countymask.r<-crop(template.r,e)
r<-raster(nrows=nrow(countymask.r), ncols=ncol(countymask.r), ext=extent(countymask.r),
          crs=crs(countymask.r), resolution=res(countymask.r), vals=NULL)
r<- rasterize(county, r)
plot(r)
writeRaster(r,file=paste(dir_counties,"dorsetmask.tif",sep=""),overwrite=TRUE)

sel <- which(swcounties.shp$Name == "Somerset")
county<-swcounties.shp[sel,]
e<-extent(county)
countymask.r<-crop(template.r,e)
r<-raster(nrows=nrow(countymask.r), ncols=ncol(countymask.r), ext=extent(countymask.r),
          crs=crs(countymask.r), resolution=res(countymask.r), vals=NULL)
r<- rasterize(county, r)
plot(r)
writeRaster(r,file=paste(dir_counties,"somersetmask.tif",sep=""),overwrite=TRUE)

