# Convert county shape files to OS rasters
dir_counties<-paste(root,"counties/",sep="")
print(dir_counties)

hcounties.shp<- readOGR(dsn = path.expand(paste(root,"OSdata/bdline_essh_gb/Data/Supplementary_Historical",sep="")), layer = "Boundary-line-historic-counties_region")

# extract county polygons - use historic counties data 
sel <- which(hcounties.shp$Name == "Cornwall")
county<-hcounties.shp[sel,]
e<-extent(county)
template.r<-crop(dem,e)
r<-raster(nrows=nrow(template.r), ncols=ncol(template.r), ext=extent(template.r),
          crs=crs(template.r), resolution=res(template.r), vals=NULL)
cornwall.r<- rasterize(county, r)
plot(cornwall.r)
writeRaster(cornwall.r,file=paste(dir_counties,"cornwall.tif",sep=""),overwrite=TRUE)


sel <- which(hcounties.shp$Name == "Devon")
county<-hcounties.shp[sel,]
e<-extent(county)
template.r<-crop(dem,e)
r<-raster(nrows=nrow(template.r), ncols=ncol(template.r), ext=extent(template.r),
          crs=crs(template.r), resolution=res(template.r), vals=NULL)
devon.r<- rasterize(county, r)
plot(devon.r)
writeRaster(devon.r,file=paste(dir_counties,"devon.tif",sep=""),overwrite=TRUE)

sel <- which(hcounties.shp$Name == "Dorset")
county<-hcounties.shp[sel,]
e<-extent(county)
template.r<-crop(dem,e)
r<-raster(nrows=nrow(template.r), ncols=ncol(template.r), ext=extent(template.r),
          crs=crs(template.r), resolution=res(template.r), vals=NULL)
dorset.r<- rasterize(county, r)
plot(dorset.r)
writeRaster(dorset.r,file=paste(dir_counties,"dorset.tif",sep=""),overwrite=TRUE)

sel <- which(hcounties.shp$Name == "Somerset")
county<-hcounties.shp[sel,]
e<-extent(county)
template.r<-crop(dem,e)
r<-raster(nrows=nrow(template.r), ncols=ncol(template.r), ext=extent(template.r),
          crs=crs(template.r), resolution=res(template.r), vals=NULL)
somerset.r<- rasterize(county, r)
plot(somerset.r)
writeRaster(somerset.r,file=paste(dir_counties,"somerset.tif",sep=""),overwrite=TRUE)


counties<-c("cornwall","devon","dorset","somerset")
for (n in 1:length(counties)){
  print(counties[n])
  county.r<-raster(paste(dir_counties,counties[n],".tif",sep=""))
  plot(county.r)
}