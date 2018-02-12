# Write topographical rasters for each county in projection for app/leaflet
# SLope, elevation and aspect by county
#######################################
# Libraries and Directories
library(rgdal)
library(raster)
library(ggplot2)

root<-"/home/ISAD/jm622/Data2015/"   # Source data and output data
in.root<-"/data/jm622/ModelData/" # model input data 
root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/"

# Output dirs
dir_shinydata<-paste(root,"shinydata/",sep="")
dir_rasters<-paste(root,"shinydata/rasters/",sep="")

# Input dirs
dir_dem<-paste(in.root,"DEM/",sep="")
dir_terrain<-paste(in.root,"Terrain/",sep="")
#dir_results<-paste(root,"Outputs/",sep="") # dir of results rasters for whole of SW
dir_results<-paste(root,"proxyt100/riskmaps/",sep="")
dir_counties<-paste(dir_shinydata,"countymasks/",sep="")

# MAC TESTING ONLY
#root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/"
# dir_counties<-paste(root,"counties/",sep="")
#######################################

# Define projections
latlong = "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"
new.crs<-"+init=epsg:3857"
#new.crs<-'+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'

# Read county mask files
cornwall.r<-raster(paste(dir_counties,"cornwallmask.tif",sep=""))
devon.r<-raster(paste(dir_counties,"devonmask.tif",sep=""))
dorset.r<-raster(paste(dir_counties,"dorsetmask.tif",sep=""))
somerset.r<-raster(paste(dir_counties,"somersetmask.tif",sep=""))

# County list and labels
crs(cornwall.r)<-new.crs
crs(devon.r)<-new.crs
crs(dorset.r)<-new.crs
crs(somerset.r)<-new.crs
county.txt<-c("cornwall","devon","dorset","somerset")
county<-list(cornwall.r,devon.r,dorset.r,somerset.r) 

# Link terrain map data to files
dem<-raster(paste(dir_dem,"dem.tif",sep=""))
slope<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem)
aspect<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem)
crs(dem)<-ukgrid 
crs(slope)<-ukgrid 
crs(aspect)<-ukgrid

# 1 Reproject terrain data then mask to extract county rasters
for (c in 1:length(county)){
  template.r<-county[[c]]
  crs(template.r)<-new.crs
  print (template.r)
  e<-extent(template.r)
  dem.out<-projectRaster(dem,to=template.r)
  dem.out<-mask(crop(dem.out,e), template.r)
  slope.out<-projectRaster(slope,to=template.r)
  slope.out<-mask(crop(slope.out,e), template.r)
  aspect.out<-projectRaster(aspect,to=template.r)
  aspect.out<-mask(crop(aspect.out,e), template.r)
  
  writeRaster(dem.out,file=paste(dir_rasters,"elevation_",county.txt[c],".tif",sep=""),overwrite=TRUE)
  writeRaster(slope.out,file=paste(dir_rasters,"slope_",county.txt[c],".tif",sep=""),overwrite=TRUE)
  writeRaster(aspect.out,file=paste(dir_rasters,"aspect_",county.txt[c],".tif",sep=""),overwrite=TRUE)
  # Plot and write histograms
  #dem.h<-plothist(getValues(dem.out))
  #slope.h<-plothist(getValues(slope.out))
  #aspect.h<-plothist(getValues(aspect.out))
}
