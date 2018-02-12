library(shiny)
library(raster)
library(leaflet)
library(magrittr)
dem<-raster(paste(dir_dem,"dem.tif",sep=""))
slope<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem)
aspect<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem)

e<-extent(120000,240000,10000,140000)
dem<-crop(dem,e)
slope<-crop(slope,e)
aspect<-crop(aspect,e)

epsg<-"+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
#epsg<-"+init=epsg:3857" # mercator
dem.epsg<-projectRaster(dem,crs=epsg,res=c(100,100),filename=paste(dir_dem,"demepsg.tif"),overwrite=TRUE)
slope.epsg<-projectRaster(slope,crs=epsg,res=c(100,100),filename=paste(dir_dem,"slopeepsg.tif"),overwrite=TRUE)
aspect.epsg<-projectRaster(aspect,crs=epsg,res=c(100,100),filename=paste(dir_dem,"aspectepsg.tif"),overwrite=TRUE)

plot(dem.epsg)

leaflet() %>% addTiles() %>%
  addRasterImage(dem.epsg,  opacity = 0.7)

pal <-colorQuantile(c("#B2FF66","#66CC00","#4C9900","#336600","#193300"),
                    NULL, n = 5, na.color="#FFFFFF")
leaflet() %>% addTiles() %>%
  addRasterImage(dem.epsg,  opacity = 0.7) %>%
  addLegend(pal = pal, values = values(dem.epsg),
            title = "Legend")

##############
# Using mapview package
library(devtools)
library(mapview)

mapviewOptions()
  
map.s<-stack(dem.epsg,slope.epsg,aspect.epsg)
mapview(map.s,maxpixels=4300000,alpha.regions=0.75,legend=F) # with underlay 

plainview(dem)
plainview(stack(dem,slope,aspect))
