# Libraries
library(R.utils) # IMPORTANT - hides raster functions requiring raster::
library(ncdf4)
library(raster)
library(rgdal)


# Directories
# Templates - 5km cells
dir_grids<-"C:/Data2015/Templates/"


# Temperature
#dir_temp<-"C:/Data2015/Temp5km/"
#dir_hrtemp<-"C:/Data2015/Temp5km/hourly/"

# Coast
dir_percland<-"C:/Data2015/CoastEffect/percland/"

# SST
dir_sst<-"C:/Data2015/SST/"  
#gunzip(paste(dir_sst,"HadISST_sst.nc.gz",sep=""))
dir_sstm<-paste(dir_sst,"monthly/",sep="")
dir_ssth<-paste(dir_sst,"hourly/",sep="")
dir_upwindsea<-"C:/Data2015/CoastEffect/upwindmaps/"

#Wind
dir_wind<-"C:/Data2015/Wind/"
dir_wind<-"C:/Data2015/Wind/"
dir_shelter<-"C:/Data2015/Wind/Shelter/"

#dir_temp<-"C:/Data2015/Temp5km/"
#dir_hrtemp<-"C:/Data2015/Temp5km/hourly/"

#dir_zip<-"~/Documents/Exeter/Data2015/Temp5km/zip/"
dir_temp<-"~/Documents/Exeter/Data2015/Temp5km/extract/"
dir_hrtemp<-"~/Documents/Exeter/Data2015/Temp5km/hourly/" # dir for output files

dir_sst<-"~/Documents/Exeter/Data2015/sst/"
#dir_sst<-"C:/Data2015/SST/"  
dir_sstm<-paste(dir_sst,"monthly/",sep="")
dir_ssth<-paste(dir_sst,"hourly/",sep="")

dir_upwindsea<-"C:/Data2015/CoastEffect/upwindmaps/"



# CONSTANTS

# DEM
latlong <- "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"

demuk<-raster("C:/Data2015/DEM100/demoriginal.tif", crs=(ukgrid))

# Define area of interest (no buffer)
#e.dem<-extent(c(70000,420000,0,180000 )) # includes scilly isles
e.dem<-extent(c( 120000,420000,10000,180000 )) # excludes scilly isles

buffer<-30000
e.buf<-extent(xmin(dem)-buffer,xmax(dem)+buffer,ymin(dem)-buffer,ymax(dem)+buffer)