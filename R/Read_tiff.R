library(raster)
library(rgdal)
infile <-"C:Data/DEM100/demoriginal.tif"
dem<-raster (infile)
ncol(dem)
nrow(dem)
swdem<-getValueBlock(dem,nrow=,nrows=,ncol= ,ncols=)



# IMport Wind ncdf data ???
infile<-"C:Data/Wind/X144-1960-70-uwindx4.nc"
windu<-raster(infile)
infile<-"C:Data/Wind/X144-1960-70-uwindx4.nc"
windv<-raster(infile)
