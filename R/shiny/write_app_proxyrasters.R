# Prepare files required for shinyapp
# Reproject  and save in appropriate directories 
# REQUIRES: reprojected raster masks already created - in write_countymask_rasters
# REQUIRES: results rasters for each county for each year (in OS projection)

# Job variables
args <-commandArgs(trailingOnly = TRUE)
print(args)
county.txt <- as.integer(args[1])
print(paste("County=",county))

#######################################
# Libraries and Directories
library(rgdal)
library(raster)
library(ggplot2)

root<-"/home/ISAD/jm622/Data2015/"   # Source data and output data
in.root<-"/data/jm622/ModelData/" # model input data 

# Output dirs
dir_shinydata<-paste(root,"shinydata/",sep="")
dir_rasters<-paste(root,"shinydata/rasters/",sep="")

# Input dirs
dir_dem<-paste(in.root,"DEM/",sep="")
dir_terrain<-paste(in.root,"Terrain/",sep="")
#dir_results<-paste(root,"Outputs/",sep="") # dir of results rasters for whole of SW
dir_results<-paste(root,"proxyt100/riskmaps/",sep="")
dir_counties<-paste(dir_shinydata,"countymasks/",sep="")

#  TESTING ONLY
#root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/"
# dir_counties<-paste(root,"counties/",sep="")
#county.txt<-c("cornwall","devon","dorset","somerset")
#county<-list(cornwall.r,devon.r,dorset.r,somerset.r) 
#######################################

# Define result years to be processed
years<-c(1983:2013)
print(years)

# Results variables 
statnames<-c("gdd10_gs","gdd5_gs","tmean_gs","tmin_year","tmax_year",
             "t20_gsdays","t25_gsdays","t30_gsdays", 
             "lastspfr_doy", "firstautfr_doy", "frostfree_days",
             "fl_tmean", "fl_numday")
labs<-c("Degree days","Degree days","degrees C","degrees C","degrees C","Number of days","Number of days","Number of days",
        "Day of Year (1:366)","Day of Year (1:366)","Number of days",
        "degrees C", "Number of days")

# Define projections
latlong = "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"
new.crs<-"+init=epsg:3857"
#new.crs<-'+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'

# Read county mask files
countymask.r<-raster(paste(dir_counties,county,"mask.tif",sep=""))
crs(countymask.r)<-new.crs

# Link terrain map data to files
dem<-raster(paste(dir_dem,"dem.tif",sep=""))
slope<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem)
aspect<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem)
crs(dem)<-ukgrid 
crs(slope)<-ukgrid 
crs(aspect)<-ukgrid

######################################
# Extract county rastersfor each statistic
######################################

# Text below used in creation of rasters masks
# First set raster template in new projection for whole results area based on dem
  #e<-extent(dem)
  #template.r<-projectExtent(dem,crs=new.crs)
  #res(template.r)<-c(100,100)
  #print(template.r)

# 2 Crop and reproject results rasters

for (y in 1:length(years)){
  for (n in 1:length(statnames)){
    print(years[y])
    print(statnames[n])
    print (county.txt)
    results.file<-paste(dir_results,statnames[n],"_",years[y],"_",county.txt,".tif",sep="")
    print(results.file)
    r<-raster(results.file) # raster for whole area
    crs(r)<-ukgrid # old crs (OS)
    ### Reproject & write county rasters
    e<-extent(countymask.r)
    county.r<-projectRaster(r,to=countymask.r) # reproject 
    county.r<-mask(crop(county.r,e), countymask.r)
    plot(county.r,main=statnames[n])
    # Write Raster -DECIDE OUTPUT DIR
    r.filename<-paste(dir_rasters,statnames[n],"_",years[y],"_",county.txt,".tif",sep="")
    print(r.filename)
    writeRaster(county.r,file=r.filename,overwrite=TRUE)
  } # for n
} # for y


# 3. Write timeseries for each variable
dir_timeseries<-paste(root,"shinydata/timeseries/",sep="")

for (n in 1:length(statnames)){
  numcells<-ncell(raster(paste(dir_rasters,"elevation_",county.txt,".tif",sep="")))
  var.m<-matrix(data=NA,ncol=length(years),nrow=numcells)
  for (y in 1:length(years)){
    r.filename<-paste(dir_rasters,statnames[n],"_",years[y],"_",county.txt,".tif",sep="")
    r<-raster(r.filename)
    values<-getValues(r)
    var.m[,y]<-values
  } # year
  var.filename<-paste(dir_timeseries,statnames[n],"_",county.txt,"_timeseries.R",sep="")
  print(var.filename)
  save(var.m,file=var.filename)
} # stat


# 4. Write summary rasters for each variable

for (n in 1:length(statnames)){
  print(years[y])
  print(statnames[n])
  #plot(r)
  print (county.txt)
  # Define results stack from which summary stats calculated
  allyears.s<-stack()
  for (y in 1:length(years)){
    print(years[y])
    r.filename<-paste(dir_rasters,statnames[n],"_",years[y],"_",county.txt,".tif",sep="")
    print(r.filename)
    r<-raster(r.filename) # raster for whole area
    allyears.s<-stack(allyears.s,r)
  } # for y
  ### Calculate and write summary rasters
  first.yr<-min(years)
  last.yr<-max(years)
  mean.r<-calc(allyears.s,mean); plot(mean.r,main=paste("Mean ",statnames[n]))
  max.r<-calc(allyears.s,max); plot(max.r,main=paste("Max ",statnames[n]))
  min.r<-calc(allyears.s,min); plot(min.r,main=paste("Min ",statnames[n]))
  # Write Rasters -DECIDE OUTPUT DIR
  mean.file<-paste(dir_rasters,statnames[n],"_mean_",first.yr,"-",last.yr,"_",county.txt,".tif",sep="")
  print(mean.file)
  writeRaster(mean.r,file=mean.file,overwrite=TRUE)
  max.file<-paste(dir_rasters,statnames[n],"_max_",first.yr,"-",last.yr,"_",county.txt,".tif",sep="")
  print(max.file)
  writeRaster(max.r,file=max.file,overwrite=TRUE)
  min.file<-paste(dir_rasters,statnames[n],"_min_",first.yr,"-",last.yr,"_",county.txt,".tif",sep="")
  print(min.file)
  writeRaster(min.r,file=min.file,overwrite=TRUE)
    
} # for n
