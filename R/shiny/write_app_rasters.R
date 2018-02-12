# Prepare files required for shinyapp
# Reproject  and save in appropriate directories 
# REQUIRES: reprojected raster masks already created - in write_countymask_rasters
# REQUIRES: results rasters for whole of SW for year (in OS projection)

# Carson job variables
args <-commandArgs(trailingOnly = TRUE)
print(args)
year <- as.integer(args[1])
print(year)
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
dir_results<-paste(root,"Outputs/",sep="") # dir of results rasters for whole of SW
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
######################################
# Extract county rastersfor each statistic
######################################

# Define result years to be processed
years<-c(year)

# Results variables 
# Statnames used in results files from analyse_years_carson.R
statnames<-c("gdd10_gs","gdd5_gs","tmean_gs","tmin_year","tmax_year",
             "t20_gsdays","t25_gsdays","t30_gsdays", 
             "lastspfr_doy", "firstautfr_doy", "frostfree_days",
             "fl_tmean", "fl_numday")
labs<-c("Degree days","Degree days","degrees C","degrees C","degrees C","Number of days","Number of days","Number of days",
        "Day of Year (1:366)","Day of Year (1:366)","Number of days",
        "degrees C", "Number of days")

# Text below used in creation of rasters masks
# First set raster template in new projection for whole results area based on dem
  #e<-extent(dem)
  #template.r<-projectExtent(dem,crs=new.crs)
  #res(template.r)<-c(100,100)
  #print(template.r)

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

# 2 Crop and reproject results rasters

for (y in 1:length(years)){
  for (n in 1:length(statnames)){
    print(years[y])
    print(statnames[n])
    results.file<-paste(dir_results,statnames[n],"_",years[y],".tif",sep="")
    print(results.file)
    r<-raster(results.file) # raster for whole area
    crs(r)<-ukgrid # old crs (OS)
    #plot(r)
    for (c in 1:length(county)){
      print (county.txt[c])
      ### Reproject & write county rasters
      e<-extent(county[[c]])
      county.r<-projectRaster(r,to=county[[c]]) # reproject 
      county.r<-mask(crop(county.r,e), county[[c]])
      plot(county.r)
      # Write Raster -DECIDE OUTPUT DIR
      r.filename<-paste(dir_rasters,statnames[n],"_",year,"_",county.txt[c],".tif",sep="")
      print(r.filename)
      writeRaster(county.r,file=r.filename,overwrite=TRUE)
    } # for c
  } # for n
} # for y

 
