
# Files & Directories
# Set location for temporary directory
#options(rasterTmpDir='/home/ISAD/jm622/rscripts/temporary/')
#options(rasterTmpDir='/local/jm622/')

root<-"/home/ISAD/jm622/Data2015/"   # Source data and output data
in.root<-"/data/jm622/ModelData/" # model input data 

#root<-"~/Documents/Exeter/Data2015/"
#root<-"C:/Data2015/"
print(paste("Input root= ",in.root,sep=""))
print(paste("Output root= ",root,sep=""))
print(paste("Working directory= ",getwd(),sep=""))

print("Defining Directories")

# DEM files
dir_dem<-paste(in.root,"DEM/",sep="")

# Templates - 5km cells
dir_grids<-paste(root,"Templates/",sep="")

# Terrain
dir_terrain<-paste(in.root,"Terrain/",sep="")

# Temperature
dir_hrtemp<-paste(in.root,"Temp5km/hourly/",sep="") # dir for output files
dir_zip<-paste(root,"Temp5km/zip/",sep="")
dir_temp<-paste(root,"Temp5km/unzip/",sep="")

# Wind
dir_wind<-paste(in.root,"Wind/",sep="")
dir_shelter<-paste(in.root,"Wind/shelter/",sep="")

dir_windstrength<-paste(root,"Wind/strength/",sep="")
dir_winddirection<-paste(root,"Wind/direction/",sep="")
dir_windinvstr<-paste(root,"Wind/invstr/",sep="")

# Coast - landsea ratio and %land
dir_ldif<-paste(in.root,"CoastEffect/ldif/",sep="")
dir_coast<-paste(root,"CoastEffect/",sep="")
dir_percland<-paste(root,"CoastEffect/percland/",sep="")
dir_lsratio<-paste(root,"CoastEffect/lsratio/",sep="")
dir_upwindsea<-paste(root,"CoastEffect/upwindmaps/",sep="")

# SST
dir_ssth<-paste(in.root,"SST/hourly/",sep="")
dir_sst<-paste(root,"SST/",sep="")
#gunzip(paste(dir_sst,"HadISST_sst.nc.gz",sep=""))
dir_sstm<-paste(dir_sst,"monthly/",sep="")


# Radiation budget
dir_sisday<-paste(in.root,"CMSAF-SIS/day/",sep="")
dir_dniday<-paste(in.root,"CMSAF-DNI/day/",sep="")
dir_calday<-paste(in.root,"CMSAF-CAL/day/",sep="")
dir_lwr<-paste(root,"LWR/",sep="")

dir_dnitar<-paste(root,"CMSAF-DNI/tar/",sep="")
dir_dnigz<-paste(root,"CMSAF-DNI/ncgz/",sep="")
dir_sistar<-paste(root,"CMSAF-SIS/tar/",sep="")
dir_sisgz<-paste(root,"CMSAF-SIS/ncgz/",sep="")
dir_caltar<-paste(root,"CMSAF-CAL/tar/",sep="")
dir_calncgz<-paste(root,"CMSAF-CAL/ncgz/",sep="")

dir_sis<-paste(root,"CMSAF-SIS/extract/",sep="")
dir_dni<-paste(root,"CMSAF-DNI/extract/",sep="")
dir_cal<-paste(root,"CMSAF-CAL/extract/",sep="")

# Sea Pressure
dir_pressure<-paste(in.root,"SeaPressure/",sep="")

# Relative Humidity
dir_rh5km<-paste(in.root,"RelHumidity/rh5km/",sep="")
dir_rh<-paste(root,"RelHumidity/",sep="")

# Albedo
dir_albedo<-paste(in.root,"Albedo/",sep="")
dir_landsattar<-paste(root,"Albedo/landsat/tar/",sep="")
dir_landsatnc<-paste(root,"Albedo/landsat/extract/",sep="")

# Flow accumulation
dir.basinmap<-paste(in.root,"basins/",sep="")

# Other
dir_finalt<-paste(root,"Temp100m/",sep="")
dir_tmaps<-paste(root,"tmaps/",sep="")
dir_res5kmyear<-"/home/ISAD/jm622/results/year5km/"

##########################################################################################
# Outputs
#pdf(file="/home/ISAD/jm622/results/plots/Rplots.pdf",onefile=TRUE)
##########################################################################################
# Libraries
print("Defining libraries")
#install.packages(c('RNetCDF'))
#library(R.utils) # IMPORTANT - hides raster functions requiring raster::
library(RNetCDF)
library(ncdf4)
library(raster)
library(rgdal)
library(sp)
#library(chron)
library(insol) # required for julian day functions
library(mgcv) # require package for inputation (of CAL)
library(fields) # required for thin plate spline

##########################################################################################
# Common / General CONSTANTS used across other functions/programs
##########################################################################################
print("Defining Constants")
W.to.MJhr<-0.0036 # converting CMSAF satellite rad values (W/m2) to MJ/m2/hour required for calcs
adiabatic.lapserate<- 9.8/1000 # C per 1000m altitude diference 
albedo<-0.2
interval<-10 # used to determine number of wind sheltermaps and downwind land/sea maps calculated = 360/interval
radius<-10000 # used to define %land map and ldif calculations
buffer<-20000
latlong <- "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"
##########################################################################################
# Define Geographical extent of interest 
##########################################################################################
print("Defining DEM")

#  DEM files
ukdem.file<-paste(dir_dem,"demoriginal.tif",sep="")
dem.file<-paste(dir_dem,"dem.tif",sep="")
dembuf.file<- paste(dir_dem,"dembuf.tif",sep="")
grid5km.file<-paste(dir_dem,"grid5km.tif",sep="")
grid5kmbuf.file<-paste(dir_dem,"grid5kmbuf.tif",sep="")
land5km.file<-paste(dir_dem,"land5km.tif",sep="")
demland.file<-paste(dir_dem,"demland.tif",sep="")

# DEM rasters used by programs
demuk<-raster(ukdem.file)
dem<-raster(dem.file)
dembuf<-raster(dembuf.file)
grid5km.r<-raster(grid5km.file)
grid5kmbuf.r<-raster(grid5kmbuf.file)
land5km.r<-raster(land5km.file)
demland.r<-raster(demland.file)

# Slope & aspect rasters used - convert sea to 0 elevation to calculate slope/aspect
#r<-calc(dembuf,function(x){ifelse(is.na(x),0,x)})
#slope.r<-terrain(r,"slope", unit='degrees')
#slope.r<-mask(slope.r,dembuf)
#aspect.r<-terrain(r,"aspect",unit='degrees')
#aspect.r<-mask(aspect.r,dembuf)
#writeRaster(slope.r,paste(dir_terrain,"slope.tif",sep=""),overwrite=TRUE)
#writeRaster(aspect.r,paste(dir_terrain,"aspect.tif",sep=""),overwrite=TRUE)

##########################################################################################
# Common / General Functions used across other functions/programs
##########################################################################################
print("Defining Functions")
##########################################################################################
calcslope<-function(dem){
  r<-calc(dem,function(x){ifelse(is.na(x),0,x)})
  slope.r<-terrain(r,"slope")
  slope.r<-mask(slope.r,dem)
  return(slope.r)
}

##########################################################################################
calcaspect<-function(dem){
  r<-calc(dem,function(x){ifelse(is.na(x),0,x)})
  aspect.r<-terrain(r,"aspect")
  aspect.r<-mask(aspect.r,dem)
  return(aspect.r)
}
##########################################################################################
# RESAMPLING RASTER using thin plate spline
#library(fields) # required for thin plate spline
tps.resample<-function(input.r,output.r,maskoutput=TRUE){
  xy <- data.frame(xyFromCell(input.r,1:ncell(input.r)))
  v <- getValues(input.r)
  tps <- Tps(xy, v) # fit tps model (Don't worry about warning)
  result<- raster(output.r) # create blank raster at resolution of output.r
  
  # use model to predict values for all locations
  result<- interpolate(result,tps)
  if (maskoutput==TRUE) result<-mask(result,output.r)
  #plot(result,main="Thin-plate spline")
  
  return(result)
}# end function
##########################################################################################
# Julian date functions for use
#library(insol)

# Calculates julian date for 12:00 on the day of parameters
# Used in: setup, t5km_to_hourly_blocks
JDdmy<-function(day,month,year) # correct
{
  options(digits=12) 
  jdate<-insol::JD(ISOdate(year,month,day))
  return(jdate)
}

# Function to convert JD to day, month year as list - allows for jd fractions
# Gives Gregorian dates (even for julian period before 15C) 
# Inputs are vectors; Output as Data.frame
# Check with: http://aa.usno.navy.mil/faq/docs/JD_Formula.php
# Used in: setup, t5km_to_hourly_blocks
DMYjd<-function(jd) # correct 
{
  options(digits=12) 
  newdate<-as.POSIXlt(insol::JD(jd,inverse=TRUE))
  dmy<-data.frame(day=newdate$mday,month=newdate$mon+1,year=newdate$year+1900)
  return (dmy)
}

JDdoy<-function(doy,year)
{
  options(digits=12) 
  newdate<-insol::doyday(year,doy)
  jdate<-insol::JD(ISOdate(newdate$year+1900,newdate$mon+1,newdate$mday))
  return(jdate)
}
##########################################################################################
# Calculates days in month (from 1980-2015) from day, month, year
# Used to write monthly files of imputed values 
days.in.month<-function(day,month,year){
  y<-year-1979 # so 1=1980
  feb.d<-c(29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28,29,28,28,28,
           29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28) # days of Feb from 1980 to 2015
  monthdays<-c(31,feb.d[y],31,30,31,30,31,31,30,31,30,31)
  return(monthdays[month])
}
##########################################################################################
# Assign temperature values from 5km cells to 100m cells in block
# Function returns 100m raster from 5km raster for block
# Compare percentland_map_function & elevdif_map_function

ref5km.to.block100m<-function(dem.block,ref5km.r)
{
  ref.block<-resample(ref5km.r,dem.block,method="ngb")
  ref.block<-mask(ref.block,dem.block)
  #plot(ref.block,main="ref block")
  return(ref.block)
} # end function 
##########################################################################################
# PLOTTING FUNCTIONS
# Function for printing multiple maps of 1+ layers in same stack. 
# Argumennts: stack, vector of layer names, hour (for label) 
plot.stack<-function(stk,lyrs,hour) {
  # Define common colour scheme and scale for plots
  #par(mfrow=c(2,2))
  brk<-c(-50,0,100,200,300,400,500,600,700,800,900,1000)
  col<-rev(rainbow(11,start=1/6,end=4/6))
  for (map in 1:length(lyrs)){
    text<-paste(lyrs[map]," at ", hour, ":00",sep="")
    plot(x=stk,lyrs[map],col=col,breaks=brk, main="")
    title(main=text)
  } # for loop
} # function

# Plots tmp map rasters to fixed scale
temp.plot<-function(temp,hr,day,month,year) {
  t<-paste(day,"-",month,"-",year," ",hr,":00",sep="")
  brk<-c(-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)
  col<-rev(heat.colors(24))
  plot(temp,main=t,col=col,breaks=brk)
}

##########################################################################################
# Functions for identifying nearest 5km cell with historic data 
# Used by various e.g.: hourly_temperature function 

# Calculates number of 100m landcells in each 5km
landcells<-function(dem100m,ref5kmgrid.r,gridonly=TRUE)
{ 
  factor<-res(ref5kmgrid.r)[1] / res(dem100m)[1]
  land5km.r<-aggregate(dem100m,factor,fun=function(x,...)length(na.omit(x))) # all 5km cells
  #plot(land5km.r)
  if (gridonly){land5km.r<-mask(land5km.r,ref5kmgrid.r)} # only those 5km cells in ref5kmgrid.r
  #plot(land5km.r,main="Number of 100m cells")
  return(land5km.r) 
} # end function


# Find nearest coordinates in refx/y to coordinates x/y and return correspoinding ref value
# calculates mean value if several points at equal distance
nearestVal<-function(x,y,refx,refy,refvals)
{ 
  x<-rep(x,length(refx))
  y<-rep(y,length(refy))
  dist<-sqrt( abs(x-refx)^2 + abs(y-refy)^2 ) 
  sel<-which(dist==min(dist))
  if (length(sel)==0){warning(paste("Warning - no min found in nearLref: x=",x," y=",y," sel=",sel," min(dist)=",min(dist),sep=""))}
  if (length(sel)>1){refval<-mean(refvals[sel],na.rm=TRUE)
  } else {refval<-refvals[sel] }
  #print(refval)
  return(refval)
}


#Function: fill.5km.map
# USES FUNCTIONS landcells, nearestVal
# Input:
# ref5kmgrid.r = 5km frame covering same extent as 
# refdata5km.r = existing(incomplete) 5km data raster
# dem100m = 100m resolution land/sea raster (sea = NA, land >0)
# Output:
#   ref5km.filled.r  = filled ref 5km data taking values from nearest 5km cell to fill missing cells that contain 100m land cells

fill.5km.map<-function(refdata5km.r,ref5kmgrid.r) 
{
  if(compareRaster(refdata5km.r,ref5kmgrid.r)!=TRUE){warning("!!! 5km ref data and 5km grid do not match !!!")}
  # A. identify missing 5km cells (containing 100m land cells but without values)
  missing.r<-overlay(refdata5km.r,ref5kmgrid.r, fun=function(x,y){ifelse(is.na(x) & y>0,0,x)})
  #plot(missing.r,main="Missing cells (grey)")
  
  # Select xy and val of  cells without historic cover
  xyvals<-rasterToPoints(missing.r)
  sel<-which(xyvals[,3]==0)
  missingxyv<-xyvals[sel,1:3]
  
  # Find value from nearest cell with historic temp data
  sel<-which(!is.na(xyvals[,3]) & xyvals[,3]!=0) # selects land cells with data
  refvals<-xyvals[sel,1:3]
  
  # set vectors
  x<-missingxyv[,1]
  y<-missingxyv[,2]
  refx<-refvals[,1]
  refy<-refvals[,2]
  refvals<-refvals[,3]     
  for (i in 1: length(x)){
    missingxyv[i,3]<-nearestVal(x[i],y[i],refx,refy,refvals)
  }    
  # B. Assign missing %land values (0)  to missingxyv values  
  ref5km.filled.r<-rasterize(missingxyv[,1:2], refdata5km.r, missingxyv[,3], fun=max, update=TRUE)
  #plot(ref5km.filled.r,main="% land 5km cells - cells missing historic data set to nearest values")
  
  return (ref5km.filled.r)
} # end function

# FUNCTION to calculate number of NA in raster
rasterNA<-function(r){
  vals<-getValues(r)
  length(which(is.na(vals)))
}
} 
  
