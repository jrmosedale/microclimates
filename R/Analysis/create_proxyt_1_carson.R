# CARSON - create proxy results
args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- as.integer(args[1])
start.month<-as.integer(args[2])
start.year<-as.integer(args[3] )
end.day<-as.integer(args[4] )
end.month<-as.integer(args[5] )
end.year<-as.integer(args[6] )

#######TEST VALUES
#end.day<-5
#end.month<-1
print(paste("Start: ",start.day,"/",start.month,"/",start.year,sep=""))
print(paste("End: ",end.day,"/",end.month,"/",end.year,sep=""))

root<-"/home/ISAD/jm622/Data2015/"   # Source data and output data
in.root<-"/data/jm622/ModelData/" # model input data 

dir_dem<-paste(in.root,"DEM/",sep="")
dir_proxyresults <-paste(root,"proxyt100/",sep="")
dir_terrain<-paste(in.root,"Terrain/",sep="")
dir_temp<-paste(root,"Temp5km/unzip/",sep="")
dir.basinmap<-paste(root,"basins/",sep="")
dir_results<-paste(root,"proxyt100/riskmaps/")
dir_counties<-paste(root,"counties/",sep="")
print(dir_counties)
# Output dir for risk rasters 
dir_results<-paste(dir_proxyresults,"riskmaps/",sep="")

# Create proxy t maps for cornwall at 100m resolution
library(raster)
library(insol)
library(fields) # required for thin plate spline
library(rgdal)
#library(rgeos)

##########################################################################################
# Common / General Functions used across other functions/programs
##########################################################################################
print("Defining Functions")
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
# Set values for all counties
lapserate<- 7/1000 # C per 1000m altitude diference 

start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
print(start.jd);print(end.jd)


dem.sw<-raster(paste(dir_dem,"dem.tif",sep=""))
slope.sw<-raster(paste(dir_terrain,"slope.tif",sep=""))
aspect.sw<-raster(paste(dir_terrain,"aspect.tif",sep=""))
twi.sw<-raster(paste(dir.basinmap,"topidx.tif",sep=""))
altdif.sw<-raster(paste(dir.basinmap,"altdif.tif",sep="")) # for cold air drainage
elevdif.sw<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # elevation effect


# for each COUNTY
counties<-c("cornwall","devon","dorset","somerset")
print(counties)

for (n in 1:length(counties)){
  print(counties[n])
  county.r<-raster(paste(dir_counties,counties[n],".tif",sep=""))
  
  # Crop all files to county
  dem<-crop(dem.sw,county.r)
  slope<-crop(slope.sw,county.r)
  aspect<-crop(aspect.sw,county.r)
  twi<-crop(twi.sw,county.r)
  altdif<-crop(altdif.sw,county.r)
  elevdif<-crop(elevdif.sw,county.r)
  
  valleyef<-(altdif*lapserate)*twi ; # plot(valleyef)


for (year in start.year:end.year){
  print(paste("Year is: ",year))
  numdays<-JDdmy(end.day,end.month,end.year)-JDdmy(start.day,start.month,start.year) +1
  tmax.filename<-paste(dir_proxyresults,counties[n],"-tmax-",year,".tif",sep="")
  tmin.filename<-paste(dir_proxyresults,counties[n],"-tmin-",year,".tif",sep="")
  tmean.filename<-paste(dir_proxyresults,counties[n],"-tmean-",year,".tif",sep="")
  print(tmax.filename); print(tmin.filename); print(tmean.filename)
  year.tmax<-brick(dem,values=FALSE,nl=numdays)
  year.tmin<-brick(dem,values=FALSE,nl=numdays)
  year.tmean<-brick(dem,values=FALSE,nl=numdays)
  
  for (jd in JDdmy(1,1,year) :JDdmy(end.day,end.month,year)){
    print(DMYjd(jd))
    # load day files 
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
    print(max.infile); print(min.infile)
    day.tmax<-crop(raster(max.infile, crs="+init=epsg:27700"),county.r)
    day.tmin<-crop(raster(min.infile, crs="+init=epsg:27700"),county.r)
    
    # interpolate to 100m
    day.tmax<-tps.resample(day.tmax,dem)
    day.tmin<-tps.resample(day.tmin,dem)
    
    # effect of cold air drainage  
    final.tmin<-day.tmin
    if ((DMYjd(jd)$month<5 | DMYjd(jd)$month>10) & cellStats(day.tmin,min)<3) {
      inversion<-runif(1)/3
      final.tmin<-overlay(day.tmin,valleyef,fun=function(x,y) { x-(inversion*y)   })
    }
    
    # correct for aspect
    final.tmax<-overlay(day.tmax,aspect,fun=function(x,y) {ifelse((y<210 & y>150),x+((x/3.5)*1/(1+exp(6-y)) ),x)  }) # south facing slopes
    final.tmax<-overlay(final.tmax,aspect,fun=function(x,y) {ifelse((y<30 | y>330),x-((x/5)*1/(1+exp(6-y)) ),x)  }) # north facing slopes
    
    # correct for elevation difference
    final.tmax<-final.tmax+(elevdif*lapserate)  # ; plot(day.tmax,main="day.max")
    final.tmin<-final.tmin+(elevdif*lapserate)
    final.tmean<-(final.tmin+final.tmax)/2
    
    # plot(final.tmax,main="finaltmax")
    # plot(final.tmin,main="tmin")
    # plot(final.tmax-day.tmax,main="dif")
    # plot(final.tmin-day.tmin)
    # plot(final.tmean,main="Tavg")
    
    # Save to  year brick
    layer<-jd-JDdmy(1,1,year)+1
    year.tmax[[layer]]<-final.tmax
    year.tmin[[layer]]<-final.tmin
    year.tmean[[layer]]<-final.tmean
    
  } # end doy
  
  writeRaster(year.tmax,filename=tmax.filename,format="GTiff",overwrite=TRUE)
  writeRaster(year.tmin,filename=tmin.filename,format="GTiff",overwrite=TRUE)
 
} # end year
} # end counties