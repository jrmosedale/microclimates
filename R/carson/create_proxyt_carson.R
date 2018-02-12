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
start.day<-1;start.month<-1; start.year<-2011
end.day<-31; end.month<-12; end.year<-2011
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
dir_results<-paste(root,"proxyt100/riskmaps/",sep="")
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
slope.sw<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem.sw)
aspect.sw<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem.sw)
twi.sw<-raster(paste(dir.basinmap,"topidx.tif",sep=""))
altdif.sw<-raster(paste(dir.basinmap,"altdif.tif",sep="")) # for cold air drainage
elevdif.sw<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # elevation effect
compareRaster(dem.sw,slope.sw,aspect.sw,twi.sw,altdif.sw,elevdif.sw)

# for each COUNTY
counties<-c("cornwall","devon","dorset","somerset")
####   TEST ONLY   ###
counties<-c("test")
print(counties)

for (n in 1:length(counties)){
  print(counties[n])
  # TEST ONLY
  #county.r<-raster(paste(dir_counties,counties[n],".tif",sep=""))
  county.r<-crop(dem.sw,extent(180000,200000,65000,85000)); plot(county.r)
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
  year.tmax<-stack(dem,values=FALSE,nl=numdays)
  year.tmin<-stack(dem,values=FALSE,nl=numdays)
  year.tmean<-stack(dem,values=FALSE,nl=numdays)
  
  for (jd in JDdmy(1,1,year) :JDdmy(end.day,end.month,year)){
    print(DMYjd(jd))
    # load day files 
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
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
  
  #writeRaster(year.tmax,filename=tmax.filename,format="GTiff",overwrite=TRUE)
  #writeRaster(year.tmin,filename=tmin.filename,format="GTiff",overwrite=TRUE)
  #writeRaster(year.tmean,filename=tmean.filename,format="GTiff",overwrite=TRUE)
  
  ####################################################################################
  # Calculate Seasonal Rasters 
  ####################################################################################
  # 1. Calculate last spring and first fall frost - not limited to growing season
  # Calc last spring frost <= 1 C 
  # Assumes no frost between end of May & early September - so reduces vector to 1-150 doy
  spfrdata.s<-subset(year.tmin,1:150)
  start.v<-rep(1,(nlayers(spfrdata.s)))
  spfrost.r<-calc(spfrdata.s,fun=function(x){ifelse(length(which(x<=1))>0,tail(which(x<=1),1)+start.v,1)}) # extract layer of last frost day
  spfrost.r<-mask(spfrost.r,dem,maskvalue=NA)
  plot(spfrost.r,main=paste("Last spring frost day ",year,sep=""))
  
  # Calculate first autumn frost (after early sept)
  autfrdata.s<-subset(year.tmin,240:nlayers(year.tmin))
  start.v<-rep(240,(nlayers(autfrdata.s)))
  autfrost.r<-calc(autfrdata.s,fun=function(x){ifelse(length(which(x<=1))>0,head(which(x<=1),1)+start.v,nlayers(year.tmin))}) # extract layer of last frost day
  autfrost.r<-mask(autfrost.r, dem,maskvalue=NA)
  plot(autfrost.r,main=paste("First autumn frost day ",year,sep=""))
  
  # Calculate frost free period
  frostfree.r<-overlay(spfrost.r,autfrost.r,fun=function(x,y){return(y-x)})
  plot(frostfree.r,main=paste("Frost free period of year ",year,sep=""))
  
  # 2. Calculate growing season stats - temperature extremes, gdd etc
  start.gs<-90; end.gs<-305
  
  Tbase<-10;  tbase.v<-rep(Tbase,(end.gs-start.gs+1))
  gdd10.r<-calc(subset(year.tmean,start.gs:end.gs),fun=function(x){sum(x-tbase.v)})
  gdd10.r<-mask(gdd10.r, dem,maskvalue=NA)
  plot(gdd10.r,main=paste("GDD10 ",year," Tbase= ",Tbase,sep=""))
  
  Tbase<-5;  tbase.v<-rep(Tbase,(end.gs-start.gs+1))
  gdd5.r<-calc(subset(year.tmean,start.gs:end.gs),fun=function(x){sum(x-tbase.v)})
  plot(gdd5.r,main=paste("GDD5 ",year," Tbase= ",Tbase,sep=""))
  
  # Calc mean T
  meangst.r<-calc(subset(year.tmean,start.gs:end.gs),fun=function(x){mean(x)})
  plot(meangst.r,main=paste("MeanT ",year,sep=""))
  
  # Calc max T
  maxgst.r<-calc(subset(year.tmax,start.gs:end.gs),fun=function(x){max(x)})
  plot(maxgst.r,main=paste("Max T ",year,sep=""))
  
  # Calc #days where max temp>  20C, 25C, 30C from April-Oct
  days20.r<-calc(subset(year.tmax,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=20))>0,length(which(x>=20)),0)} )
  days20.r<-mask(days20.r,dem,maskvalue=NA)
  
  days25.r<-calc(subset(year.tmax,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=25))>0,length(which(x>=25)),0)} )
  days25.r<-mask(days25.r, dem,maskvalue=NA)
  
  days30.r<-calc(subset(year.tmax,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=30))>0,length(which(x>=30)),0)} )
  days30.r<-mask(days30.r, dem,maskvalue=NA)
  
  plot(days20.r,main=paste("Days=>20 ",DMYjd(jd)$year,sep=""))
  plot(days25.r,main=paste("Days=>25 ",DMYjd(jd)$year,sep=""))
  plot(days30.r,main=paste("Days=>30 ",DMYjd(jd)$year,sep=""))
  
  # Calc min T in growing season
  mingst.r<-calc(subset(year.tmin,start.gs:end.gs),fun=function(x){min(x)})
  plot(mingst.r,main=paste("Min T ",DMYjd(jd)$year,sep=""))
  
  # Calculate flowering risks
  # Calculate days of flowering
  fl.start<-daydoy(year,6,22) # doy
  fl.end<-daydoy(year,7,5) # doy
  # define critical mean t
  t.fl<-15
  fltmean.r<- calc(subset(year.tmean,fl.start:fl.end),fun=function(x){mean(x)})
  plot(fltmean.r)
  
  # Number of flowering days where mean T>15C
  flday<-fl.start
  fl.numdays.r<-raster(year.tmean,layer=1)
  fl.numdays.r<-calc(fl.numdays.r,function(x){ifelse(!is.na(x),0,NA)})
  for(day in fl.start:fl.end){
    print(day)
    #plot(subset(year.tmean,day:day))
    good.day<- calc(raster(year.tmean,layer=day),fun=function(x){ifelse(mean(x)>t.fl,1,0)} )
    fl.numdays.r<-fl.numdays.r+good.day
    #plot(raster(apply(tmodel.r[,,flhr:(flhr+23)], c(1,2), function(x) {mean(x)} ),template=dem.block),main="Good day?")
  } # end for
  plot(fl.numdays.r)
  
  ####################################################################################
  # Output raster files for year
  ####################################################################################
  county<-counties[n]
  
  gdd10.fileout<-paste(dir_results,"gdd10_gs_",year,"_",county,".tif" ,sep="")
  gdd5.fileout<-paste(dir_results,"gdd5_gs_",year,"_",county,".tif" ,sep="")
  meant.fileout<-paste(dir_results,"tmean_gs_",year,"_",county ,".tif" ,sep="")
  maxt.fileout<-paste(dir_results,"tmax_year_",year,"_",county ,".tif" ,sep="")
  mint.fileout<-paste(dir_results,"tmin_year_",year,"_",county ,".tif" ,sep="")
  days20.fileout<-paste(dir_results,"t20_gsdays_",year,"_",county ,".tif" ,sep="")
  days25.fileout<-paste(dir_results,"t25_gsdays_",year,"_",county ,".tif" ,sep="")
  days30.fileout<-paste(dir_results,"t30_gsdays_",year,"_",county ,".tif" ,sep="")
  spfrost.fileout<-paste(dir_results,"lastspfr_doy_",year,"_",county ,".tif" ,sep="")
  autfrost.fileout<-paste(dir_results,"firstautfr_doy_",year,"_",county ,".tif" ,sep="")
  frostfree.fileout<-paste(dir_results,"frostfree_days_",year,"_",county ,".tif" ,sep="")
  fltmean.fileout<-paste(dir_results,"fl_tmean_",year,"_",county ,".tif" ,sep="")
  flnumday.fileout<-paste(dir_results,"fl_numday_",year,"_",county ,".tif" ,sep="")
  
  writeRaster(gdd10.r,file=gdd10.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(gdd5.r,file=gdd5.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(meangst.r, file=meant.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(maxgst.r,file=maxt.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(days20.r,file=days20.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(days25.r,file=days25.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(days30.r,file=days30.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(mingst.r,file=mint.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(spfrost.r,file=spfrost.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(autfrost.r,file=autfrost.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(frostfree.r,file=frostfree.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(fltmean.r,file=fltmean.fileout,format="GTiff",overwrite=TRUE)
  writeRaster(flnumdays.r,file=flnumday.fileout,format="GTiff",overwrite=TRUE)
  
} # end year
} # end counties