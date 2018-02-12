# Create proxy t maps for cornwall at 100m resolution
library(insol)
library(fields) # required for thin plate spline
library(rgdal)
library(rgeos)
root<-"~/Documents/Exeter/Data2015/"
dir_proxyresults <-"~/Documents/Exeter/Data2015/proxyt100/"
dir_terrain<-paste(root,"Terrain/",sep="")
dir_temp<-paste(root,"Temp5km/unzip/",sep="")
dir.basinmap<-paste(root,"basins/",sep="")

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

##########################################################################################
dem<-raster(paste(dir_dem,"dem.tif",sep="")); 
slope<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem)
aspect<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem)
twi<-raster(paste(dir.basinmap,"topidx.tif",sep=""))
altdif<-raster(paste(dir.basinmap,"altdif.tif",sep="")) # for cold air drainage
elevdif<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # elevation effect
# set terrrain to cornwall incl dem

hcounties.shp<- readOGR(dsn = path.expand(paste(root,"OSdata/bdline_essh_gb/Data/Supplementary_Historical",sep="")), layer = "Boundary-line-historic-counties_region")
# extract county polygons - use historic counties data 
sel <- which(hcounties.shp$Name == "Cornwall")
cornwall<-hcounties.shp[sel,]
e<-extent(xmin(cornwall),xmax(cornwall),ymin(cornwall),ymax(cornwall))


dem<-mask(crop(dem,e), cornwall); #plot(dem)
slope<-mask(crop(slope,e), cornwall); #plot(slope)
aspect<-mask(crop(aspect,e), cornwall); #plot(aspect)
twi<-mask(crop(twi,e), cornwall); #plot(twi)
altdif<-mask(crop(altdif,e), cornwall); #plot(altdif)
elevdif<-mask(crop(elevdif,e), cornwall); #plot(elevdif)

lapserate<- 7/1000 # C per 1000m altitude diference 
valleyef<-(altdif*lapserate)*twi ; # plot(valleyef)

start.day <- 1
start.month<-1
start.year<-2011
print(paste("Start: ",start.day,"/",start.month,"/",start.year,sep=""))
end.day<-31
end.month<-12
end.year<-2011
print(paste("End: ",end.day,"/",end.month,"/",end.year,sep=""))

start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
print(start.jd);print(end.jd)
#  year<-start.year
#  jd<-JDdmy(1,1,year)

for (year in start.year:end.year){
  print(paste("Year is: ",year))
  numdays<-JDdmy(31,12,year)-JDdmy(1,1,year) +1
  year.tmax<-brick(dem,values=FALSE,nl=numdays, filename=paste(dir_proxyresults,"tmax_cornwall_",year,".tif",sep=""))
  year.tmin<-brick(dem,values=FALSE,nl=numdays, filename=paste(dir_proxyresults,"tmin_cornwall_",year,".tif",sep=""))
  year.tmean<-brick(dem,values=FALSE,nl=numdays, filename=paste(dir_proxyresults,"tmean_cornwall_",year,".tif",sep=""))
  
  for (jd in JDdmy(1,1,year) :JDdmy(31,12,year)){
     print(DMYjd(jd))
     # load day files 
      max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
      min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
      print(max.infile); print(min.infile)
      day.tmax<-crop(raster(max.infile, crs="+init=epsg:27700"),e)
      day.tmin<-crop(raster(min.infile, crs="+init=epsg:27700"),e)
      
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
      final.tmax<-final.tmax+(elevdif.r*lapserate)  # ; plot(day.tmax,main="day.max")
      final.tmin<-final.tmin+(elevdif.r*lapserate)
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
  #writeRaster(year.tmax,filename=path.expand(paste(dir_proxyresults,"tmax_cornwall_",year,".tif",sep="")),format="GTiff",overwrite=TRUE)
  #writeRaster(year.tmin,filename=paste(dir_proxyresults,"tmin_cornwall_",year,".tif",sep=""),overwrite=TRUE)
  #writeRaster(year.tmean,filename=paste(dir_proxyresults,"tmean_cornwall_",year,".tif",sep=""),overwrite=TRUE)
  
} # end year



####################################################################################
# Calculate Seasonal Rasters 
####################################################################################
# Define Growing Season (or part of year of interest)
#start.gs<-1; end.gs<-nlayers(year.tmin)

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



####################################################################################
# Output raster files for year - seasonal characteristics by block
####################################################################################
dir_results<-"~/Documents/Exeter/Data2015/proxyt100/riskmaps/"

gdd10.fileout<-paste(dir_results,"gdd10-",year,".tif" ,sep="")
gdd5.fileout<-paste(dir_results,"gdd5-",year,".tif" ,sep="")
meant.fileout<-paste(dir_results,"meant-",year ,".tif" ,sep="")
maxt.fileout<-paste(dir_results,"maxt-",year ,".tif" ,sep="")
days20.fileout<-paste(dir_results,"days20-",year ,".tif" ,sep="")
days25.fileout<-paste(dir_results,"days25-",year ,".tif" ,sep="")
days30.fileout<-paste(dir_results,"days30-",year ,".tif" ,sep="")
mint.fileout<-paste(dir_results,"mint-",year ,".tif" ,sep="")
spfrost.fileout<-paste(dir_results,"spfrost-",year ,".tif" ,sep="")
autfrost.fileout<-paste(dir_results,"autfrost-",year ,".tif" ,sep="")
frostfree.fileout<-paste(dir_results,"frostfree-",year ,".tif" ,sep="")


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

