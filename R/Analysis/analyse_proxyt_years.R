



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

# Calculate flowering risks
# Calculate days of flowering
fl.start<-daydoy(year,6,22) # doy
fl.end<-daydoy(year,7,5) # doy
# define critical mean t
t.fl<-15
fltmean.r<- calc(subset(year.tmean,start.fl:end.fl),fun=function(x){mean(x)})

# Number of flowering days where mean T>15C
flday<-fl.start
fl.numday<-array(rep(0,ncell(year.tmean)),dim=c(nrow(year.tmean),ncol(year.tmean)) )
for(day in fl.start:fl.end){
  good.day<- calc(subset(year.tmean,day:day),fun=function(x){ifelse(mean(x)>t.fl,1,0)} )
  fl.numday<-fl.numday+good.day
  #plot(raster(apply(tmodel.r[,,flhr:(flhr+23)], c(1,2), function(x) {mean(x)} ),template=dem.block),main="Good day?")
} # end for
flnumdays.r<-raster(as.matrix(fl.numday),template=year.tmean)

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

# List of statistics
statistics<-list(gdd10_gs,gdd5_gs,tmean_gs,tmin_year,tmax_year,
                 t20_gsdays,t25_gsdays,t30_gsdays, 
                 lastspfr_doy, firstautfr_doy, frostfree_days,
                 fl_tmean, fl_numday)
statnames<-c("gdd10_gs","gdd5_gs","tmean_gs","tmin_year","tmax_year",
             "t20_gsdays","t25_gsdays","t30_gsdays", 
             "lastspfr_doy", "firstautfr_doy", "frostfree_days",
             "fl_tmean", "fl_numday")
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