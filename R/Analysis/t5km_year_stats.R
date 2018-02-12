# Calculate key temperature parameters for CV vineyard
# Carson job variables
args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- args[1] 
start.month<-args[2] 
start.year<-args[3] 
end.day<-args[4] 
end.month<-args[5] 
end.year<-args[6] 
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file (jd functions, dem etc)
dir_results<-dir_res5kmyear

## ASSUME start end end correspond to start and end of single YEAR or GROWING SEASON
#start.day<-1; start.month<-1; start.year<-2008
#end.day<-31; end.month<-12; end.year<-2008
# dir_temp<-"C:/Data2015/Temp5km/extract/"
# dir_results<-"C:/Results2015/year_stats_5km/"
# dir_allyr<-"C:/Results2015/allyear_stats_5km"
# Uses 5km grid mask: grid5km.r
# dir_grids<-"C:/Data2015/Templates/"


start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)
print(start.jd)
print(end.jd)

if (DMYjd(start.jd)$year!=DMYjd(end.jd)$year) print("Data from more than one year !!!!!!") else year<-DMYjd(start.jd)$year
print(paste("Analysing data for Year= ",year,sep=""))

# CREATE YEAR STACK OF TEMP DATA FOR WHOLE AREA - layers 
# create empty stack
tmax.s<-stack()
tmin.s<-stack()
tmean.s<-stack()

# Load 5km data for entire time period
for (jd in start.jd:end.jd){
# Read day temperature data
max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
day.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
day.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
day.tmax<-crop(x=day.tmax,y=dem)
day.tmin<-crop(x=day.tmin,y=dem)
day.tmean<-overlay(day.tmax,day.tmin,fun=function(x,y){return((x+y)/2)})
tmax.s<-stack(tmax.s,day.tmax)
tmin.s<-stack(tmin.s,day.tmin)
tmean.s<-stack(tmean.s,day.tmean)
}

####################################################################################
# Calculate Seasonal Rasters 
####################################################################################

# 1. Calculate last spring and first fall frost - not limited to growing season
# Calc last spring frost <= 1 C 
# Assumes no frost between end of May & early September - so reduces vector to start.jd - last day of jd =150 
spfrdata.s<-subset(tmin.s,1:150)
start.v<-rep(start.gs,(nlayers(spfrdata.s)))
spfrost.r<-calc(spfrdata.s,fun=function(x){ifelse(length(which(x<=2))>0,tail(which(x<=2),1)+start.v,1)}) # extract layer of last frost day
spfrost.r<-mask(spfrost.r,grid5km.r,maskvalue=NA)
plot(spfrost.r,main=paste("Last spring frost day ",DMYjd(jd)$year,sep=""))

# Calculate first autumn frost (after early sept)
autfrdata.s<-subset(tmin.s,240:nlayers(tmin.s))
start.v<-rep(240,(nlayers(autfrdata.s)))
autfrost.r<-calc(autfrdata.s,fun=function(x){ifelse(length(which(x<=2))>0,head(which(x<=2),1)+start.v,nlayers(tmin.s))}) # extract layer of last frost day
autfrost.r<-mask(autfrost.r,grid5km.r,maskvalue=NA)
plot(autfrost.r,main=paste("First autumn frost day ",DMYjd(jd)$year,sep=""))

# Calculate frost free period
frostfree.r<-overlay(spfrost.r,autfrost.r,fun=function(x,y){return(y-x)})
plot(frostfree.r,main=paste("Frost free period of year ",DMYjd(jd)$year,sep=""))


# 2. Calculate growing season stats - temperature extremes, gdd etc
# Define Growing Season (or part of year of interest)
start.gs<-1; end.gs<-nlayers(tmean.s)
#start.gs<-90; end.gs=305
# correct start and end dates to reflect zone of interest
end.jd<-start.jd+end.gs-1
start.jd<-start.jd+start.gs-1
# OR SET GS to first and last frosts??

# Calc gdd 
Tbase<-10;  tbase.v<-rep(Tbase,(end.gs-start.gs+1))
gdd10.r<-calc(subset(tmean.s,start.gs:end.gs),fun=function(x){sum(x-tbase.v)})
plot(gdd10.r,main=paste("GDD10 ",DMYjd(jd)$year," Tbase= ",Tbase,sep=""))

Tbase<-5;  tbase.v<-rep(Tbase,(end.gs-start.gs+1))
gdd5.r<-calc(subset(tmean.s,start.gs:end.gs),fun=function(x){sum(x-tbase.v)})
plot(gdd5.r,main=paste("GDD5 ",DMYjd(jd)$year," Tbase= ",Tbase,sep=""))

# Calc mean T
meangst.r<-calc(subset(tmean.s,start.gs:end.gs),fun=function(x){mean(x)})
plot(meangst.r,main=paste("MeanT ",DMYjd(jd)$year,sep=""))

# Calc max T
maxgst.r<-calc(subset(tmax.s,start.gs:end.gs),fun=function(x){max(x)})
plot(maxgst.r,main=paste("Max T ",DMYjd(jd)$year,sep=""))

# Calc #days where max temp>  20C, 25C, 30C from April-Oct
days20.r<-calc(subset(tmax.s,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=20))>0,length(which(x>=20)),0)} )
days20.r<-mask(days20.r,grid5km.r,maskvalue=NA)

days25.r<-calc(subset(tmax.s,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=25))>0,length(which(x>=25)),0)} )
days25.r<-mask(days25.r,grid5km.r,maskvalue=NA)

days30.r<-calc(subset(tmax.s,start.gs:end.gs),fun=function(x){ifelse(length(which(x>=30))>0,length(which(x>=30)),0)} )
days30.r<-mask(days30.r,grid5km.r,maskvalue=NA)

plot(days20.r,main=paste("Days=>20 ",DMYjd(jd)$year,sep=""))
plot(days25.r,main=paste("Days=>25 ",DMYjd(jd)$year,sep=""))
plot(days30.r,main=paste("Days=>30 ",DMYjd(jd)$year,sep=""))

# Calc min T in growing season
mingst.r<-calc(subset(tmin.s,start.gs:end.gs),fun=function(x){min(x)})
plot(mingst.r,main=paste("Min T ",DMYjd(jd)$year,sep=""))

####################################################################################
# Output raster files for year - seasonal characteristics for every 5km grid cell
####################################################################################
gdd10.fileout<-paste(dir_results,"gdd10_5km_",year,".tif" ,sep="")
gdd5.fileout<-paste(dir_results,"gdd5_5km_",year,".tif" ,sep="")
meant.fileout<-paste(dir_results,"meant_5km_",year ,".tif" ,sep="")
maxt.fileout<-paste(dir_results,"maxt_5km_",year ,".tif" ,sep="")
days20.fileout<-paste(dir_results,"days20_5km_",year ,".tif" ,sep="")
days25.fileout<-paste(dir_results,"days25_5km_",year ,".tif" ,sep="")
days30.fileout<-paste(dir_results,"days30_5km_",year ,".tif" ,sep="")
mint.fileout<-paste(dir_results,"mint_5km_",year ,".tif" ,sep="")
spfrost.fileout<-paste(dir_results,"spfrost_5km_",year ,".tif" ,sep="")
autfrost.fileout<-paste(dir_results,"autfrost_5km_",year ,".tif" ,sep="")
frostfree.fileout<-paste(dir_results,"frostfree_5km_",year ,".tif" ,sep="")

writeRaster(gdd10.r,file=gdd10.fileout,format="GTiff")
writeRaster(gdd5.r,file=gdd5.fileout,format="GTiff")
writeRaster(maxt.r,file=maxt.fileout,format="GTiff")
writeRaster(days20.r,file=days20.fileout,format="GTiff")
writeRaster(days25.r,file=days25.fileout,format="GTiff")
writeRaster(days30.r,file=days30.fileout,format="GTiff")
writeRaster(mint.r,file=mint.fileout,format="GTiff")
writeRaster(spfrost.r,file=spfrost.fileout,format="GTiff")
writeRaster(autfrost.r,file=autfrost.fileout,format="GTiff")
writeRaster(frostfree.r,file=frostfree.fileout,format="GTiff")



