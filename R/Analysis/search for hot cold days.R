# Identify particular weather conditions from 5km data etc

start.day <- 1
start.month<-1
start.year<-1992
print(paste("Start: ",start.day,"/",start.month,"/",start.year,sep=""))
end.day<-23
end.month<-06
end.year<-1992
print(paste("End: ",end.day,"/",end.month,"/",end.year,sep=""))

start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
print(start.jd);print(end.jd)


for (jd in start.jd:end.jd){
  year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day ;# values for t
  
 
  # Load daily files of hourly values of temperature, RH, lwr
  
  max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
  min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
  day.tmax<-raster(max.infile, crs="+init=epsg:27700")
  day.tmin<-raster(min.infile, crs="+init=epsg:27700")
  day.tmax<-crop(x=day.tmax,y=dem)
  day.tmin<-crop(x=day.tmin,y=dem)
  #day.tmean<-overlay(day.tmax,day.tmin,fun=function(x,y){return((x+y)/2)})
  tmax.mean<-cellStats(day.tmax,stat='mean',na.rm=TRUE)
  tmin.mean<-cellStats(day.tmin,stat='mean',na.rm=TRUE)
  if (tmax.mean>25) print(paste("T>28C: ",day,"/",month,"/",year," Tmean= ",tmax.mean,sep=""))
  #if ( tmin.mean<=0) print(paste("T<=0C: ",day,"/",month,"/",year," Tmean= ",tmin.mean,sep=""))
  
  
}



filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
dnr.24h<-stack(filein1) ; projection(dnr.24h)<-CRS("+init=epsg:4326")
sis.24h<-stack(filein2) ; projection(sis.24h)<-CRS("+init=epsg:4326")
dnr.24h<-projectRaster(dnr.24h,crs="+init=epsg:27700")
sis.24h<-projectRaster(sis.24h,crs="+init=epsg:27700")


# Load and reproject CAL data to OS (downscaling by resampling takes place in hourly interval)
cal.filein<-paste(dir_calday,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
print(paste("CAL file in: ",cal.filein,sep=""))
cal.day<-brick(cal.filein)
cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")

# Cold days 2009
day <- 7; month<-1 ; year<-2009 # 1-11 all below 0
day <- 3; month<-2 ; year<-2009 # 2-4/2
day <- 19; month<-12 ; year<-2009 # 18-23 all below 0
# Cold days 2010
day <- 7; month<-1 ; year<-2010 # 1-10 all below 0 RAD

# Late cold days 1992
day <- 4; month<-4 ; year<-1992 # 4-5 all below 0 RAD also
17/10/

### Hot days >25C
day <- 1; month<-7 ; year<-2009 # 1-11 all below 0



###################################################
# CHECK COld inversion for 4/41992 - CAL and Wstr
###################################################
day <- 4; month<-4 ; year<-1992 # 4-5 all below 0 RAD also
jd<-JDdmy(day,month,year)
ukcpcell<-810

gridmask.r<-land5km.r # land5km.r defined in setup
vals<-getValues(gridmask.r)
xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
sel<-which(vals==1)
landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell
print(paste("UKCP cell= ", ukcpcell,sep=""))
# Crop dem to fit cell and buffered cell
x<-landcells[ukcpcell,1]
y<-landcells[ukcpcell,2]
e.block<-extent(x-2500,x+2500,y-2500,y+2500)
dem.block<-crop(demuk,e.block) ; print (dem.block)
e.buffer<-extent(x-2500-buffer,x+2500+buffer,y-2500-buffer,y+2500+buffer)
dem.buffer<-crop(demuk,e.buffer)
# plot(dem.buffer,main=ukcpcell);plot(dem.block,main=ukcpcell)
print(paste("Plotting location of cell ",ukcpcell,sep=""))
cell_location(gridmask.r,ukcpcell) # plot map of cell location
 
  
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))
# Load and reproject CAL data to OS (downscaling by resampling takes place in hourly interval)
cal.filein<-paste(dir_calday,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
print(paste("CAL file in: ",cal.filein,sep=""))
cal.day<-brick(cal.filein)
cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")

# Load daily files of hourly values of temperature, RH, lwr
t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_100m.r", sep="")
load(t.filein) #t100m.day

# load 5km RH data
rh.filein<-paste(dir_rh5km,"RH_100m_",year,"_",month,"_",day,".R",sep="")
load(rh.filein) # rh.day

par(mfrow=c(2,2))
for (hr in 0:23){
  tref.block<-crop(raster(t100m.day[,,hr+1],template=dem),dem.block)
  rhref.block<-crop(raster(rh.day[,,hr+1],template=dem),dem.block)
  cal.block<-tps.resample(crop(raster(cal.day,layer=hr+1),dem.buffer),dem.block)
  cal.block<-calc(cal.block,fun=function(x)ifelse(x<0,0,x))
  
  wind.results<-wind.tpsdownscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval,print.res=FALSE,write.files=FALSE )
  wdir.block<-wind.results[[1]]
  wstr.block<-wind.results[[2]]
  
  lwr.block<-calc_lwr_block(jd,cal.block,rhref.block,tref.block,plotlwr=FALSE,writefile=FALSE)
  plot(tref.block,title=paste("Tref ",hr,"h00",sep=""))
  plot(cal.block,title=paste("CAL ",hr,"h00",sep=""))
  plot(lwr.block,title=paste("LWR ",hr,"h00",sep=""))
  plot(wstr.block,title=paste("WSTR ",hr,"h00",sep=""))
  
}
