# PURPOSE: Load existing 5km daily stack of hourly data and interpolate to 100m

args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- as.integer(args[1])
start.month<-as.integer(args[2])
start.year<-as.integer(args[3] )
end.day<-as.integer(args[4] )
end.month<-as.integer(args[5] )
end.year<-as.integer(args[6] )

#start.day<-1; start.month<-7; start.year<-1992
#end.day<-4; end.month<-7; end.year<-1992

print("Calling set up...")
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

start.jd<-JDdmy(start.day,start.month,start.year); print(paste("Start JD= ",start.jd,sep=""))
end.jd<-JDdmy(end.day,end.month,end.year); print(paste("End JD= ",end.jd,sep=""))

##########################################################################################
# FUNCTIONS USED TO GENERATE HOURLY TEMPERATURE DATA AT 5KM RESOLUTION
#######################################################################################

# Calculates sunrise and sunset (set Timezone and DST to a vector of 0)
# uses julian day input rather than doy and year
# CHECK - how jd varies from start to end of day
# Input / output as vectors of same length
suntimes.grid<-function(JD,Lat,Long,Timezone,DST){
  J<-JD
  lw<-Long*-1
  n<-J-2451545-0.0009-(lw/360)
  n<-floor(n)+0.5
  sn<-2451545+0.0009+(lw/360)+n
  msa<-(357.5291+0.98560028*(sn-2451545))%%360
  eoc<-1.9148*sin(msa*pi/180)+0.02*sin(2*msa*pi/180)+0.0003*sin(3*msa*pi/180)
  ecl<-(msa+102.9372+eoc+180); ecl<-ecl%%360
  st<-sn+(0.0053*sin(msa*pi/180))-(0.0069*sin(2*ecl*pi/180))
  d<-asin(sin(ecl*pi/180)*sin(23.45*pi/180))
  cos.has<-((sin(-0.83*pi/180)-sin(Lat*pi/180)*sin(d))/(cos(Lat*pi/180)*cos(d)))
  h.set<-vector(length=length(JD)); h.rise<-h.set; dl<-h.set
  # next three lines may provoke warnings in some cases?
  has<-acos(cos.has)
  J.set<-2451545+0.0009+(((has*180/pi+lw)/360)+n+0.0053*sin(msa*pi/180))-0.0069*sin(2*ecl*pi/180)
  J.rise<-st-(J.set-st)
  
  ifelse(cos.has^2<1,
         h.set<-(J.set%%1)*24+Timezone+DST,
         ifelse(cos.has>1,
                h.set<-12,
                ifelse(cos.has<(-1),
                       h.set<-0,
                       warning("Cos.has case not found") ) ) )
  
  ifelse(cos.has^2<1,
         h.rise<-(J.rise%%1)*24+Timezone+DST,
         ifelse(cos.has>1,
                h.rise<-12,
                ifelse(cos.has<(-1),
                       h.rise<-0,
                       warning("Cos.has case not found") ) ))
  
  ifelse(cos.has^2<1, dl<-(J.set-J.rise)*24,
         ifelse(cos.has>1,
                dl<-0,
                ifelse(cos.has<(-1),
                       dl<-24,
                       print("Cos.has case not found") ) ) )    
  
  if(any(dl==0)) {warning("sun below horizon for 24 hours")}   
  if(any(dl==24)) {warning("sun above horizon for 24 hours")}    
  
  sun.vars<-data.frame(sunrise=h.rise,sunset=h.set,daylight=dl)
  sun.vars
}
# calculates sunrise
sunrise.grid<-function(JD,Lat,Long,Timezone,DST){
  sun.rise<-suntimes.grid(JD-1,Lat,Long,Timezone,DST)[,1]
  sun.rise
}
# calculates sunset
sunset.grid<-function(JD,Lat,Long,Timezone,DST){
  sun.set<-suntimes.grid(JD,Lat,Long,Timezone,DST)[,2]
  sun.set
}
# generates hourly values for use between dawn and hotest part of day (~13:35)
generate.hrtemps.grid.day1<-function(min.temp,max.temp,day.length)
{
  x<-matrix(rep(0:23,NROW(min.temp)),nrow=NROW(min.temp), ncol=24, byrow=TRUE)
  A=(max.temp-min.temp)/2
  fr=1/(day.length*1.5)
  phase=13.58989
  of=A+min.temp
  y<-A*cos(2*pi*fr*(x-phase))+of
  y
}
# generates hourly values for use between hotest part of day (~13:35) and midnight
generate.hrtemps.grid.day2<-function(min.temp,max.temp,sun.rise)
{
  x<-matrix(rep(0:23,NROW(min.temp)),nrow=NROW(min.temp), ncol=24, byrow=TRUE)
  A=(max.temp-min.temp)/2
  lengt<-(24-13.58989)+sun.rise
  fr=1/(lengt*1.5)
  phase=13.58989
  of=A+min.temp
  y<-A*cos(2*pi*fr*(x-phase))+of
  y
}
# generates hourly values between midnight and dawn. Coolest part of day ~12 mins before dawn
generate.hrtemps.grid.night<-function(min.temp,max.temp,sun.rise,sun.set)
{
  x<-matrix(rep(0:23,NROW(min.temp)),nrow=NROW(min.temp), ncol=24, byrow=TRUE)
  A=(max.temp-min.temp)/2
  day.length<-sun.set-sun.rise
  lengt<-(24-13.58989)+sun.rise
  fr=1/(lengt*1.5)
  phase=sun.rise-0.2044346
  of=A+min.temp
  y<-A*sin(2*pi*fr*(x-phase)-(2*pi/4))+of
  y
}
# generates hourly values
#inputs:
# min.temp = minimum daily temperature
# max.temp = maximum daily temperature
# next.min = minimum daily temperature the next day
# prev.max = maximum daily temperature the previous day
# sun.rise = sunrise time expressed a decimal hour (24hrs) (see functions above)
# sun.set  = sunset time expressed a decimal hour (24hrs)  (see functions above)
# output:
# matrix of 24 values(cols) corresponding to estimated temperature in each hour for each grid cell (row)  
#(first value=midnight, last value= 23:00 hrs)
generate.hrtemps.grid<-function(min.temp,max.temp,next.min,prev.max,sun.rise,sun.set)
{
  x<-matrix(rep(0:23,NROW(min.temp)),nrow=NROW(min.temp), ncol=24, byrow=TRUE)
  day.length<-sun.set-sun.rise
  day1<-generate.hrtemps.grid.day1(min.temp,max.temp,day.length)
  day2<-generate.hrtemps.grid.day2(next.min,max.temp,sun.rise)
  night<-generate.hrtemps.grid.night(min.temp,prev.max,sun.rise,sun.set)
  x<-day1
  x[,15:24]<-day2[,15:24]
  # Replace x when column<= sun.rise and replace with night
  x<-ifelse(col(x)<=sun.rise,night,x) 
  return(x)
}


##########################################################################################
# Main Function
# Writes daily 5km temperature files containing hourly temperature data for whole 5km area
# USES: fill.5km.map and relted FUNCTIONS from elevdif_map_functions
# i.e. records from nearest 5km historic data grid cell
#######################################################################################
hourly_temperatures<-function(start.jd,end.jd,dir_temp,dir_hrtemp,grid5km.r){ 
  plothrs=FALSE
  e.dem<-extent(grid5km.r)  
  jd<-start.jd
  
  # Read data for jd and jd-1
  tmaxfile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
  tminfile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
  prev.tmaxfile<-paste(dir_temp,"MaxTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_Actual.txt", sep="")
  prev.tminfile<-paste(dir_temp,"MinTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_Actual.txt", sep="")
  
  # Check if files exist and assign current day file if no previous or next file
  if(file.exists(tmaxfile)==FALSE|file.exists(tminfile)==FALSE) stop(paste("No 5km data file for jd=",jd,sep=""))
  if(file.exists(prev.tmaxfile)==FALSE) { 
    print("Warning: no prev data file")
    prev.tmaxfile<-tmaxfile 
    prev.tminfile<-tminfile 
  }  
  # Load data from files
  prev.tmax<-raster(prev.tmaxfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
  prev.tmin<-raster(prev.tminfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
  day.tmax<-raster(tmaxfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
  day.tmin<-raster(tminfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
   
  # Resample all daily rasters  using tps method and dem as model 
  day.tmax<-tps.resample(crop(day.tmax,dembuf),dem)
  day.tmin<-tps.resample(crop(day.tmin,dembuf),dem)
  prev.tmax<-tps.resample(crop(prev.tmax,dembuf),dem)
  prev.tmin<-tps.resample(crop(prev.tmin,dembuf),dem)
  
  for (jd in start.jd:end.jd) {  # Loop to calculate hourly files for each day
    ptm<-proc.time(); 
    
    # Define output matrix
    t100m.day<-array(0,dim=c(nrow(dem),ncol(dem),24))
    
    # Define files of next day's temperature data
    next.tmaxfile<-paste(dir_temp,"MaxTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_Actual.txt", sep="")
    next.tminfile<-paste(dir_temp,"MinTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_Actual.txt", sep="")
    
    if(file.exists(next.tminfile)==FALSE|file.exists(next.tmaxfile)==FALSE) {
      print(paste("Warning: no next day's data file for jd= ",jd,sep=""))
      next.tminfile<-tminfile
      next.tmaxfile<-tmaxfile
    }

    next.tmax<-raster(next.tmaxfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    next.tmin<-raster(next.tminfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    next.tmax<-tps.resample(crop(next.tmax,dembuf),dem)
    next.tmin<-tps.resample(crop(next.tmin,dembuf),dem)
    
    compareRaster(prev.tmax,prev.tmin,day.tmax,day.tmin,next.tmax,next.tmin)
    
    # 1. Calculate sunrise, sunset and day length for each grid cell
    # If first day then create lat/lon grid for daylength etc calculations from cropped temp data
    if (jd==start.jd){
      osgrid<-SpatialPoints(coordinates(day.tmax), proj4string=CRS("+init=epsg:27700"), bbox = NULL)
      os.m<-coordinates(osgrid)
      llgrid<-spTransform(osgrid,CRS("+init=epsg:4326"))
      ll.m<-coordinates(llgrid) # creates vector of coordinates from top left  grid cell to bottom right by rows
    }
    # Define vectors to calculate sunrise/set and daylength from lat/lon 
    lat<-ll.m[,"y"]
    long<-ll.m[,"x"]
    numcells<-length(long) # number of grid cells in map = length of vector
    
    sunup<-rep(0,numcells)
    sundown<-rep(0,numcells)
    daylength<-rep(0,numcells)
    jd.v<-rep(jd,numcells)
    #year.v<-rep(DMYjd(jd)$year,numcells)
      
    # CHECK JD calculations !!
    sunup<-sunrise.grid(jd.v,lat,long,rep(0,numcells),rep(0,numcells))
    sundown<-sunset.grid(jd.v,lat,long,rep(0,numcells),rep(0,numcells))
    daylength<-sundown-sunup
    
    ##### 2. Generate hourly temperatures using matrices
    # Create vectors of grid cell daily min/max temperature values
    tmax<-getValues(day.tmax) # from top left to bottom right by row
    tmin<-getValues(day.tmin)
    p.tmax<-getValues(prev.tmax)
    n.tmin<-getValues(next.tmin)
    
    # Create hourly temperatures for each grid cell
    temp.hr<-generate.hrtemps.grid(tmin,tmax,n.tmin,p.tmax,sunup,sundown) 
    
    # Plot hourly values
    for (hr in 0:23) {
      print(paste("Date: ", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep="")," ",hr,":00",sep=""))
      t100m.day[,,hr+1]<-matrix(temp.hr[,hr+1],nrow=NROW(t100m.day),ncol=NCOL(t100m.day),byrow=TRUE) # convert vector to matrix
      # plot raster option
      if (plothrs==TRUE){
        hr.r<-raster(t100m.day[,,hr+1],template=dem)
        temp.plot(hr.r,hr,DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year) }
    }
    
    # Write day of hrly data as 3D matrix
    file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),"_100m.r", sep="") # define file name from year,month,day,hr
    print(paste("Writing file: ",file.out,sep=""))
    save(t100m.day, file=file.out)
    
    # Reset data for next jd
    prev.tmax<-day.tmax
    prev.tmin<-day.tmin
    day.tmax<-next.tmax
    day.tmin<-next.tmin
    remove(next.tmax,next.tmin)
    
    print(proc.time()-ptm) 
    
  } # end for day loop
  
}# end function

##########################################################################################
# Main CODE
##########################################################################################
hourly_temperatures(start.jd,end.jd,dir_temp,dir_hrtemp,grid5km.r)
  