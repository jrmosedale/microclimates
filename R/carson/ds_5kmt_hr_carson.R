##########################################################################################
# TEMPORAL DOWNSCALE OF TEMPERATURE DATA
# OUtput: hourly 5km data in daily files as raster stack
##########################################################################################
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- args[1] 
start.month<-args[2] 
start.year<-args[3] 
end.day<-args[4] 
end.month<-args[5] 
end.year<-args[6] 

## Downscale daily TEMPERATURE data to hourly  Prog: t5km_to_hrmatrix ###
#start.day<-1; start.month<-1; start.year<-1983
#end.day<-31; end.month<-12; end.year<-2013
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)

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
  plothrs=TRUE
  fill<-FALSE
  e<-extent(grid5km.r)  
  # Import Max temperature file available for day before start date (start.jd-1) 
  infile<-paste(dir_temp,"MaxTemp_", DMYjd(start.jd-1)$year, "-",sprintf("%02d",DMYjd(start.jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(start.jd-1)$day,sep=""),"_Actual.txt", sep="")
  tmpdata.r<-raster(infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
  tmpdata.r<-crop(x=tmpdata.r,y=e) # crop to geographical extent of DEM raster
  
  #Create lat/lon grid for daylength etc calculations from cropped temp data
  osgrid<-SpatialPoints(coordinates(tmpdata.r), proj4string=CRS("+init=epsg:27700"), bbox = NULL)
  os.m<-coordinates(osgrid)
  llgrid<-spTransform(osgrid,CRS("+init=epsg:4326"))
  ll.m<-coordinates(llgrid) # creates vector of coordinates from top left  grid cell to bottom right by rows
  
  for (jd in start.jd:end.jd) {  # Loop to calculate hourly files for each day
    # Define output matrix
    t5km.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
    # Read day temperature data
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
    day.tmax<-raster(max.infile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    day.tmin<-raster(min.infile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Read previous day temperature data
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_Actual.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_Actual.txt", sep="")
    prev.tmax<-raster(max.infile, xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    prev.tmin<-raster(min.infile, xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Read NEXT day file 
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_Actual.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_Actual.txt", sep="")
    next.tmax<-raster(max.infile, xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    next.tmin<-raster(min.infile, xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Crop all daily rasters  using dem as model - USE BRICK?? - COMBINE WITH BELOW?
    day.tmax<-crop(x=day.tmax,y=e)
    day.tmin<-crop(x=day.tmin,y=e)
    next.tmax<-crop(x=next.tmax,y=e)
    next.tmin<-crop(x=next.tmin,y=e)
    prev.tmax<-crop(x=prev.tmax,y=e)
    prev.tmin<-crop(x=prev.tmin,y=e)
    
    # OPTION - Fill missing 5km cells without temperature values but containing land cells from nearest neighbour
    if (fill){
        day.tmax<-fill.5km.map(day.tmax,dem,grid5km.r) 
        day.tmin<-fill.5km.map(day.tmin,dem,grid5km.r) 
        next.tmax<-fill.5km.map(next.tmax,dem,grid5km.r) 
        next.tmin<-fill.5km.map(next.tmin,dem,grid5km.r) 
        prev.tmax<-fill.5km.map(prev.tmax,dem,grid5km.r) 
        prev.tmin<-fill.5km.map(prev.tmin,dem,grid5km.r) 
        }
    
    # 1. Calculate sunrise, sunset and day length for each grid cell
    
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
    tmax<-getValues(crop(x=day.tmax,y=e)) # from top left to bottom right by row
    tmin<-getValues(crop(x=day.tmin,y=e))
    p.tmax<-getValues(crop(x=prev.tmax,y=e))
    n.tmin<-getValues(crop(x=next.tmin,y=e))
    
    # Create hourly temperatures for each grid cell
    temp.hr<-generate.hrtemps.grid(tmin,tmax,n.tmin,p.tmax,sunup,sundown) 
    print(paste("Date: ", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),sep=""))
    
    # Write one hourly temperatures as single R matrix cellxhour for each day at 5km resolution 
    # file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),".r", sep="") # define file name from year,month,day,hr
    # save(temp.5km.hr, file=file.out)
    
    # Plot hourly values
    for (hr in 0:23) {
      t5km.day[,,hr+1]<-matrix(temp.hr[,hr+1],nrow=NROW(day.tmax),ncol=NCOL(day.tmax),byrow=TRUE) # convert vector to matrix
      # plot raster option
      if (plothrs==TRUE){
        hr.r<-raster(nrow=nrow(day.tmax),ncol=ncol(day.tmax))
        hr.r<-setValues(hr.r,temp.hr[,hr+1]) # convert matrix to raster
        plot(hr.r,main=paste(DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year,":",hr) )
        }
    }
    
    # Write day of hrly data as 3D matrix
    file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),".r", sep="") # define file name from year,month,day,hr
    save(t5km.day, file=file.out)
    
  } # end for day loop
  
}# end function

##########################################################################################
# Call Function
##########################################################################################

# includes data for previous day to start.day
print(paste("Start: ",start.day,"/",start.month,"/",start.year,sep=""))
print(paste("End:",end.day,"/",end.month,"/",end.year,sep=""))
print (paste("Start JD= ",start.jd,sep=""))
print(paste("End JD ",end.jd,sep=""))

hourly_temperatures(start.jd,end.jd,dir_temp,dir_hrtemp,grid5km.r) 
