library(ncdf4)
library(raster)
library(rgdal)
library(RAtmosphere)
library(insol)
library(sp)

# INPUTS to correct:
#       Longitude and Lattitude in sunrise/set Function calls (and check timezone)
#       NOBS in daily weather data in hourly generation loop code (currently set as 333096)
# TIME to run
#       ~5-6 mins for one year of data - about 1GB of files written per 10 years 

##########################################################################################
### FUNCTIONS
### Require: change from DOY to JD or d/m/y??
##########################################################################################
# Plots tmp map rasters to fixed scale
temp.plot<-function(temp,hr,day,month,year) {
t<-paste(day,"-",month,"-",year," ",hr,":00",sep="")
brk<-c(-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)
col<-rev(heat.colors(19))
plot(temp,main=t,col=col,breaks=brk)
}

# Calculates DOY (from 1980-2015) from day, month, year
CalcDOY<-function(day,month,year){
  y=year-1979 # so 1=1980
  feb.d<-c(29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28,29,28,28,28,
           29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28) # days of Feb from 1980 to 2015
  monthdays<-c(31,feb.d[y],31,30,31,30,31,31,30,31,30,31)
  doy<-sum(monthdays[1:month-1])+day
  return(doy)
}

# Function to convert JD to day, month year as list
# Gives Gregorian dates (even for julian period before 15C) 
# Inputs and outputs are vectors
# Based on: https://en.wikipedia.org/wiki/Julian_day 
# Check if correction to Year (-1 when month is Mar-Dec)
# Check with: http://aa.usno.navy.mil/faq/docs/JD_Formula.php
DMYjd<-function(JD) { #assumes 12:00 for JD
  y<-4716; v<-3
  j<-1401; u<-5
  m<-2; s<-153
  n<-12; w<-2
  r<-4; B<-274277
  p<-1461; C<--38
  f<-JD+j+(((4*JD+B)%/%146097)*3)%/%4+C
  e<-r*f+v #
  g<-(e%%p)%/%r
  h<-u*g+w
  D <- (h%%s) %/% u+1
  M<- ( (h%/%s + m) %% n)+1
  Y <- ifelse(M<3,(e%/%p) - y + (n+m-M) %/%n, (e%/%p) - (y) + (n+m-M) %/%n) 
  dmy<-list(day=D,month=M,year=Y)
  return(dmy)
}

# Calculate JD from DMY - input/output as vectors
# Based on http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
# confirmed using https://www.aavso.org/jd-calculator  and  http://aa.usno.navy.mil/data/docs/JulianDate.php
JDdmy<-function(day,month,year){ # JD for 00:00 on date given (ie JD.5) 
  year<-ifelse(month<3,year-1,year)
  month<-ifelse(month<3,month+12,month)
  A <-floor(year/100);#print(A)
  B <-floor(A/4);#print(B)
  C<- 2-A+B; #print(C)
  E <- floor(365.25*(year+4716)); #print(E)
  F <- floor(30.6001*(month+1)); #print(F)
  JD<- C+day+E+F-1524.5; #print (JD)
  return(JD)
}

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
  sun.rise<-suntimes.grid(JD,Lat,Long,Timezone,DST)[,1]
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
# vector of 24 values corresponding to estimated temperature in each hour 
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
  x
}

##########################################################################################
# Define directories
#dir_temp<-"C:/Data2015/Temp5km/"
#dir_hrtemp<-"C:/Data2015/Temp5km/hourly/"

dir_zip<-"~/Documents/Exeter/Data2015/Temp5km/zip/"
dir_temp<-"~/Documents/Exeter/Data2015/Temp5km/extract/"
dir_hrtemp<-"~/Documents/Exeter/Data2015/Temp5km/hourly/" # dir for output files

#######################################################################################
#### Read in Digital Eelevation data - DECISION: what area we want to cover ####
#######################################################################################
#dem<-raster("C:/Data2015/DEM100/dem_sw_x60-420k_y-10-180k.tif")
dem<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif")
plot(dem,main="DEM-full")
extent(dem)
#e.dem<-extent(c(70000,420000,0,180000)) # includes scilly isles
e.dem<-extent(c(120000,420000,0,180000)) # excludes scilly isles
dem<-crop(dem,e.dem)
plot(dem,main="DEM-sw")
e.dem <-extent(dem)

#######################################################################################

# Define DATES for which to create hourly files and caulculate JDs
# Designed to analyse by YEAR of data
start.d<-1; start.m<-1; start.y<-2010
end.d<-31; end.m<-1; end.y<-2010
print (paste("Start: ",start.d,"/",start.m,"/",start.y,"  End: ",end.d,"/",end.m,"/",end.y,sep=""))

start.jd<-JDdmy(start.d,start.m,start.y) 
end.jd<-JDdmy(end.d,end.m,end.y)
print (paste("Start JD: ",start.jd,"  End JD: ",end.jd,sep=""))

#######################################################################################

# Unzip year files of temperature data here
zip.file<-paste(dir_zip,"MinTemp_", DMYjd(start.jd)$year,".zip", sep="")
print (paste("Unzipping ",zip.file,sep=""))
unzip(zip.file, exdir=dir_temp)
zip.file<-paste(dir_zip,"MaxTemp_", DMYjd(start.jd)$year,".zip", sep="")
print (paste("Unzipping ",zip.file,sep=""))
unzip(zip.file, exdir=dir_temp)


# Check if start or end dates are 1st or last of a year and unzips data for these years as well
if ( start.d==1&start.m==1) { # then extract previous year of data
  zip.file<-paste(dir_zip,"MinTemp_", DMYjd(start.jd)$year-1,".zip", sep="")
  print (paste("Unzipping ",zip.file,sep=""))
  unzip(zip.file, exdir=dir_temp) 
  zip.file<-paste(dir_zip,"MaxTemp_", DMYjd(start.jd)$year-1,".zip", sep="")
  print (paste("Unzipping ",zip.file,sep=""))
  unzip(zip.file, exdir=dir_temp)
}
if ( end.d==31&end.m==12) { # then extract next year of data
  zip.file<-paste(dir_zip,"MinTemp_", DMYjd(start.jd)$year+1,".zip", sep="")
  print (paste("Unzipping ",zip.file,sep=""))
  unzip(zip.file, exdir=dir_temp) 
  zip.file<-paste(dir_zip,"MaxTemp_", DMYjd(start.jd)$year+1,".zip", sep="")
  print (paste("Unzipping ",zip.file,sep=""))
  unzip(zip.file, exdir=dir_temp)
}
#######################################################################################

# Import temperature file available for day before start date (start.jd-1) 
infile<-paste(dir_temp,"MaxTemp_", DMYjd(start.jd-1)$year, "-",sprintf("%02d",DMYjd(start.jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(start.jd-1)$day,sep=""),"_ACTUAL.txt", sep="")
tmpdata.r<-raster(infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
tmpdata.r<-crop(x=tmpdata.r,y=dem) # crop to geographical extent of DEM raster

#Create lat/lon grid for daylength etc calculations from cropped temp data
osgrid<-SpatialPoints(coordinates(tmpdata.r), proj4string=CRS("+init=epsg:27700"), bbox = NULL)
os.m<-coordinates(osgrid)
llgrid<-spTransform(osgrid,CRS("+init=epsg:4326"))
ll.m<-coordinates(llgrid) # creates vector of coordinates from top left  grid cell to bottom right by rows

# Read start day temperature data
max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(start.jd)$year, "-",sprintf("%02d",DMYjd(start.jd)$month,sep=""),"-", sprintf("%02d",DMYjd(start.jd)$day,sep=""),"_ACTUAL.txt", sep="")
min.infile<-paste(dir_temp,"MinTemp_", DMYjd(start.jd)$year, "-",sprintf("%02d",DMYjd(start.jd)$month,sep=""),"-", sprintf("%02d",DMYjd(start.jd)$day,sep=""),"_ACTUAL.txt", sep="")
day.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
day.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")

# Read previous day temperature data
max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(start.jd-1)$year, "-",sprintf("%02d",DMYjd(start.jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(start.jd-1)$day-1,sep=""),"_ACTUAL.txt", sep="")
min.infile<-paste(dir_temp,"MinTemp_", DMYjd(start.jd-1)$year, "-",sprintf("%02d",DMYjd(start.jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(start.jd-1)$day-1,sep=""),"_ACTUAL.txt", sep="")
prev.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
prev.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")


for (jd in start.jd:end.jd) {  # Loop to calculate hourly files for each day
  
      # Read NEXT day file 
      max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="")
      min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="")
      next.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
      next.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
      
      # Crop all rasters individually using dem as model
      day.tmax<-crop(x=day.tmax,y=dem)
      day.tmin<-crop(x=day.tmin,y=dem)
      next.tmax<-crop(x=next.tmax,y=dem)
      next.tmin<-crop(x=next.tmin,y=dem)
      prev.tmax<-crop(x=prev.tmax,y=dem)
      prev.tmin<-crop(x=prev.tmin,y=dem)
      
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
      daylength=sundown-sunup
      
      ##### 2. Generate hourly temperatures using matrices
      # Create matrices of daily min/max temperatures
      tmax<-getValues(crop(x=day.tmax,y=dem)) # from top left to bottom right by row
      tmin<-getValues(crop(x=day.tmin,y=dem))
      p.tmax<-getValues(crop(x=prev.tmax,y=dem))
      n.tmin<-getValues(crop(x=next.tmin,y=dem))
      
      # Create hourly temperatures for each grid cell
      hr.temps<-generate.hrtemps.grid(tmin,tmax,n.tmin,p.tmax,sunup,sundown) 
      print(paste("Date: ", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),sep=""))
      
      # Write 24 hourly temp files at 5km resolution (as downscaling will depend upon other variables)
      for (hr in 0:23) {
            # convert vector to matrix and raster and plot (use osgrid as template)
            hr.m<-matrix(hr.temps[,hr+1],nrow=NROW(day.tmax),ncol=NCOL(day.tmax),byrow=TRUE) # convert array to matrix
            hr.r<-raster(hr.m, template=day.tmax) # convert matrix to raster

            # Hourly plots
            # temp.plot(hr.r,hr,DMYjd(jd)$day[1],DMYjd(jd)$month[1],DMYjd(jd)$year[1])

            # Write matrix as R file
            data.out<-hr.m
            file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),"-",sprintf("%02d",hr,sep=""),"00.r", sep="") # define file name from year,month,day,hr
            write(data.out, file=file.out)
      } # end for hr 
      
      # Set prev to equal day and day to equal next 
      prev.tmax<-day.tmax
      prev.tmin<-day.tmin
      day.tmax<-next.tmax
      day.tmin<-next.tmin

} # end for day loop

#######################################################################################
# Delete unzipped files used
for (jd in start.jd:end.jd) {
  max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="")
  min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="") 
  file.remove(max.infile)
  file.remove(min.infile)
}

#######################################################################################


# CUT OUTS
# Create stack and crop to extent of dem
#daily.stk<-stack(day.tmax,day.tmin, prev.tmax,prev.tmin,next.tmax,next.tmin)
#names(daily.stk)<-c("day.tmax","day.tmin","prev.tmax", "prev.tmin","next.tmax","next.tmin")
#daily.stk<-crop(x=daily.stk,y=dem)

# create matrix  of 24 x numcells (i.e. temp map for each hr) 
# to test set i to vector range
# parameters<-paste("i=",i," tmin=",tmin[i]," tmax=",tmax[i]," n.tmin=",n.tmin[i]," p.tmax=",p.tmax[i]," sunup=", sunup[i]," sunset=", sundown[i],sep="")
# print(parameters)
# hr.temps<-generate.hrtemps.grid(tmin[i],tmax[i],n.tmin[i],p.tmax[i],sunup[i],sundown[i]) 

# resample to 100m - do this before or after creation of hourly temps? 
# include elevation correction based on 100m cell difference from 5km mean
# ?? coastal effects etc??
