
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
hourly_temperatures<-function(jd,dir_temp,dir_hrtemp,grid5km.r, plothrs=FALSE){ 
  fill<-FALSE
  e<-extent(grid5km.r)  

  #Create lat/lon grid for daylength etc calculations from cropped temp data
  osgrid<-SpatialPoints(coordinates(grid5km.r), proj4string=CRS("+init=epsg:27700"), bbox = NULL)
  os.m<-coordinates(osgrid)
  llgrid<-spTransform(osgrid,CRS("+init=epsg:4326"))
  ll.m<-coordinates(llgrid) # creates vector of coordinates from top left  grid cell to bottom right by rows
  
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
    
    return(t5km.day)
    
    # Write day of hrly data as 3D matrix
    # file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),".r", sep="") # define file name from year,month,day,hr
    #vsave(t5km.day, file=file.out)
    
}# end function


#######################################################################################
# FUNCTIONS TO CALCULATE RELATIVE HUMIDITY AT 5KM RESOLUTION
#######################################################################################
#dir_rh5km<-"~/Documents/Exeter/Data2015/RelHumidity/rh5km/"
add.zero<-function(x)
{
  y<-x
  if (y<9) y<-paste("0",x,sep="")
  y
}

rel.to.abs<-function(rh,t)
{
  e0.1<-0.6108*exp(17.27*t/(t+237.3))
  e<-e0.1*(rh/100)
  abso<-(2165*e)/(t+273.16) # grams per metre cubed
  abso
}

abs.to.rel<-function(ab,t)
{
  e<-ab*(t+273.16)/2165
  e0.1<-0.6108*exp(17.27*t/(t+237.3))
  rh<-100*(e/e0.1)
  rh
}

# Centres raster on GML 
reslice.raster<-function(r)
{
  arr<-array(0,dim=c(73,144))
  e1<-extent(1.25,181.25,-91.25,91.25)
  e2<-extent(181.25,358.75,-91.25,91.25)
  e3<-extent(-1.25,1.25,-91.25,91.25)
  arr[,1:71]<-getValues(crop(r,e2),format="matrix")
  arr[,72]<-getValues(crop(r,e3),format="matrix")
  arr[,73:144]<-getValues(crop(r,e1),format="matrix")
  r2<-raster(arr,xmn=-178.25,xmx=181.25,ymn=-91.25,ymx=91.25)
  res(r2)<-2.5
  r2
}

reproject.r<-function(a,e,grid5km.r)
{
  r<-raster(a,xmn=xmin(e),xmx=xmax(e),ymn=ymin(e),ymx=ymax(e))
  projection(r)<-"+init=epsg:4326"
  r2<-projectRaster(r,crs="+init=epsg:27700")
  r3<-resample(r2,grid5km.r)
  #r4<-crop(r3,e2)
  r3
}

rh.plot<-function(rh,hr,day,month,year) 
{
  t<-paste("RH on ",day,"-",month,"-",year," ",hr,":00",sep="")
  brk<-c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,160,170,180,190,200)
  col<-rev(heat.colors(21))
  plot(rh,main=t,col=col,breaks=brk)
}
###########################################################################
# Main Function: Write daily files of hourly 5km Rel Humidity 
# RH data in = 4xdaily global data
###########################################################################
rh.hourly<-function(jd,dir_rh,dir_rh5km,hr.temp1,hr.temp2, grid5km.r,hourplot=FALSE)
{
  e<-extent(-5.75,-0.75,48.75,53.75)  # lon/lat extent of cropped rasters
  #e<-extent(353.75,359.25,48.75,53.75)
  
  # Define output matrix to hold one day of hourly relative humidity values at 5km resolution for area covered by grid5km.r
  rh.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
  
  # Reads in Relative Humidity data from yearly file if new year
  infile.rh<-paste(dir_rh,"rhum.sig995.",DMYjd(jd)$year,".nc",sep="") # reads yearly data file
  print(paste("This year's data file= ",infile.rh,sep=""))
  rh.s<-stack(infile.rh) # all times
  # Add next day's data from next year file if necessary
  if ((DMYjd(jd)$day==31 & DMYjd(jd)$month==12)){
    nextfile.rh<-paste(dir_rh,"rhum.sig995.",DMYjd(jd)$year+1,".nc",sep="") # reads yearly data file
    print(paste("Adding next year's data file= ",nextfile.rh,sep=""))
    nextrh.s<-stack(nextfile.rh)
    rh.s<-addLayer(rh.s,nextrh.s[[1]]) # Adds first layer of next year to currnet year
    remove(nextrh.s)
    # get times
    netRH<-nc_open(infile.rh)
    tm<-ncvar_get(netRH,"time")
    tm<-c(tm,tm[length(tm)]+6) # =adds time of first reading of next year
   }
    
  for (period in 0:3)
  { print(paste("Period= ",period,sep=""))
    hr<-(period*6)
    # Identify rh layers and crop
    Jul.base<-JDdmy(1,1,1800)
    #Jul.actual<-JD(ISOdate(year,month,day))
    hr.val<-(jd-Jul.base)*24+hr
    sel<-which(tm==hr.val)
    rh.all0<-subset(rh.s,sel)
    rh.all6<-subset(rh.s,(sel+1))
    rh.all0<-reslice.raster(rh.all0)
    rh.all6<-reslice.raster(rh.all6)
    projection(rh.all0)<-"+init=epsg:4326"
    projection(rh.all6)<-"+init=epsg:4326"
    rh0<-crop(rh.all0,e)
    rh6<-crop(rh.all6,e) 
    #plot(rh0,main="rh0")
    #plot(rh6,main="rh6")
    
    # read in temperature rasters
    t0.5km<-raster(hr.temp1[,,hr+1],template=grid5km.r)
    if (period<3) t6.5km<-raster(hr.temp1[,,hr+7],template=grid5km.r)
    if (period==3) t6.5km<-raster(hr.temp2[,,1],template=grid5km.r) # 0hr00 for following day
    
    # reproject and calculate aggregate temperatures to match humidity data resolution
    projection(t0.5km)<-"+init=epsg:27700"
    projection(t6.5km)<-"+init=epsg:27700"
    t0.ll<-projectRaster(t0.5km,crs="+init=epsg:4326")
    t6.ll<-projectRaster(t6.5km,crs="+init=epsg:4326")
    t0.ll<-raster::extend(t0.ll,e) # required to allow aggregate values of rh cells only partly covered by temp data
    t6.ll<-raster::extend(t6.ll,e)
    tem<-raster(array(0,dim=c(2,2)),xmn=xmin(e),xmx=xmax(e),ymn=ymin(e),ymx=ymax(e))
    t0<-raster::resample(t0.ll,rh6)
    t6<-raster::resample(t6.ll,rh6)
    
    #########################################
    
    # Convert to absolute humidity
    abs.hum0<-rel.to.abs(getValues(rh0,format="matrix"),getValues(t0,format="matrix"))
    abs.hum6<-rel.to.abs(getValues(rh6,format="matrix"),getValues(t6,format="matrix"))
    
    # interpolate for missing hours
    abs.hum1<-(abs.hum0*5+abs.hum6*1)/6
    abs.hum2<-(abs.hum0*4+abs.hum6*2)/6
    abs.hum3<-(abs.hum0*3+abs.hum6*3)/6
    abs.hum4<-(abs.hum0*2+abs.hum6*4)/6
    abs.hum5<-(abs.hum0*1+abs.hum6*5)/6
    
    # Convert absolute humidity to 5km grid cell resolution 
    e2<-raster::extent(grid5km.r)
    abs.h1<-reproject.r(abs.hum1,e,grid5km.r)
    abs.h2<-reproject.r(abs.hum2,e,grid5km.r)
    abs.h3<-reproject.r(abs.hum3,e,grid5km.r)
    abs.h4<-reproject.r(abs.hum4,e,grid5km.r)
    abs.h5<-reproject.r(abs.hum5,e,grid5km.r)
    abs.h6<-reproject.r(abs.hum6,e,grid5km.r)
    
    # Load remaining hrs of 5km temperature data 
    t1.5km<-raster(hr.temp1[,,hr+2],template=grid5km.r)
    t2.5km<-raster(hr.temp1[,,hr+3],template=grid5km.r)
    t3.5km<-raster(hr.temp1[,,hr+4],template=grid5km.r)
    t4.5km<-raster(hr.temp1[,,hr+5],template=grid5km.r)
    t5.5km<-raster(hr.temp1[,,hr+6],template=grid5km.r)
    
    # Convert to 5km relative humidity
    hr<-period*6
    rh.day[,,hr+1]<-abs.to.rel(getValues(abs.h1,format="matrix"),getValues(t1.5km,format="matrix"))
    rh.day[,,hr+2]<-abs.to.rel(getValues(abs.h2,format="matrix"),getValues(t2.5km,format="matrix"))
    rh.day[,,hr+3]<-abs.to.rel(getValues(abs.h3,format="matrix"),getValues(t3.5km,format="matrix"))
    rh.day[,,hr+4]<-abs.to.rel(getValues(abs.h4,format="matrix"),getValues(t4.5km,format="matrix"))
    rh.day[,,hr+5]<-abs.to.rel(getValues(abs.h5,format="matrix"),getValues(t5.5km,format="matrix"))
    rh.day[,,hr+6]<-abs.to.rel(getValues(abs.h6,format="matrix"),getValues(t6.5km,format="matrix"))
    
    
  } # end period
  if (hourplot==TRUE){
    for (t in 1:24) {
      rh.plot(raster(rh.day[,,t],template=grid5km.r),t-1,DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year) 
    } }
  
  return(rh.day)
  # WRITE DAILY files of 5km RELATIVE humidity
  # fileout<-paste(dir_rh5km,"RH_5km_",DMYjd(jd)$year,"_",DMYjd(jd)$month,"_",DMYjd(jd)$day,".R",sep="")
  #print(fileout)
  # save(rh.day,file=fileout)

} # end function



####################################################
# Downscale SST
# Input: single ncdf file of Met Office Hadley Centre HadISST1 data defined in sst.spdownsc function
# Output: Spatially and temporally ownscaled SST eg daily at 5km
# Issues/Problems: interpolation/resampling - in effect across N/S coasts
# Main function is sst.time.int - extracts monthly data and downscales
####################################################

# Functions to return year/month from sst t value
ttoyr<-function(t){
  yr<-ceiling(t/12)-1+1960
  return(yr)
}

ttomonth<-function(t){
  month<-t-((ttoyr(t)-1960)*12)
  return(month)
}

# Function to spatially downscale sst data for month t to dsgrid  
sst.spdownsc<-function(year,month,dsgrid.r,dir_sst,dir_sstm) {
  t<-(year-1960)*12+month
  
  in.file<-paste(dir_sst,"HadISST_sst.nc",sep="") # in.file required for sst.spdownsc
  ncfile<-nc_open(paste(dir_sst,"HadISST_sst.nc",sep=""))
  
  # Read ncdf direct to raster
  sst.world.r<-raster(in.file,band=t) # read all of level t of ncdf file
  #plot(sst.world.r)
  # Extract relevant 1 deg cells 353-360E, 48:53N, 1960-2015 
  long.mn<--7;long.mx<--1 # -7 to -1
  lat.mn<-49; lat.mx<-52 # 49 to 52
  e.sst<-extent(c(long.mn,long.mx,lat.mn,lat.mx))
  sst.r<-crop(sst.world.r,e.sst)
  #plot(sst.r,main="t")
  #plot(sst.r,main=paste("SST long: ",long.mn," to ",long.mx,", lat: ",lat.mn," to ",lat.mx,sep=""))
  
  # Set land cells to have same temperature as nearby sea
  #sst.r<-focal(sst.r,w=matrix(1,3,3),fun=mean,na.rm=TRUE,NAonly=TRUE)
  sst.m<-getValues(sst.r,format="matrix")
  sst.m[1,5]<-(sst.m[1,4]) # Bristol channel - count as sea = to adjacent cell (NOT mean including CHannel!)
  #sst.m[1,6]<-(sst.m[2,6]+sst.m[1,5])/2 # keep as NA all land cells
  sst.r<-raster(sst.m,template=sst.r)
  #plot(sst.r)
  
  # Convert projection to OSGB
  projection(sst.r)<-"+init=epsg:4326"
  sst.os<-projectRaster(sst.r,crs="+init=epsg:27700")
  
  # resample to 5km resolution - cf with ngb method
  sst5km.r<-raster::resample(sst.os,dsgrid.r,method="bilinear",na.rm=TRUE)
  # plot(sst5km.r)
  
  nc_close(ncfile)
  
  return(sst5km.r)
}# end of function

##########################################################################################
# 2. TEMPORAL DOWNSCALING FUNCTIONS to DAILY  time period (could also be used for hourly period if necessary)
##########################################################################################

nday.in.month <- function(date)
{
  m<-format(date, format="%m")
  while(format(date, format="%m") == m) date <- date + 1
  return(as.integer(format(date - 1, format="%d")))
}

nday.in.month.jd <- function(jd)
{
  m<-DMYjd(jd)$month
  jd1<-jd-DMYjd(jd)$day+1 # jd for 1st of the month
  while(DMYjd(jd1)$month == m) jd1 <- jd1 + 1
  return(DMYjd(jd1-1)$day)
}

# Uses julian date functions and function nday.in.month.jd
sst.time.int<-function(jd,dsgrid.r, dir_sstm,dir_ssth,hr=12,msk=FALSE)
{
  # convert julian date to day/month/year
  day<-DMYjd(jd)$day ; month<-DMYjd(jd)$month ; year<-DMYjd(jd)$year
  
  # How many days in month
  dimth<-nday.in.month.jd(jd)
  
  # Calculates SSTs for the two months that lie either side of the date in question
  # reads in month before if date is in first half of month - using jd allows for month being in different year
  mth1<-DMYjd(jd-(dimth/2))$month  
  yr1<-DMYjd(jd-(dimth/2))$year
  mth2<-DMYjd(jd+(dimth/2))$month 
  yr2<-DMYjd(jd-(dimth/2))$year
  
  r1<-sst.spdownsc(yr1,mth1,dsgrid.r,dir_sst,dir_sstm) 
  r2<-sst.spdownsc(yr2,mth2,dsgrid.r,dir_sst,dir_sstm) 
  
  v1<-getValues(r1,format="matrix")
  v2<-getValues(r2,format="matrix")
  
  # work out weighting to attach
  
  # How many days in both months
  dimth1<-nday.in.month.jd(jd-(dimth/2))
  dimth2<-nday.in.month.jd(jd+(dimth/2))
  print(paste(day,"/",month, ".  dimth1= ",dimth1," dimth2=", dimth2,sep=""))
  
  # time after dimth1
  tt<-day+hr/24+dimth1/2
  if (day+hr/24>dimth/2) tt<-day+hr/24-dimth1/2
  wgt<-tt/((dimth1+dimth2)/2)
  v<-v2*wgt+v1*(1-wgt)
  sstds.r<-raster(v,template=r1)
  
  if (msk) {
    seagrid.r<-calc(dsgrid.r,function(x) ifelse(is.na(x),0,NA))
    sstds.r<-mask(sstds.r,seagrid.r)
  }
  #ssthr.r<-mask(sst.hr.r,dsgrid.r) # to set land cells to NA
  
  # write daily file
  tl<-paste("Year: ",year," Month: ",month," Day: ",day," Hour: ",hr,sep="")
  #plot(sstds.r,main=tl)
  
  return (sstds.r)
  
} # end function



#######################################################################################
# FUNCTIONS TO CALCULATE LWR from CAL, RH and T
#######################################################################################
# CALCULATE LW RADIATION from CAL , Temp, RH
# Units of output = MJ/m2
# Inputs: rh cropped to block, tref requires cropping, CAL - requires downscaling and cropping
# Requires jd functions

# FUNCTION - calculate long wave radiation from Temp, RH & CAL
lwr<-function(Temp,RH,CAL)
{
  e0<-0.6108*exp(17.27*Temp/(Temp+237.3)) # saturated vapour pressure
  ea<-e0*(RH/100) # actual vapour pressure
  moisture.absorb<-0.34-0.14*sqrt(ea)
  cloud.absorb<-1.35*(1-CAL)-0.35
  rnl<-2.043*10^-10*(Temp+273.16)^4*moisture.absorb*cloud.absorb
  rnl
}

calc_lwr_block<-function(jd,cal.block,rhref.block,tref.block,plotlwr=FALSE,writefile=FALSE)
{
  # Define output file - one day of hourly 100m data
  #proc.time()->ptm
  year<-DMYjd(jd)$year
  month<-DMYjd(jd)$month
  day<-DMYjd(jd)$day
  #print(paste("Date: ",day,"/",month,"/",year,sep=""))
  compareRaster(cal.block,rhref.block,tref.block)
  
  # Calculate longwave radiation
  lwr.m<-lwr(getValues(tref.block,format="matrix"),
             getValues(rhref.block,format="matrix"),
             getValues(cal.block,format="matrix"))
  lwr.block<-raster(lwr.m,template=tref.block)
  
  if (plotlwr==TRUE) {
    par(mfrow=c(2,2))
    dayhr<-paste(day,"/",month,"/",year," ",hr+1,"h00",sep="")
    plot(cal.block,main=paste("Effective cloud albedo ",dayhr,sep=""))
    plot(rhref.block,main=paste("Relative humdidity ",dayhr,sep=""))
    plot(tref.block,main=paste("Temperature ",dayhr,sep=""))
    plot(lwr.block,main=paste("Long-wave radiation ",dayhr,sep=""))
  }
  #print(proc.time()-ptm)
  if (writefile==TRUE){
    file.out<-paste(dir_lwr,"lwr_",year,"_",month,"_",day,"_100m.R",sep="")
    print(paste("File out: ",file.out,sep=""))
    save(lwr.day,file=file.out)
    print(proc.time()-ptm)
  } # if writefile
  return(lwr.block)
} # end function

#######################################################################################
# LATENTHEAT.BLOCK
# This function calculates crop reference evapotranspiration (using the Penman-Monteith equation).
# Inputs and outputs are as rasters excpet for dn which may be raster or a single value
# Details of algorithm here: http://www.fao.org/docrep/x0490e/x0490e00.htm
# Input variables:
# Temp is the temperature at the site in degrees C. You will need to use the anomoly between the
# 5 km grid and the 100m cell in the previous time-step to estimate this (as the CRE function is
# used to derive the local temperature)
# Net radiation  - importantly this is in MJ m-2 hour-1 and may require conversion of units from the
# satellite derived estimates (typical value ~0.2)
# RH - relative humidity expressed as a percentage (typical value 80%)
# P - Atmposheric rressure - in millibars (typical value ~ 1000)
# currently not assumed to vary by location within a 5km grid cell, but we could do an altitude correction:
# http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric pressure (p)
# dn  # a binary variable specifying whether it is day (1)  or   night (0) needed as the equation for Soil Heat Flux
# changes depending on whether it is night or day - can be a Raster or single value
# u2 is wind speed at at 2 m height [m s-1] (typical value ~5)
# Output variable:
#Crop reference evapotranspiration [mm m-2 hr-1]
CRE<-function(Temp,Rn,RH,P,dn,u2){
  e0<-0.6108*exp(17.27*Temp/(Temp+237.3)) # saturated vapour pressure
  ea<-e0*(RH/100) # actual vapour pressure
  delta<-4098*(0.6108*exp(17.27*Temp/(Temp+237.3)))/((Temp+237.3)^2)  # slope vapour pressure curve
  if(class(dn)=="RasterLayer") G<-overlay(Rn,dn,fun = function(x, y) ifelse(y>0, 0.1*Rn, 0.5*Rn)) else {
    if (dn>0) G<-0.1*Rn else G<-0.5*Rn } # Soil heat flux
  gamma<-0.000665*(P/10)  # psychrometric constant
  ET0<-(0.408*delta*(Rn-G)+gamma*37/(Temp+273)*u2*(e0-ea))/(delta+gamma*(1+0.34*u2))
  ET0<-calc(ET0,fun=function(x){ifelse(x>0,x,0)})
  ET0
}

# This function calculates the change in relative humidity between two locations as a result of the temperature change
# Input variables:
# t1 is the reference temperature in degrees C - i.e. the value for each 5 km grid cell
# t2 is the is the temperature at the site in degrees C. You will need to use the anomoly between the
#5 km grid and the 100m cell in the previous time-step to estimate this
# rh is the reference relative humdity, expressed as a percentage - i.e. the value for each 5 km grid cell
# Output variable:
# the relative humdity at the site, expressed as a percentage - i.e. the value for each 100m grid cell
rh.change<-function(t1,t2,rh){
  # absolute humidity at t1
  e0.1<-0.6108*exp(17.27*t1/(t1+237.3))
  e<-e0.1*(rh/100)
  a.hum.1<-(2165*e)/(t1+273.16) # grams per metre cubed
  # rel humidity at t2
  s.e<-(a.hum.1*(t2+273.16))/2165
  s.e0<-0.6108*exp(17.27*t2/(t2+237.3))
  rhs<-(s.e/s.e0)*100
  rhs
}

# This function calculates the amount of water that can be expected to condense,
# as either a result of a change in temperature from one place to anotehr or through time
# Input variables:all as RASTERS
# t1 is the reference temperature in degrees C - i.e. the value for each 5 km grid cell
# t2 is the is the temperature at the site in degrees C. You will need to use the anomoly between the
# 5km grid and the 100m cell in the previous time-step to estimate this
# rh is the reference relative humdity, expressed as a percentage - i.e. the value for each 5km grid cell
# Output variable:
# the amount of water condensed [mm m-2 hr-1]. Zero if relative humidity is less than 100% - as RASTER
Water.conden<-function(t1,t2,rh){
  # absolute humidity at t2, rh=100
  e0.100<-0.6108*exp(17.27*t2/(t2+237.3))
  e100<-e0.100*(100/100)
  a.hum.100<-(2165*e100)/(t2+273.16) # grams per metre cubed
  # absolute humidity at rh2
  rh2<-rh.change(t1,t2,rh)
  e0.2<-0.6108*exp(17.27*t2/(t2+237.3))
  e2<-e0.2*(rh2/100)
  a.hum.2<-(2165*e2)/(t2+273.16) # grams per metre cubed
  a.hum.2<-a.hum.2-a.hum.100
  #a.hum.2<-ifelse(rh2>100,a.hum.2,0)
  a.hum.2<-overlay(a.hum.2,rh2,fun=function(x,y) ifelse(y>100,x,0) )
  a.hum.2
}
##########################################################################################
# PREPARE_PRESSURE
# Prepare Sea level Pressure data - extract data at 100m resolution for time jd
# Requires JD functions
# Works for whole area if block=dembuf or similar or for small block
# Write raster for area grid

downscale.pressure<-function(dem.block,p.ncfile,jd,write.file=FALSE)
{
  # unzip .gz file if necessary
  #gzfile<-paste(dir_pressure025,"pp_0.25deg_reg_v11.0.nc.gz",sep="")
  #ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")
  #gunzip(filename=gzfile, destname=ncfile, overwrite=TRUE)
  #ncdf_test<-nc_open(ncfile)
  
  # Select time bands 
  # Base jd value = 1/1/1950
  jd.base<-JDdmy(1,1,1950)
  jd.band<-jd-jd.base
  
  # Load raster from ncdf file 
  pressure.r<-raster(p.ncfile,band=jd.band)
  #plot(pressure.r)
  
  # Reproject cropped version to OSGB projection 
  pressure.r<-projectRaster(crop(pressure.r,c(-7,0,49,52)),crs="+init=epsg:27700")
  #plot(pressure.r)
  
  # Resample to 100m for block
  p.block<-tps.resample(pressure.r,dem.block)
  #p5km.r<-resample(pressure.r,land5km.r)
  #p5km.r<-fill.5km.map(p5km.r,land5km.r)
  #p.block<-ref5km.to.block100m(dem.block,p5km.r)

  # Crop and write pressure raster (all times) for area of interest
  if(write.file){
    plot(p.block,main=paste("Sea level Pressure 100m cell for ",DMYjd(jd)$day,"/",DMYjd(jd)$month, "/",DMYjd(jd)$year,sep=""))
    outfile<-paste(dir_pressure,"pressure.tif",sep="")  
    writeRaster(file=outfile,p.block,overwrite=TRUE )
  }
  return(p.block)
} # end function

# Correct pressure for elevation according to hypsometric formula
# Inputs: sea level pressure, tmp at location, elevation above sea level
# References:  http://keisan.casio.com/exec/system/1224579725, http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric%20pressure
correct.pressure<-function(sl.pressure,tmp,elevation)
{
  kelvins<-273+tmp
  pressure<-sl.pressure*((kelvins-(0.0065*elevation)) / kelvins )^5.25588
  return(pressure)
}
##########################################################################################
# Calculates single raster  ldif.block
# Input: ldif.stack of ldif for each wind direction
#        wind direction for block
#        direction interval for which ldif calculated

calc.ldif.block<-function(ldif.stack,wdir.block,interval=10){
  ldif.block<-raster
  wdir.layer<-round(wdir.block/interval) # creates raster holding layer in ldif.stack for wind direction
  wdir.layer<-calc(wdir.layer,fun=function(x){ifelse(x==0,36,x)}) # corrects layer 0 to 36 (angle 0 to 360)
  ldif.block<-stackSelect(ldif.stack,wdir.layer)
  return(ldif.block)
} # end function
##########################################################################################

# Functions to calculate Upwind SST and SST-Tref for each 100m cell in block given wind direction and sst map
# Input:  100m dem grid -  block and buffered regions
#         raster of sst for buffer attime  t
#         raster of wind direction for block at time t
# Output: raster of SST and SST-Tref at t

#####################################################################
# FUNCTIONS 
# Used by upwind.sst
findsst<-function(x)
{
  is.sea<-ifelse(x==-999,NA,1)
  nearest.sea<-match(1,is.sea)
  if (is.na(nearest.sea)){
    seacell<-NA 
  } else {seacell<-x[nearest.sea]}  
  return(seacell)
}
#####################################################################
# UPWIND>SST BLOCK 
# Finds nearest upwind sea cell within 'buffer' (20km) for a block of cells for a SINGLE wind direction 
# If no sea cell found returns value of...0 sea cells returned = NA
# Otherwise records upwind sst held in sst.r
# Input:  sstbuffer.r - sea surface temperatures for buffer region
#         gridbuffer.r (including buffer) - could be dem
#         gridblock.r - defines area block of interest within buffer - could be dem
#         direction assumed to be constant across area
#         distance = max distance a search for nearest sea cell (= buffer of 10km)
# called by: upwind.sst.block
upwind.sst<-function(sst.buffer,dem.block,direction,distance=10000)
{
  #print(paste("Distance = ",distance))
  x<-dim(sst.buffer)[1]
  y<-dim(sst.buffer)[2]
  step<-distance/res(sst.buffer)[1] # max number of cells from focal cell to be searched
  
  # create matrix holding sea temperature for sea cells and -999 for landcells
  cells<-getValues(sst.buffer,format="matrix")
  seacells<-ifelse(is.na(cells),-999,cells)
  #plot(raster(seacells,template=sst.buffer))
  
  store<-array(0,dim=c(dim(seacells)[1]-(2*step),dim(seacells)[2]-(2*step),step+1)) # 3d array to hold values for non-buffered region for every 'step'
  store[,,1]<-seacells[(step+1):(dim(seacells)[1]-step),(step+1):(dim(seacells)[2]-step)]
  
  for (i in 1:step)
  {
    xshift<-round(i*sin(direction*(pi/180)),0) 
    yshift<-round(i*cos(direction*(pi/180)),0)
    #print(paste("i: ",i, " xshift: ",xshift," yshift: ",yshift,sep=""))
    yshift<-yshift*(-1)
    store[,,(i+1)]<-seacells[(step+1+yshift):(dim(seacells)[1]-step+yshift),(step+1+xshift):(dim(seacells)[2]-step+xshift)]
    
  } # end
  
  # use first/last to find nearest sea cell??
  storev<-array(store,dim=c((dim(store)[1]*dim(store)[2]),step+1))
  
  upwind.sst<-apply(storev,1,findsst) # if no sea within buffer then = NA?
  upwind.sst<-matrix(upwind.sst,nrow=nrow(store),ncol=ncol(store)) 
  
  e<-extent(xmin(dem.block)-distance,xmax(dem.block)+distance,ymin(dem.block)-distance, ymax(dem.block)+distance)
  upwind.sst.r<-crop(raster(upwind.sst,template=crop(dem.buffer,e)),dem.block)
  upwind.sst.r<-raster::mask(upwind.sst.r,dem.block)
  #plot(upwind.sst.r,main=paste("nearest seacell where direction= ",direction,sep=""))
  
  return(upwind.sst.r)
}# end function

#####

#####

# Function to return upwind.sst for block given w.dir for block
upwind.sst.block<-function(wdir.block,sst.buffer,dem.buffer,dem.block, plotresult=FALSE)
{
  # define output raster
  sst.vals<-array(NA,dim=c(nrow(wdir.block),ncol(wdir.block)))
  # Round wind direction to nearest degree
  wdir.block<-round(wdir.block,0)
  # Calc number of unique wind dir - create raster stack of upwind.sst for each unique w.direction
  wdir.vals<-unique(round(getValues(wdir.block)[which(!is.na(getValues(wdir.block)))])) # every unique wind.val excluding NA and rounding to nearest degree
  
  # Create stack of upwind sst for block for every unique wind direction during time t
  upwind.vals<-array(0,dim=c(nrow(wdir.block),ncol(wdir.block),length(wdir.vals)))  
  for (i in 1:length(wdir.vals)){
    upwind.r<-upwind.sst(sst.buffer,dem.block,wdir.vals[i],distance=10000)
    upwind.vals[,,i]<-getValues(upwind.r, format="matrix")
  }
  # assign sst according to w.dir of cell
  wdir.class<-cbind(wdir.vals,1:length(wdir.vals)) # create 2col vector for reclassifying
  wdir.class.r<-reclassify(wdir.block,wdir.class)  #reclass so that values = layer of stack holding sst values 
  class.vals<-getValues(wdir.class.r, format="matrix")
  
  for (i in 1:length(wdir.vals)){
    sel<-which(class.vals==i, arr.ind=TRUE)
    upwind.i<-upwind.vals[,,i]
    sst.vals[sel]<-upwind.i[sel]
    #plot(raster(sst.vals,template=wdir.block))
  }
  result<-raster(sst.vals,template=wdir.block)
  if (plotresult) plot(result)
  return(result)
} # end function

#######
# Function to calculate sst-tref 
sst.tref<-function(sst.r,tref.r)
{
  if (compareRaster(sst.r,tref.r)!=TRUE){warning("!!sst-tref rasters not comparable!!")}
  sst.tref.r<-overlay(sst.r,tref.r,fun=function(x,y){x-y})
  sst.tref.r<-calc(sst.tref.r,fun=function(x){ifelse(is.na(x),0,x)})
  return(sst.tref.r)
} # end function

#####################################################################
# WIND_DOWNSCALE FUNCTIONS
#####################################################################
# REQUIRES: Function JDdmy for computing Julian data to work out number of days after 1st Jan 1960 (1st wind data observation),
# so that correct element of array can be extracted

# Based on hour, day, month and year, extracts the required value from the array of values stored by wind_downscale1.R
# Inputs:
# hr: the hour (0-23) 0 = midnight
# day: the day of the month (0-31)
# month: the month of the year (numeric: 0-12)
# year: any year from 1960 to 2014
# Output: the element of the array stored by wind_downscale1.R that corresponds to either that hour, or the latest
# period immediatly before that hour (data only available 6-hourly)
array.val<-function(hr,day,month,yr)
{
  jd.base=JDdmy(1,1,1960)
  jd<-JDdmy(day,month,yr)
  dval<-(jd-jd.base)*4
  hval<-floor(hr/6)
  val<-dval+hval+1
  val
}

####################################################
# A. Prepare and load data for all time periods - Requires dem.block and e.block
####################################################
# 1. Load shelter maps for block using interval to which wind direction is rounded - eg 10 = every 10 degrees
# using interval, e.block, dem.block
block.sheltermap<-function(dem.block,dir_shelter,interval){
  dem.m<-getValues(dem.block,format="matrix")
  shelter<-array(NA, dim=c((360/interval),nrow(dem.m),ncol(dem.m)))
  for (i in 1:(360/interval)) {
    dir<-i*interval
    in.file<-paste(dir_shelter,"Shelter_",sprintf("%03d",dir,sep=""),"_deg.tif",sep="")
    print(in.file)
    wcoef.r<-raster(in.file) 
    projection(wcoef.r)<-CRS("+init=epsg:27700")
    wcoef.r<-crop(wcoef.r,dem.block)
    shelter[i,,]<-getValues(wcoef.r,format="matrix") # fills shelter[i,1:end,1] to shelter[i,1:10,end]
    print(paste("i= ",i," dir= ",dir))
  }   
  return(shelter)
} # end function


# 2. Load wind data (single file for whole time period - created by wind_downscale1)
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v
#load(file=paste(dir_wind,"wind_u.r",sep=""))
#load(file=paste(dir_wind,"wind_v.r",sep=""))


####################################################
# # # This Function downscales the wind to 100m cells for a specific hour
####################################################
# Specify hour,day, month and year for which data are required
#year=2010;month=6;day<-10;hr<-11
# Stages:
# (1) get wind values for a given hour, day, month and year
# (2) convert to 100m resolution raster OSGB grid reference using THIN SPLINET INTERPOLATION - see sepaate function
# (3) adjust based on altitude
# (4) adjust based on shelter coefficient
# OUTPUTS: wstr, wdir,invwstr and refwstr = wstr without elev/sheleter correction

wind.tpsdownscale<-function(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter,interval=10,print.results=FALSE,write.files=FALSE)
{
  # Can be quite slow. Allows you to keep tabs on progress by printing hour, day, month & year
  tp<-paste("year=",year," month=",month," day=",day," hour=",hr,sep="")
  print(tp)
  
  #############
  # Stage 1: get wind values for a given day month and year
  #############
  # As original data are 4x daily, but data are required for each hour,
  # this bit reads in the data for the periods immediatly before after for which there are data and calculates
  # weighted mean
  av1<-array.val(hr,day,month,year)
  av2<-av1+1
  rem<-hr/6-floor(hr/6)
  uwind1<-wind_u[,,av1]
  uwind2<-wind_u[,,av2]
  vwind1<-wind_v[,,av1]
  vwind2<-wind_v[,,av2]
  uwind<-(1-rem)*uwind1+rem*uwind2
  vwind<-(1-rem)*vwind1+rem*vwind2
  
  #############
  # Stage 2: convert to 100m resolution raster OSGB grid reference using tps resampling
  #############
  # Convert to raster (original lat long format and resolution - CHECK LAT/LON MAX/MIN
  uwind.r<-raster(uwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
  vwind.r<-raster(vwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
  # Reproject in OSGB projection
  crs(uwind.r)<-latlong
  crs(vwind.r)<-latlong
  u_osgb<-projectRaster(uwind.r,crs=ukgrid)
  v_osgb<-projectRaster(vwind.r,crs=ukgrid)
  # tps resampling to 100m using FUNCTION tps.resample
  u_100<-tps.resample(u_osgb,dem.block)
  v_100<-tps.resample(v_osgb,dem.block)
  
  #############
  # Stage 3: adjust based on altitude of terrain
  # NB Height adjustment based on wind spped values at different pressures downloaded  Earth System Research Lab
  # Typical heights at different pressures calculated from Allen et al 1998 http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric pressure (p)
  # Quadratic function fitted - NB this works well for heights up to ~1800m. IT won't work above ~2000m
  # Function was first derived by comparing values at different pressures (heights) over the course of a year (2014)
  #############
  # adjust wind speeds by height of dem
  # convert to matrices
  uwind.m<-getValues(u_100,format="matrix")
  vwind.m<-getValues(v_100,format="matrix")
  dem.m<-getValues(dem.block,format="matrix")
  # adjust wind by height
  ustr<-sqrt(uwind.m^2) # wind strength
  vstr<-sqrt(vwind.m^2) # wind strength
  udir<-ifelse(uwind.m>0,1,-1) # positive or negative
  vdir<-ifelse(vwind.m>0,1,-1) # positive or negative
  u.adj<-ustr*((-0.000000108025)*dem.m^2+0.000408692*dem.m+0.956139) # NB don't worry about warnings. Calculation assigns NAs to the sea
  v.adj<-vstr*((-0.000000108025)*dem.m^2+0.000408692*dem.m+0.956139)  # NB don't worry about warnings. Calculation assigns NAs to the sea
  # adjust values to correspond to wind speed 1m above theground
  # rescaling factor first derived by comparing values ot Culdrose wind data using wind_downscale3
  # however, in line wiht what you'd expect from: http://www.fao.org/docrep/x0490e/x0490e07.htm#wind profile relationship
  u.adj<-u.adj*0.373686439
  v.adj<-v.adj*0.373686439
  u.m<-u.adj*udir
  v.m<-v.adj*vdir
  u.r<-raster(u.m,template=dem.block) #only for checking
  v.r<-raster(v.m,template=dem.block)# only for checking
  # plot(u.r,main="u.r");plot(v.r,main="v.r")
  
  #############
  # Stage 4: height adjustments done using shelter coefficient maps based on topography and wind direction
  ############# 
  # Calculate Wind Direction
  dir.m <- (180/pi)*(atan2(u.m,v.m))  # NB this is direction in which wind blows to
  dir.m<-ifelse(dir.m<=180,dir.m+180,dir.m-180) # NB this direction from which wind originates (360 deg)
  dir.r<-raster(dir.m,template=dem.block)
  
  # Calculate Wind Strength from u.m and v.m components
  str.m<-matrix(NA,nrow=(nrow(dem.m)),ncol=(ncol(dem.m)) )# matrix for storing all values
  str.m<-sqrt(u.m^2+v.m^2)
  refwstr.r<-raster(str.m,template=dem.block)# USED FOR LATENT HEAT REF VALUES - WSTR WITHOUT SHELTER CORRECTION
  
  #  plot(wstr.r,main="wstr") ; plot(dir.r,main="dir")
  
  # Uses Rounded wind direction of each cell to select correct shelter coefficient
  dir.m<-round(dir.m/interval)*interval # round to interval used for shelter maps
  dir.m<-ifelse(dir.m==0,360,dir.m) # converts 0 to 360 degree direction 
  # Applies shelter coefficient to calculate wind strength if land cell (else NA)
  mxrws<-nrow(dem.m)
  mxcls<-ncol(dem.m)
  for (rws in 1:mxrws) {
    for (cls in 1:mxcls) {    
      str.m[rws,cls]<-str.m[rws,cls]*shelter[(dir.m[rws,cls]/interval),rws,cls] 
    }
  }  
  
  #############
  # Stage 5: Format outputs
  #############
  # Calculate inverse wind strength
  invstr.m<-1/(sqrt(str.m+1))
  
  # Convert to raster and save tif files 
  wstr.r<-raster(str.m,template=dem.block)
  wdir.r<-raster(dir.m,template=dem.block)
  invwstr.r<-raster(invstr.m,template=dem.block)
  if (write.files==TRUE){
    dir_windstrength<-paste(dir_wind,"strength/",sep="")
    dir_winddirection<- paste(dir_wind,"direction/",sep="")
    dir_windinvstr<- paste(dir_wind,"invstr/",sep="")
    fileout.1<-paste(dir_windstrength,"strength_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    fileout.2<-paste(dir_winddirection,"direction_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    fileout.3<-paste(dir_windinvstr,"invstr_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    #print(fileout.1); print(fileout.2); print(fileout.3)
    writeRaster(wstr.r,file=fileout.1,overwrite=TRUE)
    writeRaster(wdir.r,file=fileout.2,overwrite=TRUE)
    writeRaster(invwstr.r,file=fileout.3,overwrite=TRUE)
  }
  #Create raster stack and print if requested
  if (print.results==TRUE){
    shelterblock<-raster(shelter[(mean(dir.m,na.rm=TRUE)/interval),,],template=dem.block)
    result.stack<-stack(u.r,v.r,shelterblock,wstr.r,wdir.r,invwstr.r)
    names(result.stack)<-c("u.r","v.r","shelter block","wind strength","wind direction","inv wind str")
    par(mfrow=c(2,3))
    plot(result.stack)  
  } # end if
  
  wind.results<-c(wdir.r,wstr.r,invwstr.r,refwstr.r)
  return(wind.results)
  
} # end function wind.downscale



####################################################
# RADPROG FUNCTIONS
####################################################
# Converts a raster to a matrix for use with solar index functions NOT USED - USE getValues(r,format="matrix")
use.raster<-function(r)
{
  xr<-dim(r)[1]
  xc<-dim(r)[2]
  m<-array(getValues(r),dim=c(xc,xr))
  m<-t(m) # transpose
  m
}

# OS GB grid ref to Lat and Long
OSGBtolatlong<-function(x,y)
{
  pt = data.frame(x,y)
  coordinates(pt)=~x+y
  proj4string(pt)=CRS("+init=epsg:27700")
  latlong<-spTransform(pt,CRS("+init=epsg:4326"))
  ll<-as.data.frame(latlong)
  ll
}

# Needed for solar index function
solartime <- function(localtime,Long,Julian,merid=0,dst=0)
{
  Bn <- 2 * 3.141 * (Julian - 81) / 364
  eot <- 9.87 * sin(2 * Bn) - 7.53 * cos(Bn) - 1.5 * sin(Bn)
  solartime <- localtime + (4 / 60) * (2 * 3.141 * (merid - Long) / 360) + (eot / 60) - dst
  solartime
}

solalt <- function(localtime,Lat,Long,Julian,merid=0,dst=0)
{
  stime<-solartime(localtime,Long,Julian,merid,dst)
  tt <- 0.261799 * (stime - 12)
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((Julian - 171) / 365.25))
  Sinh = sin(declin) * sin(Lat * pi / 180) + cos(declin) * cos(Lat * 3.141 / 180) * cos(tt)
  solalt = (180 * atan(Sinh / sqrt(1 - Sinh * Sinh))) / pi
  solalt
}

solazi <- function(localtime,Lat,Long,Julian,merid=0,dst=0)
{
  stime<-solartime(localtime,Long,Julian,merid,dst)
  tt = 0.261799 * (stime - 12)
  declin = (pi * 23.5 / 180) * cos(2 * pi * ((Julian - 171) / 365.25))
  Sinh = sin(declin) * sin(Lat * pi / 180) + cos(declin) * cos(Lat * pi / 180) * cos(tt)
  hh = (atan(Sinh / sqrt(1 - Sinh * Sinh)))
  Sinazi = cos(declin) * sin(tt) / cos(hh)
  cosazi = (sin(Lat * pi / 180) * cos(declin) * cos(tt) - cos(pi * Lat / 180) * sin(declin)) / sqrt((cos(declin) *
                                                                                                       sin(tt)) ^ 2 + (sin(pi * Lat / 180) * cos(declin) * cos(tt) - cos(pi * Lat / 180) * sin(declin)) ^ 2)
  solazi = 180 + (180 * atan(Sinazi / sqrt(1 - Sinazi * Sinazi))) / pi
  if (cosazi < 0) {
    if (Sinazi < 0) {
      solazi = 180 - solazi
    } else {
      solazi = 540 - solazi
    }
  }
  solazi
}

horizonangle <- function(dtm,azimuth,res=100)
{
  dtm<-(dtm*5)/res
  azimuth<-azimuth-90
  azi <- azimuth * (pi/180)
  horizon <- array(0,dim(dtm))
  dtm3 <- array(0,dim(dtm)+200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)] <- dtm
  for (step in 1:10) {
    horizon[1:x,1:y] <- pmax(horizon[1:x,1:y], (dtm3[(101+sin(azi)*step^2):(x+100+sin(azi)*step^2),(101+cos(azi)*step^2):(y+100+cos(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(5*step^2))
  }
  horizon
}

solarindex <- function(slope,aspect,localtime,Lat,Long,Julian,dtm=array(0,dim=c(1,1)),res=100,merid=0,dst=0,shadow=TRUE)
{
  saltitude<-array(solalt(localtime,Lat,Long,Julian,merid,dst),dim(dtm))
  alt <- saltitude * (pi/180)
  zen <- pi/2 - alt
  sazimuth<-array(solazi(localtime,Lat,Long,Julian,merid,dst),dim(dtm))
  azi <- sazimuth * (pi/180)
  sl <- slope * (pi/180)
  asp <- aspect * (pi/180)
  shadowmask <- array(1,dim(dtm))
  horangle<-horizonangle(dtm,sazimuth)
  #plot(raster(horangle,template=dem.buffer))
  if(shadow) {
    shadowmask[horizonangle(dtm,sazimuth)>tan(alt)] <- 0
  }
  index <- array(0,dim(dtm))
  index <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - asp)
  index[index<0] <- 0
  index <- index * shadowmask
  #plot(raster(index,template=dem.buffer),main="Solar Index")
  index
}

skyview <-function(dtm,steps=36)
{
  sky <- array(1,dim(dtm))
  for (s in 1:steps) {
    sky <- sky-atan(horizonangle(dtm,s*360/steps))/((pi/2)*steps)
  }
  sky
}

####################################################################################
# Direct normal radiation: to downscaled direct radiation
####################################################################################

radiation_downscale_stack<-function(day,month,year,hr,sis.r,dnr.r, dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=FALSE)
{ 
  # Calculate lat and long of centre of grid
  ll<-OSGBtolatlong(xmin(dem.buffer)+(0.5*(xmax(dem.buffer)-xmin(dem.buffer))) , ymin(dem.buffer)+(0.5*(ymax(dem.buffer)-ymin(dem.buffer))) )
  lat<-as.numeric(ll[2])
  long<-as.numeric(ll[1])
  # Calculate if day or night from sunvector
  jd.h<-JDdmy(day,month,year)+((-12+hr)/24 )# Julian day with hour fraction
  if (sunvector(jd.h,lat,long,0)[3]>=0) # if DAYTIME 
  {
    # resample to 100m and extent of dem.buffer
    #dnr.buffer<-raster::resample(dnr.r,dem.buffer) # change this to different interpolation method?
    #sis.buffer<-raster::resample(sis.r,dem.buffer)
    # crop to 3* dem.buffer then tps model and interpolate
    #xdim<-(xmax(dem.buffer)-xmin(dem.buffer))
    #ydim<-(ymax(dem.buffer)-ymin(dem.buffer))
    #e.tps<-extent(xmin(dem.buffer)-xdim,xmax(dem.buffer)+xdim,ymin(dem.buffer)-ydim,ymax(dem.buffer)+ydim)# Run setup programs for creating constant raster maps etc
    e.tps<-extent(dem.buffer)
    sis.r<-crop(sis.r,e.tps); 
    dnr.r<-crop(dnr.r,e.tps); 
    sis.buffer<-tps.resample(sis.r,dem.buffer,FALSE)
    dnr.buffer<-tps.resample(dnr.r,dem.buffer,FALSE)
    print("Completed Rad resampling ") 

    
    #plot(sis.buffer); plot(dnr.buffer)
    
    # work out Julian day and time - ncdf time in hrs since 1/1/1983 0:00
    #jul.base<-JDdoy(1,1983)
    # hrs<-ncvar_get(ncdf_dnr,"time")
    #days<-floor(hrs/24)
    # jul.day<-as.numeric(jul.base+days) #JD at 12:00 on ncdf day
    # h<-as.numeric(hrs%%24)
    # if (h!=hr) {print("WARNING h^= hr!!!")}     
    
    # creates raster template for storing direct and diffuse radiation values with cell values of -9
    direct.r<-dem.buffer*0-9
    diffuse.r<-dem.buffer*0-9
    total.r<-dem.buffer*0-9
    reftotal.r<-dem.buffer*0-9
    
    # convert values to matrices for use with solar index function
    m.dem<-getValues(dem.buffer,format="matrix") # use.raster function above
    m.slope<-getValues(slope.buffer,format="matrix")
    m.aspect<-getValues(aspect.buffer,format="matrix")
    
    # converts NAs to zeros - NB: buffer of boundary blocks = NA->0
    sel<-which(is.na(m.dem)==T); m.dem[sel]<-0
    sel<-which(is.na(m.slope)==T); m.slope[sel]<-0
    sel<-which(is.na(m.aspect)==T); m.aspect[sel]<-0
    
    # Calculate solar index for dem 
    print("call solarindex")
    si<-solarindex(slope=m.slope,aspect=m.aspect,localtime=hr,
                   Lat=lat,Long=long,Julian=jd,dtm=m.dem)
    # Calculate single value solar index for flat slope
    slope.flat<-0;aspect.flat<-0
    si.flat<-solarindex(slope=slope.flat,aspect=aspect.flat,localtime=hr,
                        Lat=lat,Long=long,Julian=jd,shadow=F)[1,1]
    plot(raster(si,template=dem.buffer),main="Solar index")
    
    # Direct normal radiation:
    dnr.m<-getValues(dnr.buffer,format="matrix")
    # Direct radiation: flat
    dir.m<-dnr.m*si.flat # i.e. SID=Cos(sza)*DNR
    # Diffuse radiation flat - 
    sis.m<-getValues(sis.buffer,format="matrix")
    dif.m<-sis.m-dir.m
    #plot(raster(dir.m,template=dem.buffer),main="Direct - flat")
    #plot(raster(dif.m,template=dem.buffer),main="Diffuse - flat")
   
    
    # Calculate matrix coords that define central 'block'
    xmn<-1+(xmin(dem.block)-xmin(dem.buffer))/res(dem.block)[1]
    xmx<-dim(dem.buffer)[1]-(xmax(dem.buffer)-xmax(dem.block))/res(dem.block)[1]
    ymn<-1+(ymin(dem.block)-ymin(dem.buffer))/res(dem.block)[1]
    ymx<-dim(dem.buffer)[2]-(ymax(dem.buffer)-ymax(dem.block))/res(dem.block)[1]
    
    reftotal<-dir.m[xmn:xmx,ymn:ymx]+dif.m[xmn:xmx,ymn:ymx] # ref total rad without terrain effects
    
    # downscaled direct radiation
    direct<-dnr.m*si;   #  direct<-dir.m*si ?????
    direct<-direct[xmn:xmx,ymn:ymx] # Extract block
    #plot(raster(direct,template=dem.block),main="Old Direct")
    # downscaled diffuse radiation
    sv<-skyview(m.dem)
    diffuse<-dif.m*sv
    diffuse<-diffuse[xmn:xmx,ymn:ymx] # Extract block
    # add mask back in so that SIs are only produced for land
    mask<-getValues(dem.block,format="matrix")*0
    direct<-direct+mask
    diffuse<-diffuse+mask
    total<-direct+diffuse
    
    # Extract central Block and Convert to rasters
    direct.r<-raster(direct,template=dem.block)
    diffuse.r<-raster(diffuse,template=dem.block)
    total.r<-raster(total,template=dem.block) 
    reftotal.r<-raster(reftotal,template=dem.block)   #=total rad without terrain effects
      
  } else {         # IF NIGHTIME         
    direct.r<-dem.block*0
    diffuse.r<-dem.block*0
    total.r<-dem.block*0
    reftotal.r<-dem.block*0
  }
  
  #Create raster stack and print if requested
  if (print.results){
    result.stack<-stack(direct.r,diffuse.r,total.r,reftotal.r)
    names(result.stack)<-c("direct","diffuse","total","reftotal")
    par(mfrow=c(2,2))
    plot(result.stack)  
  } # end if
  
  results<-c(direct.r,diffuse.r,total.r,reftotal.r)
  return(results)
} # end function


###################################################################################



