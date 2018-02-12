##########################################################################################
# DOWNSCALE INPUT DATA
# 
##########################################################################################

## Downscale daily TEMPERATURE data to hourly  Prog: t5km_to_hrmatrix ###
start.day<-1; start.month<-1; start.year<-1983
end.day<-1; end.month<-1; end.year<-1983

# includes data for previous day to start.day
hourly_temperatures(start.jd,end.jd,dir_temp,dir_hrtemp,grid5km.r) 

### Downscale direct and diffuse RADIATION - writes daily files
radiation_daystack(start.jd,end.jd,dir_sisday,dir_dniday)

#Downscale REL HUMIDITY to hourly 5km Prog: rel.hum.v2
# Writes rh.day as dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R"
rh.hourly(start.jd-1,end.jd,dir_rh,dir_rh5km,dir_hrtemp,grid5km.r)

# Downscale and impute CLOUD ALBEDO - WARNING - LONG TIME - RUN SEPERATELY??  Prog: calc.cal.hrly
# WRITES calimp.day as dir_calimp,"CALimp_5km_",year,"_",month,"_",day,".R
cal.daystack(start.jd,end.jd,dir_cal,dir_calday,plotcal=FALSE)
  
# Calculate LONG WAVE RADITAION at 5km hourly res ffrom CAL, RH T  Prog: longwav_grids
write_lwr_dayfiles(start.jd,end.jd,grid5km.r)

# Pressure - assign file - already unpacked
p.ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")


##########################################################################################

#####################################################################
# DATA PROCESSING FUNCTIONS
#####################################################################

# GENERATE HOURLY TEMPERATURE DATA AT %KM RESOLUTION

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



#######################################################################################
# Main Function
# Writes daily 5km temperature files containing hourly temperature data for whole 5km area
# USES: fill.5km.map and relted FUNCTIONS from elevdif_map_functions
# i.e. records from nearest 5km historic data grid cell
#######################################################################################
hourly_temperatures<-function(start.jd,end.jd,dir_temp,dir_hrtemp,grid5km.r){ 
  plothrs=FALSE
  e.dem<-extent(grid5km.r)  
  # Import temperature file available for day before start date (start.jd-1) 
  infile<-paste(dir_temp,"MaxTemp_", DMYjd(start.jd-1)$year, "-",sprintf("%02d",DMYjd(start.jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(start.jd-1)$day,sep=""),"_ACTUAL.txt", sep="")
  tmpdata.r<-raster(infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
  tmpdata.r<-crop(x=tmpdata.r,y=e.dem) # crop to geographical extent of DEM raster
  
  #Create lat/lon grid for daylength etc calculations from cropped temp data
  osgrid<-SpatialPoints(coordinates(tmpdata.r), proj4string=CRS("+init=epsg:27700"), bbox = NULL)
  os.m<-coordinates(osgrid)
  llgrid<-spTransform(osgrid,CRS("+init=epsg:4326"))
  ll.m<-coordinates(llgrid) # creates vector of coordinates from top left  grid cell to bottom right by rows
  
  for (jd in start.jd:end.jd) {  # Loop to calculate hourly files for each day
    # Define output matrix
    t5km.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
    # Read day temperature data
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
    day.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    day.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Read previous day temperature data
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_ACTUAL.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_ACTUAL.txt", sep="")
    prev.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    prev.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Read NEXT day file 
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="")
    next.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    next.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Crop all daily rasters  using dem as model - USE BRICK?? - COMBINE WITH BELOW?
    day.tmax<-crop(x=day.tmax,y=e.dem)
    day.tmin<-crop(x=day.tmin,y=e.dem)
    next.tmax<-crop(x=next.tmax,y=e.dem)
    next.tmin<-crop(x=next.tmin,y=e.dem)
    prev.tmax<-crop(x=prev.tmax,y=e.dem)
    prev.tmin<-crop(x=prev.tmin,y=e.dem)
    
    # OPTION - Fill missing 5km cells without temperature values but containing land cells from nearest neighbour
    day.tmax<-fill.5km.map(day.tmax,dem,grid5km.r) 
    day.tmin<-fill.5km.map(day.tmin,dem,grid5km.r) 
    next.tmax<-fill.5km.map(next.tmax,dem,grid5km.r) 
    next.tmin<-fill.5km.map(next.tmin,dem,grid5km.r) 
    prev.tmax<-fill.5km.map(prev.tmax,dem,grid5km.r) 
    prev.tmin<-fill.5km.map(prev.tmin,dem,grid5km.r) 
    
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
    tmax<-getValues(crop(x=day.tmax,y=dem)) # from top left to bottom right by row
    tmin<-getValues(crop(x=day.tmin,y=dem))
    p.tmax<-getValues(crop(x=prev.tmax,y=dem))
    n.tmin<-getValues(crop(x=next.tmin,y=dem))
    
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
        temp.plot(hr.r,hr,DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year) }
    }
    
    # Write day of hrly data as 3D matrix
    file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),".r", sep="") # define file name from year,month,day,hr
    save(t5km.day, file=file.out)
    
  } # end for day loop
  
}# end function


#######################################################################################
# FUNCTION GENREATES DAY STACK FILE OF HOURLY SIS AND DNI
# loads hourly files od dni and sis for one day
# Re-projection and crop but NO resampling
# Interpolates missing hours NA layers
# Writes day file of hourly layers
#######################################################################################
radiation_daystack<-function(start.jd,end.jd,dir_sisday,dir_dniday)
{  
  #dir_sisday<-"~/Documents/Exeter/Data2015/CMSAF-SIS/day/"
  #dir_dniday<-"~/Documents/Exeter/Data2015/CMSAF-DNI/day/"
  
  for (jd in start.jd:end.jd){
    day<-DMYjd(jd)$day
    month<-DMYjd(jd)$month
    year<-DMYjd(jd)$year
    
    dnr.stack<-stack()
    sis.stack<-stack()
    par(mfrow=c(3,4))
    
    # Read hourly data and create daily stack
    for (hr in 0:23){
      datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
      print(datetime)
      
      infile.dnr<-paste(dir_dni,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
      infile.sis<-paste(dir_sis,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
      print(paste(infile.dnr,"  ",infile.sis,sep=""))
      
      # Read in data from ncdf file and add to stack
      dnr.stack<-stack(dnr.stack,raster(infile.dnr))   
      sis.stack<-stack(sis.stack,raster(infile.sis))   
      #plot(raster(infile.dnr))
      
    }# end hr loop creating stack
    
    # Reproject and crop to UK area - keep original resolution
    projection(dnr.stack)<-"+init=epsg:4326"
    projection(sis.stack)<-"+init=epsg:4326"
    # reproject to OSGB and set extent to same as DEM
    dnr.stack<-projectRaster(dnr.stack,crs="+init=epsg:27700")
    sis.stack<-projectRaster(sis.stack,crs="+init=epsg:27700")
    dnr.stack<-crop(dnr.stack,demuk)
    sis.stack<-crop(sis.stack,demuk)
    
    # Interpolate NA layers
    dnr.int<-approxNA(dnr.stack,method="linear",rule=2)
    sis.int<-approxNA(sis.stack,method="linear",rule=2)   
    
    # To test loop
    for (lyr in 1:24){
      datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",(lyr-1),sep=""),":00",sep="")
      plot(raster(dnr.stack,lyr),main=datetime)
      if (compareRaster(raster(dnr.stack,lyr),raster(dnr.int,lyr), values=TRUE, stopiffalse=FALSE)==FALSE) print(paste("Interpolated values for: ",datetime,sep=""))
    }
    
    # Write files - raster stack by day  
    fileout1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
    fileout2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
    print(fileout1)
    print(fileout2)
    writeRaster(dnr.int,file=fileout1,overwrite=TRUE)
    writeRaster(sis.int,file=fileout2,overwrite=TRUE)
    
  } # end jd loop
  
} # end function

#####################################################################
# WRITE DAILY FILES OF RELATIVE HUMIDITY - rel.hum.v2 
#####################################################################
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
# Write daily files of hourly 5km Rel Humidity 
# RH data in = 4xdaily global data

rh.hourly<-function(start.jd,end.jd,dir_rh,dir_rh5km,dir_hrtemp, grid5km.r,hourplot=FALSE)
{
  e<-extent(-5.75,-0.75,48.75,53.75)  # lon/lat extent of cropped rasters
  #e<-extent(353.75,359.25,48.75,53.75)
  
  for (jd in start.jd:end.jd)
  {
    # Define output matrix to hold one day of hourly relative humidity values at 5km resolution for area covered by grid5km.r
    rh.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
    
    # Reads in Relative Humidity data from yearly file if new year
    if (jd==start.jd | (DMYjd(jd)$day==1 & DMYjd(jd)$month==1)){
      infile.rh<-paste(dir_rh,"rhum.sig995.",DMYjd(jd)$year,".nc",sep="") # reads yearly data file
      #print(in.file)
      rh.brick<-brick(infile.rh) # all times
      # get times
      netRH<-nc_open(infile.rh)
      tm<-ncvar_get(netRH,"time")
    }
    
    # Read in daily Temperature file and resample to low res lat/lon
    filein1<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),".r", sep="") # define file name from year,month,day,hr
    print(filein1)
    load(filein1)
    hr.temp1<-t5km.day
    
    filein2<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),".r", sep="") # define file name from year,month,day,hr
    print(filein2)
    load(filein2)
    hr.temp2<-t5km.day
    
    for (period in 0:3)
    {
      hr<-(period*6)
      # Identify rh layers and crop
      Jul.base<-JDdmy(1,1,1800)
      #Jul.actual<-JD(ISOdate(year,month,day))
      hr.val<-(jd-Jul.base)*24+hr
      sel<-which(tm==hr.val)
      rh.all0<-subset(rh.brick,sel)
      rh.all6<-subset(rh.brick,(sel+1))
      rh.all0<-reslice.raster(rh.all0)
      rh.all6<-reslice.raster(rh.all6)
      projection(rh.all0)<-"+init=epsg:4326"
      projection(rh.all6)<-"+init=epsg:4326"
      rh0<-crop(rh.all0,e)
      rh6<-crop(rh.all6,e)     
      
      # read in temperature rasters
      t0.5km<-raster(hr.temp1[,,hr+1],template=grid5km.r)
      if (period<3) t6.5km<-raster(hr.temp1[,,hr+7],template=grid5km.r)
      if (period==4) t6.5km<-raster(hr.temp2[,,1],template=grid5km.r) # 0hr00 for following day
      
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
    
    # WRITE DAILY files of 5km RELATIVE humidity
    fileout<-paste(dir_rh5km,"RH_5km_",DMYjd(jd)$year,"_",DMYjd(jd)$month,"_",DMYjd(jd)$day,".R",sep="")
    print(fileout)
    save(rh.day,file=fileout)
  } # end jd
  
} # end function

#####################################################################
# Imputes and write daily files of hourly data of CAL
#####################################################################
# Convert OS map extent to lat long extent
#### USES OSGBtolatlong FUNCTION
OStoLL.extent<-function(e){ 
  #e<-extent(dembuf)
  corner.xy<-c(OSGBtolatlong(xmin(e),ymin(e)) , OSGBtolatlong(xmax(e),ymin(e)) , OSGBtolatlong(xmin(e),ymax(e)) , OSGBtolatlong(xmax(e),ymax(e)) )
  xmn<-min(corner.xy[1]$x,corner.xy[3]$x,corner.xy[5]$x,corner.xy[7]$x)
  xmx<-max(corner.xy[1]$x,corner.xy[3]$x,corner.xy[5]$x,corner.xy[7]$x)
  ymn<-min(corner.xy[2]$y,corner.xy[4]$y,corner.xy[6]$y,corner.xy[8]$y)
  ymx<-max(corner.xy[2]$y,corner.xy[4]$y,corner.xy[6]$y,corner.xy[8]$y)
  e.ll<-extent(xmn,xmx,ymn,ymx)
  return(e.ll)
}

# FUNCTION WRITE DAILY STACK OF HOURLY IMPUTED CAL - NIGHT AND DAY
# Cropped to dembuf
# NO resampling or reprojection
# Output maps as 
cal.daystack<-function(start.jd,end.jd,dir_cal,dir_calday,plotcal=FALSE){
  dir_calday<-"~/Documents/Exeter/Data2015/CMSAF-CAL/day/"  
  par(mfrow=c(3,4))
  
  for (jd in start.jd:end.jd){
      cal.stack<-stack(); cal.int<-stack()
      # Start and end points at 12:00 (midday) of day before and after to permit interpolation of nighttime values
      # Read data for end of jd-1
      for (hr in 12:23){
        day<-DMYjd(jd-1)$day;  month<-DMYjd(jd-1)$month; year<-DMYjd(jd-1)$year
        datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
        #print(datetime)
        infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
        #print(infile.cal)   
        # Read in data from ncdf file and add to stack
        cal.stack<-stack(cal.stack,raster(infile.cal))   
        #plot(raster(infile.dnr))
      }# end hr loop creating stack
      
      # Read hourly data for jd
      for (hr in 0:23){
        day<-DMYjd(jd)$day;  month<-DMYjd(jd)$month; year<-DMYjd(jd)$year
        datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
        #print(datetime)
        infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
        #print(infile.cal)   
        # Read in data from ncdf file and add to stack
        cal.stack<-stack(cal.stack,raster(infile.cal))   
        #plot(raster(infile.dnr))
      }# end hr loop creating stack
      
      # Read data for end of jd+1
      for (hr in 0:11){
        day<-DMYjd(jd+1)$day;  month<-DMYjd(jd+1)$month; year<-DMYjd(jd+1)$year
        datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
        #print(datetime)
        infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
        #print(infile.cal)   
        # Read in data from ncdf file and add to stack
        cal.stack<-stack(cal.stack,raster(infile.cal))   
        #plot(raster(infile.dnr))
      }# end hr loop creating stack
      
      # Reproject and crop to UK area - keep original resolution
      projection(cal.stack)<-"+init=epsg:4326"
      # reproject to OSGB and set extent to same as DEM
      #cal.stack<-projectRaster(cal.stack,crs="+init=epsg:27700")
      e.ll<-OStoLL.extent(extent(dembuf))
      cal.stack<-crop(cal.stack,e.ll)
      # Interpolate NA layers
      cal.int<-approxNA(cal.stack,method="linear",rule=2)
      
      # To test loop
      if (plotcal){
        for (lyr in 1:nlayers(cal.int)){
          datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",(lyr-1),sep=""),":00",sep="")
          plot(raster(cal.int,lyr),main=lyr)
          if (compareRaster(raster(cal.stack,lyr),raster(cal.int,lyr), values=TRUE, stopiffalse=FALSE)==FALSE) print(paste("Interpolated values for: ",datetime,sep=""))
        } }
      
      cal.int<-stack(cal.int,layers=c(13:36)) # keep only one day of data
      #cal.int<-projectRaster(cal.int,crs="+init=epsg:27700") # reproject to OSGB and set extent to same as DEM
     
      #Write files - raster stack by day  
      day<-DMYjd(jd)$day;  month<-DMYjd(jd)$month; year<-DMYjd(jd)$year
      fileout1<-paste(dir_calday,"CALimp",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
      print(fileout1)
      writeRaster(cal.int,file=fileout1,overwrite=TRUE)
  }# end jd loop
} # end function

#####################################################################
# CALCULATE LW RADIATION from CAL , Temp, RH
# Units of output = MJ/m2
# Historic data: run rel.hum.v2, cal.cal.hrly & t5km_to_hrmatrix
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

write_lwr_dayfiles<-function(start.jd,end.jd,grid5km.r,plotlwr=FALSE)
{
  # Define output file - one day of hourly 5km data
  lwr.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
  for (t in start.jd:end.jd)
  {
    year<-DMYjd(t)$year
    month<-DMYjd(t)$month
    day<-DMYjd(t)$day
    
    # Reads in daily files - all assumed to be 5km OSGB gridcells with identical extents of grid5km.r
    rh.filein<-paste(dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R",sep="")
    print(rh.filein)
    load(rh.filein) # rh.day
    
    t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),".r", sep="")
    print(t.filein)
    load(t.filein) # t5km.day
    
    cal.filein<-paste(dir_calday,"CALimp",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
    print(cal.filein)
    cal.day<-stack(cal.filein)
    # Reproject to OSGB  and DOWNSCALE to 5km grid
    cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")
    cal.day<-resample(cal.day,grid5km.r)
    
    for (hr in 1:24)
    {
      # rel hum
      m.rh<-rh.day[,,hr]
      r.rh<-raster(m.rh,template=grid5km.r)
      # temperature
      m.temp<-t5km.day[,,hr]
      r.temp<-raster(m.temp,template=grid5km.r)
      # CAL 
      #m.cal<-CALimp.day[,,hr]
      r.cal<-raster(cal.day,layer=hr)
      
      # Calculate longwave radiation
      lwr.day[,,hr]<-lwr(getValues(r.temp,format="matrix"),
                         getValues(r.rh,format="matrix"),
                         getValues(r.cal,format="matrix"))
      
      if (plotlwr==TRUE) {
        par(mfrow=c(2,2))
        dayhr<-paste(day,"/",month,"/",year," ",hr-1,"h00",sep="")
        plot(r.cal,main=paste("Effective cloud albedo ",dayhr,sep=""))
        plot(r.rh,main=paste("Relative humdidity ",dayhr,sep=""))
        plot(r.temp,main=paste("Temperature ",dayhr,sep=""))
        plot(raster(lwr.day[,,hr],template=r.rh),main=paste("Long-wave radiation ",dayhr,sep=""))
      }
    } # for hr
    # WRITE DAILY lwr files
    file.out<-paste(dir_lwr,"lwr_",year,"_",month,"_",day,".R",sep="")
    print(file.out)
    save(lwr.day,file=file.out)
    
  } # for day
} # end function

#####################################################################
# DOWNSCALE  SST
#####################################################################
# Input: single ncdf file of Met Office Hadley Centre HadISST1 data 
# Output: Downscaled SST to hourly 100m data
# Issues/Problems: interpolation/resampling - in effect across N/S coasts
# set buffer eg to 10km
# Input data
#dir_sst<-"~/Documents/Exeter/Data2015/sst/"
#dir_sst<-"C:/Data2015/SST/"  
#gunzip(paste(dir_sst,"HadISST_sst.nc.gz",sep=""))
#dir_sstm<-paste(dir_sst,"monthly/",sep="")
#dir_ssth<-paste(dir_sst,"hourly/",sep="")
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
sst.spdownsc<-function(start.year,start.month,end.year,end.month,dsgrid.r,dir_sstm) {
  start.n<-(start.year-1960)*12+start.month
  end.n<-((end.year-1960)*12)+end.month
  for (n in start.n:end.n){
    # Read ncdf direct to raster
    sst.world.r<-raster(in.file,band=n) # read all of level t of ncdf file
    #plot(sst.world.r)
    # Extract relevant 1 deg cells 353-360E, 48:53N, 1960-2015 
    long.mn<--7;long.mx<--1 # -7 to -1
    lat.mn<-49; lat.mx<-52 # 49 to 52
    e.sst<-extent(c(long.mn,long.mx,lat.mn,lat.mx))
    sst.r<-crop(sst.world.r,e.sst)
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
    sst5km.r<-tps.resample(sst.os,dsgrid.r,msk=FALSE)
    #sst5km.r<-raster::resample(sst.os,dsgrid.r,method="bilinear",na.rm=TRUE)
    # sst100m.r<-resample(sst.os,dem.buffer,method="bilinear",na.rm=TRUE) # alternative to 100m
    
    #sst5km.r<-mask(sst5km.r,dsgrid.r,inverse=TRUE) # No masking to allow future downscaling to higher coastal resolutions
    # sst100m.r<-mask(sst5km.r,dsgrid.r,inverse=TRUE)  ; plot(sst100m.r)
    
    # write raster file at 5km resolution for every month
    tl<-paste("Year: ",ttoyr(n)," Month: ",ttomonth(n),sep="")
    plot(sst5km.r,main=tl)
    fileout<-paste(dir_sstm,"sst_",ttoyr(n),"_",ttomonth(n),".tif",sep="")
    print(fileout)
    writeRaster(sst5km.r,file=fileout,overwrite=T)
  }      
}# end of function

##########################################################################################
# 2. TEMPORAL DOWNSCALING FUNCTIONS to DAILY  time period (could also be used for hourly period if necessary)

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
sst.time.int<-function(jd,dir_sstm,dir_ssth,hr=12)
{
  # convert julian date to day/month/year
  day<-DMYjd(jd)$day ; month<-DMYjd(jd)$month ; year<-DMYjd(jd)$year
  
  # How many days in month
  dimth<-nday.in.month.jd(jd)
  
  # Read in SSTs for the two months that lie either side of the date in question
  # reads in month before if date is in first half of month - using jd allows for month being in different year
  mth1<-DMYjd(jd-(dimth/2))$month 
  mth2<-DMYjd(jd+(dimth/2))$month 
  
  filein1<-paste(dir_sstm,"sst_",year,"_",mth1,".tif",sep="")
  filein2<-paste(dir_sstm,"sst_",year,"_",mth2,".tif",sep="")
  r1<-raster(filein1)
  r2<-raster(filein2)
  v1<-getValues(r1,format="matrix")
  v2<-getValues(r2,format="matrix")
  
  # work out weighting to attach
  
  # How many days in both months
  dimth1<-nday.in.month.jd(jd-(dimth/2))
  dimth2<-nday.in.month.jd(jd+(dimth/2))
  #print(paste(day,"/",month, ".  dimth1= ",dimth1," dimth2=", dimth2,sep=""))
  
  # time after dimth1
  tt<-day+hr/24+dimth1/2
  if (day+hr/24>dimth/2) tt<-day+hr/24-dimth1/2
  wgt<-tt/((dimth1+dimth2)/2)
  v<-v2*wgt+v1*(1-wgt)
  ssthr.r<-raster(v,template=r1)
  
  # OPTION - downscale to 100m dembuf and mask
  #ssthr.r<-resample(ssthr.r,dembuf)
  #ssthr.r<-mask(sst.test,dembuf,inverse=TRUE)
  
  # write daily file
  tl<-paste("Year: ",year," Month: ",month," Day: ",day," Hour: ",hr,sep="")
  #plot(ssthr.r,main=tl)
  fileout<-paste(dir_ssth,"sst_",year,"_",month,"_",day,"_",hr,"h.tif",sep="")
  print(fileout)
  writeRaster(ssthr.r,file=fileout,overwrite=T)
  
} # end function

