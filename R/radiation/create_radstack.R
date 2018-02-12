

#################################
# FUNCTIONS
#################################
replace.missing<-function(r.stack,empty.r,zero.r)
{
  # Calculate earliest sunrise and latest sunset in raster
  sunup<-min(sunrise.grid(jd,ymax(empty.r),xmax(empty.r),0,0), sunrise.grid(jd,ymin(empty.r),xmax(empty.r),0,0))
  sunset<-max(sunset.grid(jd,ymax(empty.r),xmin(empty.r),0,0), sunset.grid(jd,ymin(empty.r),xmin(empty.r),0,0))
  print(paste("Sunrise= ",sunup," Sunset= ",sunset,sep=""))
  
  print(paste("Missing layers before interpolation: ",num.missing(r.stack,empty.r),sep=""))
  #if (num.missing(r.stack,empty.r)>0 ) print(missing.layers(r.stack,empty.r) )
  # Night-time defined for most easterly then westerly location in raster - replace any NA rasters with zero rasters
  for (lyr in 1:(floor(sunup)+1)) 
    if(compareRaster(r.stack[[lyr]],empty.r,values=TRUE,stopiffalse=FALSE)==TRUE) r.stack[[lyr]]<-zero.r
  for (lyr in (1+ceiling(sunset)):24) 
    if(compareRaster(r.stack[[lyr]],empty.r,values=TRUE,stopiffalse=FALSE)==TRUE) r.stack[[lyr]]<-zero.r
  
  # Apply NA interpolation if not too many missing layers
  am<-subset(r.stack,1:15)
  am<-approxNA(am,method='linear',rule=1)
  r.stack<-stack(am,r.stack[[16:24]])
  pm<-subset(r.stack,12:24)
  pm<-approxNA(pm,method='linear',rule=1)
  r.stack<-stack(r.stack[[1:11]],pm)
  return(r.stack)
}# end function replace.missing

# Create T/F vector of whether different layers = NA
missing.layers<-function(r.stack,empty.r){
  missing<-rep(FALSE,nlayers(r.stack))
  for (lyr in 1:nlayers(r.stack)) {
    if(compareRaster(r.stack[[lyr]],empty.r,values=TRUE,stopiffalse=FALSE)==TRUE)  missing[lyr]<-TRUE
  }
  return(missing)
}# end function

num.missing<-function(r.stack,empty.r){
  num.missing<-0
  for (lyr in 1:nlayers(r.stack)) {
    if(compareRaster(r.stack[[lyr]],empty.r,values=TRUE,stopiffalse=FALSE)==TRUE)  num.missing<-num.missing+1
  }
  return(num.missing)
}

# Calculate sunrise and sunset raster (azimuth for each hour?)
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



#################################
# CODE to interpolate/replace missing layers of daily raster stacks
#################################
# Interpolate NA layers
#start.day<-1; start.month<-1; start.year<-1983
#end.day<-18 ; end.month<-3; end.year<-1985
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)
#jd<-start.jd

# Uses stored ecamples of empty and zero rasters for each data type
# Define empty layer for testing from saved files as models
empty<-"~/Documents/Exeter/Data2015/CMSAF-DNI/empty/DNIempty.nc"
emptydni.r<-raster(empty)
empty<-"~/Documents/Exeter/Data2015/CMSAF-SIS/empty/SISempty.nc"
emptysis.r<-raster(empty)
#Define zero stack
zero.dni<-"~/Documents/Exeter/Data2015/CMSAF-DNI/empty/DNIzero.nc"
zero.sis<-"~/Documents/Exeter/Data2015/CMSAF-SIS/empty/SISzero.nc"
zerodni.r<-raster(zero.dni)
zerosis.r<-raster(zero.sis)
#par(mfrow=c(3,4))

for (jd in start.jd:end.jd){
  # Read hourly data and create daily stack
  day<-DMYjd(jd)$day
  month<-DMYjd(jd)$month
  year<-DMYjd(jd)$year
  
  ### 1. Load 24h stack of rasters
  dnr.24h<-stack()
  sis.24h<-stack()
  for (hr in 0:23){
    datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
    print(datetime)
    infile.dnr<-paste(dir_dni,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    infile.sis<-paste(dir_sis,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    #print(paste(infile.dnr,"  ",infile.sis,sep=""))
    # Read in data from ncdf file and add to stack
    dnr.24h<-stack(dnr.24h,raster(infile.dnr))   
    sis.24h<-stack(sis.24h,raster(infile.sis))   
  }# end hr loop creating stack  
  #plot(sis.24h[[1:12]],main="Before any replacement");  plot(sis.24h[[13:24]],main="Before any replacement")

  ### 2. Interpolate or replace missing rasters
  # Replace / interpolate missing rasters where possible
  dnr.24h<-replace.missing(dnr.24h,emptydni.r,zerodni.r)
  sis.24h<-replace.missing(sis.24h,emptysis.r,zerosis.r)
  #plot(sis.24h,main="After initial replacement")
  
  # If still missing layers then assign values of previous day
  for (lyr in 1:nlayers(dnr.24h)) {
    if(compareRaster(dnr.24h[[lyr]],empty.r,values=TRUE,stopiffalse=FALSE)==TRUE) {
      print(paste("Replacing DNI data at hr= ",lyr+1," with previous day data",sep=""))
      prev.jd<-jd-1
      p.day<-DMYjd(prev.jd)$day; p.month<-DMYjd(prev.jd)$month;p.year<-DMYjd(prev.jd)$year
      previous.dnr<-paste(dir_dni,"DNIhm",p.year,sprintf("%02d",p.month,sep=""),sprintf("%02d",p.day,sep=""),sprintf("%02d",lyr+1,sep=""),"00002UD1000101UD.nc",sep="")
      dnr.24h[[lyr]]<-raster("previous.dnr")
    }
    if(compareRaster(sis.24h[[lyr]],empty.r,values=TRUE,stopiffalse=FALSE)==TRUE) {
      print(paste("Replacing SIS data at hr= ",lyr+1," with previous day data",sep=""))
      prev.jd<-jd-1
      p.day<-DMYjd(prev.jd)$day; p.month<-DMYjd(prev.jd)$month;p.year<-DMYjd(prev.jd)$year
      previous.sis<-paste(dir_sis,"SIShm",p.year,sprintf("%02d",p.month,sep=""),sprintf("%02d",p.day,sep=""),sprintf("%02d",lyr+1,sep=""),"00002UD1000101UD.nc",sep="")
      sis.24h[[lyr]]<-raster("previous.sis")
    }
  } 
  # Final check and print warning if missing layers
  if (num.missing(dnr.24h,emptydni.r)>0) print("WARNING - missing layers remain for DNI data !!!")
  if (num.missing(sis.24h,emptysis.r)>0) print("WARNING - missing layers remain for SIS data !!!")
  
  
  ### 3. Reproject and crop to UK area - keep original resolution
  projection(dnr.stack)<-"+init=epsg:4326"
  projection(sis.stack)<-"+init=epsg:4326"
  # reproject to OSGB and set extent to same as DEM
  dnr.stack<-projectRaster(dnr.stack,crs="+init=epsg:27700")
  sis.stack<-projectRaster(sis.stack,crs="+init=epsg:27700")
  dnr.stack<-crop(dnr.stack,demuk)
  sis.stack<-crop(sis.stack,demuk)
  
  ### 4. Write files - raster stack by day - FORMAT
  fileout1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
  fileout2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
  print(fileout1)
  print(fileout2)
  writeRaster(dnr.stack,file=fileout1,format="GTiff",overwrite=TRUE)
  writeRaster(sis.stack,file=fileout2,format="GTiff",overwrite=TRUE)
  
} # end for day loop








