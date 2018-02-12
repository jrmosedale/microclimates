##########################################################################################
# TEMPORAL DOWNSCALE AND INTERPOLATION OF RADIATION DATA
# Output: Daily files of hourly interpolated data at original resolution on OS projection
# Interpolation uses: Raster approx NA function
##########################################################################################
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- as.integer(args[1])
start.month<-as.integer(args[2])
start.year<-as.integer(args[3] )
end.day<-as.integer(args[4] )
end.month<-as.integer(args[5] )
end.year<-as.integer(args[6] )

#start.day<-1; start.month<-1; start.year<-1983
#end.day<-3; end.month<-1; end.year<-1983
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)

### Uses stored ecamples of empty and zero rasters for each data type
emptydni<-"/home/ISAD/jm622/Data2015/CMSAF-DNI/empty/DNIemptyUK.nc"
#emptydni<-"~/Documents/Exeter/Data2015/CMSAF-DNI/empty/DNIemptyUK.nc"
emptydni.r<-raster(emptydni)
emptysis<-"/home/ISAD/jm622/Data2015/CMSAF-SIS/empty/SISemptyUK.nc"
#emptysis<-"~/Documents/Exeter/Data2015/CMSAF-SIS/empty/SISemptyUK.nc"
emptysis.r<-raster(emptysis)
#Define zero stack
zero.dni<-"/home/ISAD/jm622/Data2015/CMSAF-DNI/empty/DNIzeroUK.nc"
#zero.dni<-"~/Documents/Exeter/Data2015/CMSAF-DNI/empty/DNIzeroUK.nc"
zero.sis<-"/home/ISAD/jm622/Data2015/CMSAF-SIS/empty/SISzeroUK.nc"
#zero.sis<-"~/Documents/Exeter/Data2015/CMSAF-SIS/empty/SISzeroUK.nc"
zerodni.r<-raster(zero.dni)
zerosis.r<-raster(zero.sis)
#par(mfrow=c(3,4))

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



##########################################################################################
# Create daily stacks
##########################################################################################

for (jd in start.jd:end.jd){
  day<-DMYjd(jd)$day; month<-DMYjd(jd)$month; year<-DMYjd(jd)$year
  
  ### 1. Load 24h stack of rasters and crop to lat lon of UK 
  dnr.24h<-stack()
  sis.24h<-stack()
  for (hr in 0:23){
    datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
    print(datetime)
    infile.dnr<-paste(dir_dni,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    infile.sis<-paste(dir_sis,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    #print(paste(infile.dnr,"  ",infile.sis,sep=""))
    dnr.24h<-stack(dnr.24h,raster(infile.dnr))   
    sis.24h<-stack(sis.24h,raster(infile.sis))   
  }# end hr loop creating stack  
  #plot(sis.24h[[1:12]],main="Before any replacement");  plot(sis.24h[[13:24]],main="Before any replacement")
  
  # Confirm original projection and crop to area of UK in lat lon
  projection(dnr.24h)<-"+init=epsg:4326"
  projection(sis.24h)<-"+init=epsg:4326"
  e.ukll<-extent(-7,2,49,61)
  dnr.24h<-crop(dnr.24h,e.ukll)
  sis.24h<-crop(sis.24h,e.ukll)
  
  ### 2. Interpolate or replace missing rasters
  # Interpolate missing rasters where possible
  dnr.24h<-replace.missing(dnr.24h,emptydni.r,zerodni.r)
  sis.24h<-replace.missing(sis.24h,emptysis.r,zerosis.r)
  #plot(sis.24h,main="After initial replacement")
  # If still missing layers then assign values of previous day
  for (lyr in 1:nlayers(dnr.24h)) {
    if(compareRaster(dnr.24h[[lyr]],emptydni.r,values=TRUE,stopiffalse=FALSE)==TRUE) {
      print(paste("Replacing DNI data at hr= ",lyr-1," with previous day data",sep=""))
      #p.day<-DMYjd(jd-1)$day; p.month<-DMYjd(jd-1)$month;p.year<-DMYjd(jd-1)$year
      #previous.dnr<-paste(dir_dniday,"DNIhm",p.year,sprintf("%02d",p.month,sep=""),sprintf("%02d.tif",p.day,sep=""),sep="")
      #dnr.24h[[lyr]]<-raster(previous.dnr,band=lyr)
      dnr.24h[[lyr]]<-subset(prevdnr.24h,lyr)
    }
    if(compareRaster(sis.24h[[lyr]],emptysis.r,values=TRUE,stopiffalse=FALSE)==TRUE) {
      print(paste("Replacing SIS data at hr= ",lyr-1," with previous day data",sep=""))
      #p.day<-DMYjd(jd-1)$day; p.month<-DMYjd(jd-1)$month;p.year<-DMYjd(jd-1)$year
      #previous.sis<-paste(dir_sisday,"SIShm",p.year,sprintf("%02d",p.month,sep=""),sprintf("%02d.tif",p.day,sep=""),sep="")
      #sis.24h[[lyr]]<-raster(previous.sis,band=lyr)
       sis.24h[[lyr]]<-subset(prevsis.24h,lyr)
    }
  } 
  # Final check and print warning if missing layers
  if (num.missing(dnr.24h,emptydni.r)>0) print("WARNING - missing layers remain for DNI data !!!")
  if (num.missing(sis.24h,emptysis.r)>0) print("WARNING - missing layers remain for SIS data !!!")  
  
  # Save stack if needed to substitue next day's raster
  prevdnr.24h<-dnr.24h
  prevsis.24h<-sis.24h
  
  ### 4. Write files - raster stack by day - FORMAT
  fileout1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
  fileout2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
  print(paste("Writing: ",fileout1,"   ", fileout2,sep=""))
  writeRaster(dnr.24h,file=fileout1,format="GTiff",overwrite=TRUE)
  writeRaster(sis.24h,file=fileout2,format="GTiff",overwrite=TRUE)
  
} # end for day loop



