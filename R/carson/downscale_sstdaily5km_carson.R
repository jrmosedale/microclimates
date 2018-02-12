##########################################################################################
# CALCULATE HOURLY 5KM Sea Surface Temperature 
##########################################################################################
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- as.integer(args[1] )
start.month<-as.integer(args[2] )
start.year<-as.integer(args[3]) 
end.day<-as.integer(args[4] )
end.month<-as.integer(args[5] )
end.year<-as.integer(args[6] )

## Downscale daily TEMPERATURE data to hourly  Prog: t5km_to_hrmatrix ###
# start.day<-1; start.month<-1; start.year<-1983
# end.day<-31; end.month<-12; end.year<-2013
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)
print(paste("Extracting data from ",start.day,"/",start.month,"/",start.year," to ",end.day,"/",end.month,"/",end.year,sep=""))

in.file<-paste(dir_sst,"HadISST_sst.nc",sep="") # in.file required for sst.spdownsc
ncfile<-nc_open(paste(dir_sst,"HadISST_sst.nc",sep="")) 

#######################################################################################
# FUNCTIONS TO WRITE DAILY FILES
#######################################################################################
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
sst.spdownsc<-function(start.year,start.month,end.year,end.month,dsgrid.r) {
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
    sst5km.r<-tps.resample(sst.os,dsgrid.r,maskoutput=FALSE)
    #sst5km.r<-raster::resample(sst.os,dsgrid.r,method="bilinear",na.rm=TRUE)
    #sst100m.r<-resample(sst.os,dem.buffer,method="bilinear",na.rm=TRUE) # alternative to 100m
    
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
sst.time.int<-function(jd,hr=12)
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


#######################################################################################
# CODE & FUNCTION CALL
#######################################################################################
### Downscale SST to daily 5km grid - Prog: sst_downscale_function ###
# IMPORTANT: requires 20KM BUFFER REGION to allow cal of upwind sst
# 1. WRITE monthly files for buffered uk region - normally 10km buffer - NEED DATA MONTH BEFORE START!!!
# Calculate jd for start of PREVIOUS month to start.jd *
# ?? Uses land mask of grid5km.r
if (start.month==1) { 
  year1<-start.year-1
  month1<-12 } else {
    year1<-start.year
    month1<-start.month-1 }
if (end.month==12) { 
  year2<-end.year+1
  month2<-1 } else {
    year2<-end.year
    month2<-end.month+1 }
print("Extracting monthly data")
sst.spdownsc(year1,month1,year2,month2,grid5kmbuf.r)  

# 2. WRITE daily files 
print("Writing daily data")
for (jd in start.jd:end.jd) {
  sst.time.int(jd)
} # end day


nc_close(in.file)