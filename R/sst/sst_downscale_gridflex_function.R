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
  plot(sstds.r,main=tl)
  
  return (sstds.r)

} # end function

# Call example:
# sst.day<-sst.time.int(jd,dsgrid.r, dir_sstm,dir_ssth,hr=12)
