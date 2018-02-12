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
cal.daystack<-function(dir_cal,dir_calday,jd,plotcal=FALSE){
  dir_calday<-"~/Documents/Exeter/Data2015/CMSAF-CAL/day/"  
  cal.stack<-stack(); cal.int<-stack()
  par(mfrow=c(3,4))
  
  # Start and end points at 12:00 (midday) of day before and after to permit interpolation of nighttime values
  # Read data for end of jd-1
  for (hr in 12:23){
    day<-DMYjd(jd-1)$day;  month<-DMYjd(jd-1)$month; year<-DMYjd(jd-1)$year
    datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
    print(datetime)
    infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    print(infile.cal)   
    # Read in data from ncdf file and add to stack
    cal.stack<-stack(cal.stack,raster(infile.cal))   
    #plot(raster(infile.dnr))
  }# end hr loop creating stack
  
  # Read hourly data for jd
  for (hr in 0:23){
    day<-DMYjd(jd)$day;  month<-DMYjd(jd)$month; year<-DMYjd(jd)$year
    datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
    print(datetime)
    infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    print(infile.cal)   
    # Read in data from ncdf file and add to stack
    cal.stack<-stack(cal.stack,raster(infile.cal))   
    #plot(raster(infile.dnr))
  }# end hr loop creating stack
  
  # Read data for end of jd+1
  for (hr in 0:11){
    day<-DMYjd(jd+1)$day;  month<-DMYjd(jd+1)$month; year<-DMYjd(jd+1)$year
    datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
    print(datetime)
    infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    print(infile.cal)   
    # Read in data from ncdf file and add to stack
    cal.stack<-stack(cal.stack,raster(infile.cal))   
    #plot(raster(infile.dnr))
  }# end hr loop creating stack
  
  # Reproject and crop to UK area - keep original resolution
  projection(cal.stack)<-"+init=epsg:4326"
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
 
  # Keep 24hrs of data for jd
  cal.int<-stack(cal.int,layers=c(13:36)) # keep only one day of data
  # Reproject to OSGB  and DOWNSCALE to 5km grid
  cal.int<-projectRaster(cal.int,crs="+init=epsg:27700")
  cal.5km<-resample(cal.int,grid5km.r)
  
  # Write files - raster stack by day 
  fileout1<-paste(dir_calday,"CALimp",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
  print(fileout1)
  writeRaster(cal.5km,file=fileout1,overwrite=TRUE)
  
} # end function




cal_downscale<-function(cal.r,downscale.r,print.results=TRUE,write.files=FALSE)
{
    #cal.r<-raster(cal.int,layer=13)
    #downscale.r<-grid5km.r
    cal.r<-projectRaster(cal.r,crs="+init=epsg:27700")
    #cal.r<-crop(cal.r,dembuf)
    xdim<-20000
    ydim<-20000
    e.tps<-extent(xmin(downscale.r)-xdim,xmax(downscale.r)+xdim,ymin(downscale.r)-ydim,ymax(downscale.r)+ydim)# Run setup programs for creating constant raster maps etc   
    cal.r<-crop(cal.r,e.tps)
    #cal.dsc<-tps.resample(cal.r,downscale.r,FALSE)
    cal.dsc<-resample(cal.r,downscale.r) # as effective as tps for downscale to 5km
    #plot(cal.dsc)
    return(cal.dsc)
} # end function


