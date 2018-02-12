##########################################################################################
# TEMPORAL DOWNSCALE AND INTERPOLATION OF CAL DATA
# Output: Daily files of hourly interpolated data
# INterpolation uses: Raster approx NA function
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


#######################################################################################
# FUNCTION GENREATES DAY STACK FILE OF HOURLY CAL DATA
# loads hourly files od dni and sis for one day
# Re-projection and crop but NO resampling
# Interpolates missing hours NA layers
# Writes day file of hourly layers
#######################################################################################
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
  #dir_calday<-"~/Documents/Exeter/Data2015/CMSAF-CAL/day/"  
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
    
    # Read data for start of jd+1
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


##########################################################################################
# FUNCTION CALL TO GENERATE  DATA
#######################################################################################
# Downscale and impute CLOUD ALBEDO - WARNING - LONG TIME - RUN SEPERATELY??  Prog: calc.cal.hrly
# WRITES calimp.day as dir_calimp,"CALimp_5km_",year,"_",month,"_",day,".R
cal.daystack(start.jd,end.jd,dir_cal,dir_calday,plotcal=FALSE)
