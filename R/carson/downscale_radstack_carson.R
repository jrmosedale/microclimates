##########################################################################################
# TEMPORAL DOWNSCALE AND INTERPOLATION OF RADIATION DATA
# Output: Daily files of hourly interpolated data at original resolution on OS projection
# INterpolation uses: Raster approx NA function
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

#start.day<-1; start.month<-1; start.year<-1983
#end.day<-31; end.month<-12; end.year<-2013
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)

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
      
    # Write files - raster stack by day - FORMAT
    fileout1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
    fileout2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
    print(fileout1)
    print(fileout2)
    writeRaster(dnr.stack,file=fileout1,overwrite=TRUE)
    writeRaster(sis.stack,file=fileout2,overwrite=TRUE)
    
  } # end jd loop
  
} # end function

##########################################################################################
# FUNCTION CALL TO GENERATE  DATA
#######################################################################################
### Downscale direct and diffuse RADIATION - writes daily files
radiation_daystack(start.jd,end.jd,dir_sisday,dir_dniday)
