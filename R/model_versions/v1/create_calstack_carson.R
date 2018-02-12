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

max.jd<-JDdmy(31,12,2013) # = last date for which there is data
print(paste("Max JD= ",max.jd,sep=""))
#start.day<-4; start.month<-12; start.year<-2013
#end.day<-7; end.month<-12; end.year<-2013
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)
#dir_cal<-"~/Documents/Exeter/Data2015/CMSAF-CAL/extract/"
#dir_calday<-"~/Documents/Exeter/Data2015/CMSAF-CAL/day/"

### Uses stored ecamples of empty and zero rasters for each data type
emptycal<-"/home/ISAD/jm622/Data2015/CMSAF-CAL/empty/CALemptyUK"
#emptycal<-"~/Documents/Exeter/Data2015/CMSAF-CAL/empty/CALemptyUK"
emptycal.r<-raster(emptycal)

#Define zero stack
zerocal<-"/home/ISAD/jm622/Data2015/CMSAF-CAL/empty/CALzeroUK"
#zerocal<-"~/Documents/Exeter/Data2015/CMSAF-CAL/empty/CALzeroUK"
zerocal.r<-raster(zerocal)
#par(mfrow=c(3,4))

#################################
# FUNCTIONS
#################################
num.missing<-function(r.stack,empty.r){
  num.missing<-0
  for (lyr in 1:nlayers(r.stack)) {
    if(compareRaster(r.stack[[lyr]],empty.r,values=TRUE,stopiffalse=FALSE)==TRUE)  num.missing<-num.missing+1
  }
  return(num.missing)
}
##########################################################################################
# Create daily stacks
##########################################################################################

max.na<-0.35 # = maximum proportion of NA that is used to interpolate values - low because interested in SW corner of raster

for (jd in start.jd:end.jd){
  cal.s<-stack()
  
  ### 1. load last raster of previous (interpolated) day or set to NA
  if (jd>start.jd) {prev.dayfile<-paste(dir_calday,"CALhm",DMYjd(jd-1)$year,sprintf("%02d",DMYjd(jd-1)$month,sep=""),sprintf("%02d.tif",DMYjd(jd-1)$day,sep=""),sep="")
  prev.r<-raster(prev.dayfile,band=24) } else prev.r<-emptycal.r
  cal.s<-stack(cal.s,prev.r)
  
  ### 2. Load 24h stack of rasters and crop to lat lon of UK converting zero raster to NA
  for (hr in 0:23){
    datetime<-paste(DMYjd(jd)$year,"/",sprintf("%02d",DMYjd(jd)$month,sep=""),"/",sprintf("%02d",DMYjd(jd)$day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
    print(datetime)
    infile.cal<-paste(dir_cal,"CALhm",DMYjd(jd)$year,sprintf("%02d",DMYjd(jd)$month,sep=""),sprintf("%02d",DMYjd(jd)$day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    cal.r<-raster(infile.cal)
    # Confirm original projection and crop to area of UK in lat lon
    projection(cal.r)<-"+init=epsg:4326"
    e.ukll<-extent(-7,2,49,61)
    cal.r<-crop(cal.r,e.ukll)
    # Convert zero (night-time) rasters or those with high proportion of NA to missing NA rasters
   # if(compareRaster(cal.r,zerocal.r,values=TRUE,stopiffalse=FALSE)==TRUE) cal.r<-emptycal.r 
    #if( length(which(is.na(getValues(cal.r))) ) / ncell(cal.r)>max.na) cal.r<-emptycal.r 
    cal.s<-stack(cal.s,cal.r)   
  }# end hr loop creating stack  
  #plot(cal.s[[1:12]],main="Before any replacement");  plot(cal.s[[13:24]],main="Before any replacement")

  ### 3. Load additional rasters till next non NA and non-zero raster
  next.day<-jd+1
  hr<-0
  while ( length(which(is.na(getValues(cal.r))) ) / ncell(cal.r)>max.na & next.day<=end.jd){
    if (next.day-jd>31 | next.day>max.jd){print("WARNING Over 31 days without valid raster - QUITTING LOOP"); break}
    if (hr>23) { next.day<-next.day+1 ; hr<-0} # go to next day if still no valid raster
    infile.cal<-paste(dir_cal,"CALhm",DMYjd(next.day)$year,sprintf("%02d",DMYjd(next.day)$month,sep=""),sprintf("%02d",DMYjd(next.day)$day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    cal.r<-raster(infile.cal)
    # Confirm original projection and crop to area of UK in lat lon
    projection(cal.r)<-"+init=epsg:4326"
    e.ukll<-extent(-7,2,49,61)
    cal.r<-crop(cal.r,e.ukll)
    # Convert zero (night-time) rasters to missing NA rasters
    if(compareRaster(cal.r,zerocal.r,values=TRUE,stopiffalse=FALSE)==TRUE) cal.r<-emptycal.r   
    cal.s<-stack(cal.s,cal.r) 
    hr<-hr+1
  } # while loop
  print(paste("Loaded additional ",next.day-(jd+1)," days and ", hr,"hrs for interpolation"))
  
  ### 4. Interpolate missing rasters within stack - Rule allows for NA at either end to account for 1st/last days
  cal.s<-approxNA(cal.s,method='linear',rule=2:2) 
  cal.24h<-subset(cal.s,2:25) # extract core day of data  
  #plot(cal.24h[[1:12]],main=datetime); plot(cal.24h[[13:24]],main=datetime)
  # Final check and print warning if missing layers
  if (num.missing(cal.24h,emptycal.r)>0) print("WARNING - missing layers remain for DNI data !!!")
  
  ### 4. Write file for original day=jd
  fileout<-paste(dir_calday,"CALhm",DMYjd(jd)$year,sprintf("%02d",DMYjd(jd)$month,sep=""),sprintf("%02d.tif",DMYjd(jd)$day,sep=""),sep="")
  print(paste("Writing: ",fileout,sep=""))
  writeRaster(cal.24h,file=fileout,format="GTiff",overwrite=TRUE)
  
} # end for jd loop




