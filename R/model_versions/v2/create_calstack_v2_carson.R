##########################################################################################
# TEMPORAL DOWNSCALE AND INTERPOLATION OF CAL DATA
# Output: Daily files of hourly interpolated data at original resolution on OS projection
# Interpolation uses: Raster approx NA function
# Assumes ONE YEAR JOBS - calculates and interpolates CAL across whole year
# ASSUMES begins at 0h00 and ends 23h00 ie no partial days
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

#start.day<-20; start.month<-6; start.year<-1992
#end.day<-31; end.month<-12; end.year<-1992
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)
#dir_cal<-"~/Documents/Exeter/Data2015/CMSAF-CAL/extract/"
#dir_calday<-"~/Documents/Exeter/Data2015/CMSAF-CAL/day/"

# Define crop area 
e.ukll<-extent(-7,2,49,53)

### Uses stored ecamples of empty and zero rasters for each data type
emptycal<-"/home/ISAD/jm622/Data2015/CMSAF-CAL/empty/CALempty"
#emptycal<-"~/Documents/Exeter/Data2015/CMSAF-CAL/empty/CALempty"
emptycal.r<-raster(emptycal)
#################################
# FUNCTIONS
num.missing<-function(r.stack,empty.r){
  num.missing<-0
  for (lyr in 1:nlayers(r.stack)) {
    if(compareRaster(r.stack[[lyr]],empty.r,values=TRUE,stopiffalse=FALSE)==TRUE)  num.missing<-num.missing+1
  }
  return(num.missing)
}
#################################

# 1. Create daily stacks
par(mfrow=c(3,4))
cal.s<-stack()

for (jd in start.jd:end.jd){ # load all files for a year
  # Calculate last sunrise and first sunset for map area
  sunrise<-ceiling(daylength(52,-7,jd,0)[1])+1 
  sunset<-floor(daylength(52,2,jd,0)[2])-1
  
  ### Load 24h stack of rasters and crop to lat lon of UK converting zero raster to NA
  for (hr in 0:23){
    datetime<-paste(DMYjd(jd)$year,"/",sprintf("%02d",DMYjd(jd)$month,sep=""),"/",sprintf("%02d",DMYjd(jd)$day,sep="")," ",sprintf("%02d",hr,sep=""),":00",sep="")
    print(datetime)
    infile.cal<-paste(dir_cal,"CALhm",DMYjd(jd)$year,sprintf("%02d",DMYjd(jd)$month,sep=""),sprintf("%02d",DMYjd(jd)$day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
    # Load and crop daytime CAL files or set to empty (using sunrise/set values)
    if (hr<sunrise|hr>sunset) {
      cal.r<-emptycal.r 
      } else {
        print(infile.cal)
        cal.r<-raster(infile.cal)
        projection(cal.r)<-"+init=epsg:4326"
        cal.r<-crop(cal.r,e.ukll)
       # plot(cal.r,main=datetime)
      }
    cal.s<-stack(cal.s,cal.r)   
  }# end hr loop   
} # end jd
  
# 2. Set 1st and last rasters to nearest valid daytime layer
# find first valid CAL data layer and set 1st layer to equal this
n<-1
while (compareRaster(raster(cal.s,layer=n),emptycal.r,values=TRUE,stopiffalse=FALSE)==TRUE & n<nlayers(cal.s)){
  n<-n+1
  if (n>48) print("More than 48hrs without data!!!")
 if (n>240) print("More than 10 DAYS without data!!!")
}
if (n==nlayers(cal.s)) print ("WARNING No valid data") else {
  new.r<-raster(cal.s,layer=n)
  cal.s<-dropLayer(cal.s,1)
  cal.s<-stack(new.r,cal.s)
  }
  
# find last valid CAL data layer and set last layer to equal this
n<-nlayers(cal.s)
while (compareRaster(raster(cal.s,layer=n),emptycal.r,values=TRUE,stopiffalse=FALSE)==TRUE & n>1){
  n<-n-1
} 
if (n==2) print ("WARNING Lack of valid data") else {
  new.r<-raster(cal.s,layer=n)
  cal.s<-dropLayer(cal.s,nlayers(cal.s))
  cal.s<-stack(cal.s,new.r)
}

### 3. Interpolate missing rasters within stack - Rule allows for NA at either end to account for 1st/last days
cal.s<-approxNA(cal.s,method='linear',rule=2:2) 

# Final check and print warning if missing layers
if (num.missing(cal.s,emptycal.r)>0) print("WARNING - missing layers remain for DNI data !!!")
if (nlayers(cal.s)%%24!=0) print("WARNING - Number of Layers NOt divisible by 24!!!")

### 4. Write DAILY files for original day=jd
n<-0
for (jd in start.jd:end.jd){
  cal.24h<-subset(cal.s,((n*24)+1):((n*24)+24) )
  fileout<-paste(dir_calday,"CALhm",DMYjd(jd)$year,sprintf("%02d",DMYjd(jd)$month,sep=""),sprintf("%02d.tif",DMYjd(jd)$day,sep=""),sep="")
  print(paste("Writing: ",fileout," using stack layers ",(n*24)+1," to ",(n*24)+24, sep=""))
  writeRaster(cal.24h,file=fileout,format="GTiff",overwrite=TRUE)
  n<-n+1
}

# 
# for (n in 4632:nlayers(cal.s)){
# plot(raster(cal.s,layer=n),main=n) }



