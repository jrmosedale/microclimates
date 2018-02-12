## Downscale daily TEMPERATURE data to hourly  Prog: t5km_to_hrmatrix ###
#start.day<-1; start.month<-7; start.year<-1992
#end.day<-3; end.month<-7; end.year<-1992
start.jd<-JDdmy(start.day,start.month,start.year)
end.jd<-JDdmy(end.day,end.month,end.year)

##########################################################################################
# Main Function
# Writes daily 5km temperature files containing hourly temperature data for whole 5km area
# USES: fill.5km.map and relted FUNCTIONS from elevdif_map_functions - interpolates daily files then downscales
# i.e. records from nearest 5km historic data grid cell
#######################################################################################
hourly_temperatures<-function(start.jd,end.jd,dir_temp,dir_hrtemp,grid5km.r){ 
  plothrs=FALSE
  e.dem<-extent(grid5km.r)  
 
  for (jd in start.jd:end.jd) {  # Loop to calculate hourly files for each day
    ptm<-proc.time(); 
    
    # Define output matrix
    t100m.day<-array(0,dim=c(nrow(dem),ncol(dem),24))
    
    # Define files of current, previous and next day's temperature data
    tmaxfile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
    tminfile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_Actual.txt", sep="")
    prev.tmaxfile<-paste(dir_temp,"MaxTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_Actual.txt", sep="")
    #prev.tminfile<-paste(dir_temp,"MinTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_Actual.txt", sep="")
    #next.tmaxfile<-paste(dir_temp,"MaxTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_Actual.txt", sep="")
    next.tminfile<-paste(dir_temp,"MinTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_Actual.txt", sep="")
    
    # Check if files exist and assign current day file if no previous or next file
    if(file.exists(tmaxfile)==FALSE|file.exists(tminfile)==FALSE) stop(paste("No 5km data file for jd=",jd,sep=""))
    if(file.exists(prev.tmaxfile)==FALSE) {print("Warning: no prev data file"); prev.tmaxfile<-tmaxfile }
    if(file.exists(next.tminfile)==FALSE) {print("Warning: no prev data file");  next.tminfile<-tminfile}

    # Load data from files
    prev.tmax<-raster(prev.tmaxfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    #prev.tmin<-raster(prev.tminfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    day.tmax<-raster(tmaxfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    day.tmin<-raster(tminfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    #next.tmax<-raster(next.tmaxfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    next.tmin<-raster(next.tminfile,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Resample all daily rasters  using tps method and dem as model 
    day.tmax<-tps.resample(crop(day.tmax,dembuf),dem)
    day.tmin<-tps.resample(crop(day.tmin,dembuf),dem)
    #next.tmax<-tps.resample(crop(next.tmax,dembuf),dem)
    next.tmin<-tps.resample(crop(prev.tmax,dembuf),dem)
    prev.tmax<-tps.resample(crop(prev.tmax,dembuf),dem)
    #prev.tmin<-tps.resample(crop(prev.tmin,dembuf),dem)
    
    # 1. Calculate sunrise, sunset and day length for each grid cell
    # If first day then create lat/lon grid for daylength etc calculations from cropped temp data
    if (jd==start.jd){
      osgrid<-SpatialPoints(coordinates(day.tmax), proj4string=CRS("+init=epsg:27700"), bbox = NULL)
      os.m<-coordinates(osgrid)
      llgrid<-spTransform(osgrid,CRS("+init=epsg:4326"))
      ll.m<-coordinates(llgrid) # creates vector of coordinates from top left  grid cell to bottom right by rows
    }
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
    tmax<-getValues(day.tmax) # from top left to bottom right by row
    tmin<-getValues(day.tmin)
    p.tmax<-getValues(prev.tmax)
    n.tmin<-getValues(next.tmin)
    
    # Create hourly temperatures for each grid cell
    temp.hr<-generate.hrtemps.grid(tmin,tmax,n.tmin,p.tmax,sunup,sundown) 
    print(paste("Date: ", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),sep=""))
    
    # Plot hourly values
    for (hr in 0:23) {
      t100m.day[,,hr+1]<-matrix(temp.hr[,hr+1],nrow=NROW(t100m.day),ncol=NCOL(t100m.day),byrow=TRUE) # convert vector to matrix
      # plot raster option
      if (plothrs==TRUE){
        hr.r<-raster(t100m.day[,,hr+1],template=dem)
        temp.plot(hr.r,hr,DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year) }
    }
    
    # Write day of hrly data as 3D matrix
    file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),"_100m.r", sep="") # define file name from year,month,day,hr
    save(t100m.day, file=file.out)
    
    print(proc.time()-ptm) 
    
  } # end for day loop
  
}# end function

start.t<-proc.time()
hourly_temperatures(start.jd,end.jd,dir_temp,dir_hrtemp,grid5km.r) 
print(proc.time()-start.t)