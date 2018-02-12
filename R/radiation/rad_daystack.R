# FUNCTION loads hourly files od dni and sis for one day
# Re-projection and crop but NO resampling
# Interpolates missing hours NA layers
# Writes day file of hourly layers

radiation_daystack<-function(day,month,year,dem.buffer)
{  
  dir_sisday<-"~/Documents/Exeter/Data2015/CMSAF-SIS/day/"
  dir_dniday<-"~/Documents/Exeter/Data2015/CMSAF-DNI/day/"
  
  dnr.stack<-stack()
  sis.stack<-stack()
  par(mfrow=c(3,4))
  datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",(lyr-1),sep=""),":00",sep="")
  
  # Read hourly data and create stack
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
  
  # Interpolate NA layers
  dnr.int<-approxNA(dnr.stack,method="linear",rule=2)
  sis.int<-approxNA(sis.stack,method="linear",rule=2)   

  # To test loop
  for (lyr in 1:24){
    plot(raster(dnr.stack,lyr),main=datetime)
    if (compareRaster(raster(dnr.stack,lyr),raster(dnr.int,lyr), values=TRUE, stopiffalse=FALSE)==FALSE) print(paste("Interpolated values for: ",datetime,sep=""))
  }
  
  # Write files - raster stack by day  
  fileout1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
  fileout2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
  print(fileout1)
  print(fileout2)
  writeRaster(dnr.int,file=fileout1,overwrite=TRUE)
  writeRaster(sis.int,file=fileout2,overwrite=TRUE)
  
} # end function


radiation_downscale_stack<-function(day,month,year,hr,demuk,dem.buffer,dem.block,slope.buffer,aspect.buffer,dir_rad,print.results=TRUE,write.files=FALSE)
{
  filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
  filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
  dnr.24h<-stack(filein1)
  sis.24h<-stack(filein2)
  dnr.r<-raster(dnr.24h,layer=hr+1);plot(dnr.r)
  sis.r<-raster(sis.24h,layer=hr+1); plot(sis.r)

  if (cellStats(sis.r,max)>0) # if DAYTIME
  {
    # resample to 100m and extent of dem.buffer
    #dnr.buffer<-raster::resample(dnr.r,dem.buffer) # change this to different interpolation method?
    #sis.buffer<-raster::resample(sis.r,dem.buffer)
    # crop to 3* dem.buffer then tps model and interpolate
    xdim<-(xmax(dem.buffer)-xmin(dem.buffer))
    ydim<-(ymax(dem.buffer)-ymin(dem.buffer))
    e.tps<-extent(xmin(dem.buffer)-xdim,xmax(dem.buffer)+xdim,ymin(dem.buffer)-ydim,ymax(dem.buffer)+ydim)# Run setup programs for creating constant raster maps etc

    sis.r<-crop(sis.r,e.tps)
    sis.buffer<-tps.resample(sis.r,dem.buffer,FALSE)
    dnr.r<-crop(dnr.r,e.tps)
    dnr.buffer<-tps.resample(dnr.r,dem.buffer,FALSE)
    
    plot(sis.buffer); plot(dnr.buffer)
    # work out Julian day and time - ncdf time in hrs since 1/1/1983 0:00
    jul.base<-JDdoy(1,1983)
    hrs<-ncvar_get(ncdf_dnr,"time")
    days<-floor(hrs/24)
    jul.day<-as.numeric(jul.base+days) #JD at 12:00 on ncdf day
    h<-as.numeric(hrs%%24)
    if (h!=hr) {print("WARNING h^= hr!!!")}     
    
    # creates raster template for storing direct and diffuse radiation values with cell values of -9
    direct.r<-dem.buffer*0-9
    diffuse.r<-dem.buffer*0-9
    total.r<-dem.buffer*0-9
    
    #plot(dem.buffer)
    #plot(slope.buffer)
    #plot(aspect.buffer)
    
    # convert values to matrices for use with solar index function
    m.dem<-getValues(dem.buffer,format="matrix") # use.raster function above
    m.slope<-getValues(slope.buffer,format="matrix")
    m.aspect<-getValues(aspect.buffer,format="matrix")
    
    # converts NAs to zeros - NB: buffer of boundary blocks = NA->0
    sel<-which(is.na(m.dem)==T); m.dem[sel]<-0
    sel<-which(is.na(m.slope)==T); m.slope[sel]<-0
    sel<-which(is.na(m.aspect)==T); m.aspect[sel]<-0
    
    # Calculate lat and long of centre of grid
    ll<-OSGBtolatlong(xmin(dem.buffer)+(0.5*(xmax(dem.buffer)-xmin(dem.buffer))) , ymin(dem.buffer)+(0.5*(ymax(dem.buffer)-ymin(dem.buffer))) )
    lat<-as.numeric(ll[2])
    long<-as.numeric(ll[1])
    si<-solarindex(slope=m.slope,aspect=m.aspect,localtime=hr,
                   Lat=lat,Long=long,Julian=jul.day,dtm=m.dem)
    
    si.flat<-solarindex(slope=0,aspect=0,localtime=hr,
                        Lat=lat,Long=long,Julian=jul.day,shadow=F)[1,1]
    
    # Direct normal radiation:
    dnr.m<-getValues(dnr.buffer,format="matrix")
    # Direct radiation: flat
    dir.flat<-dnr.m*si.flat
    # Diffuse radiation flat
    sis.m<-getValues(sis.buffer,format="matrix")
    dif.m<-sis.m-dir.flat
    
    # Calculate matrix coords that define central 'block'
    xmn<-1+(xmin(dem.block)-xmin(dem.buffer))/res(dem.block)[1]
    xmx<-dim(dem.buffer)[1]-(xmax(dem.buffer)-xmax(dem.block))/res(dem.block)[1]
    ymn<-1+(ymin(dem.block)-ymin(dem.buffer))/res(dem.block)[1]
    ymx<-dim(dem.buffer)[2]-(ymax(dem.buffer)-ymax(dem.block))/res(dem.block)[1]
    
    # downscaled direct radiation
    direct<-dnr.m*si
    direct<-direct[xmn:xmx,ymn:ymx] # Extract block
    # downscaled diffuse radiation
    sv<-skyview(m.dem)
    diffuse<-dif.m*sv
    diffuse<-diffuse[xmn:xmx,ymn:ymx] # Extract block
    # add mask back in so that SIs are only produced for land
    mask<-getValues(dem.block,format="matrix")*0
    direct<-direct+mask
    diffuse<-diffuse+mask
    total<-direct+diffuse
    
    # Extract central Block and Convert to rasters
    direct.r<-raster(direct,template=dem.block)
    diffuse.r<-raster(diffuse,template=dem.block)
    total.r<-raster(total,template=dem.block) 
    
  } else {         # IF NIGHTIME         
    direct.r<-dem.block*0
    diffuse.r<-dem.block*0
    total.r<-dem.block*0
  }
  
  # Write files
  if (write.files){
    fileout1<-paste(dir_direct,"direct_",year,"_",sprintf("%02d",month,sep=""),"_",
                    sprintf("%02d",day,sep=""),"_",hr,"h.tif",sep="")
    fileout2<-paste(dir_diffuse,"diffuse_",year,"_",sprintf("%02d",month,sep=""),"_",
                    sprintf("%02d",day,sep=""),"_",hr,"h.tif",sep="")
    fileout3<-paste(dir_total,"total_",year,"_",sprintf("%02d",month,sep=""),"_",
                    sprintf("%02d",day,sep=""),"_",hr,"h.tif",sep="")
    writeRaster(direct.r,file=fileout1,overwrite=T)
    writeRaster(diffuse.r,file=fileout2,overwrite=T)
    writeRaster(total.r,file=fileout3,overwrite=T)
  }
  #output<-list(direct.r,diffuse.r,total.r)
  
  #direct.r;diffuse.r;total.r;dem.block;dem.buffer
  #Create raster stack and print if requested
  if (print.results){
    result.stack<-stack(direct.r,diffuse.r,total.r,dem.block)
    names(result.stack)<-c("direct","diffuse","total","dem block")
    par(mfrow=c(2,2))
    plot(result.stack)  
  } # end if
  
  results<-c(direct.r,diffuse.r,total.r)
  return(results)
} # end function
  
  
  

# CALL FUNCTION
year<-1992
month<-7
day<-1
radiation_daystack(day,month,year)
  
radiation_downscale_stack(day,month,year,hr,demuk,dem.buffer,dem.block,slope.buffer,aspect.buffer,dir_rad,print.results=TRUE,write.files=FALSE)
  
