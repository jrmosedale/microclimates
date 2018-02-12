####################################################
# # # This Function downscales the wind to 100m cells for a specific hour
####################################################
# Specify hour,day, month and year for which data are required
#year=2010;month=6;day<-10;hr<-11
# Stages:
# (1) get wind values for a given hour, day, month and year
# (2) convert to 100m resolution raster OSGB grid reference
# (3) adjust based on altitude
# (4) adjust based on shelter coefficient

wind.downscale<-function(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter,interval=10,print.results=FALSE,write.files=FALSE)
{
  # Can be quite slow. Allows you to keep tabs on progress by printing hour, day, month & year
  tp<-paste("year=",year," month=",month," day=",day," hour=",hr,sep="")
  print(tp)
  dir_windstrength<-paste(dir_wind,"strength/",sep="")
  dir_winddirection<- paste(dir_wind,"direction/",sep="")
  dir_windinvstr<- paste(dir_wind,"invstr/",sep="")
  
  #############
  # Stage 1: get wind values for a given day month and year
  #############
  # As original data are 4x daily, but data are required for each hour,
  # this bit reads in the data for the periods immediatly before after for which there are data and calculates
  # weighted mean
  av1<-array.val(hr,day,month,year)
  av2<-av1+1
  rem<-hr/6-floor(hr/6)
  uwind1<-wind_u[,,av1]
  uwind2<-wind_u[,,av2]
  vwind1<-wind_v[,,av1]
  vwind2<-wind_v[,,av2]
  uwind<-(1-rem)*uwind1+rem*uwind2
  vwind<-(1-rem)*vwind1+rem*vwind2
  
  #############
  # Stage 2: convert to 100m resolution raster OSGB grid reference using e.block and dem.block
  #############
  # Convert to raster (original lat long format and resolution - CHECK LAT/LON MAX/MIN
  uwind.r<-raster(uwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
  vwind.r<-raster(vwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
  # Reproject in OSGB projection
  crs(uwind.r)<-latlong
  crs(vwind.r)<-latlong
  u_osgb<-projectRaster(uwind.r,crs=ukgrid)
  v_osgb<-projectRaster(vwind.r,crs=ukgrid)
  # Trim to desired dem block and resample the data at a 100m resolution. Interpolation is set at default simple bilinear
  u_100<-resample(u_osgb,dem.block)
  v_100<-resample(v_osgb,dem.block)
  
  #############
  # Stage 3: adjust based on altitude of terrain
  # NB Height adjustment based on wind spped values at different pressures downloaded  Earth System Research Lab
  # Typical heights at different pressures calculated from Allen et al 1998 http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric pressure (p)
  # Quadratic function fitted - NB this works well for heights up to ~1800m. IT won't work above ~2000m
  # Function was first derived by comparing values at different pressures (heights) over the course of a year (2014)
  #############
  # adjust wind speeds by height of dem
  # convert to matrices
  uwind.m<-getValues(u_100,format="matrix")
  vwind.m<-getValues(v_100,format="matrix")
  dem.m<-getValues(dem.block,format="matrix")
  # adjust wind by height
  ustr<-sqrt(uwind.m^2) # wind strength
  vstr<-sqrt(vwind.m^2) # wind strength
  udir<-ifelse(uwind.m>0,1,-1) # positive or negative
  vdir<-ifelse(vwind.m>0,1,-1) # positive or negative
  u.adj<-ustr*((-0.000000108025)*dem.m^2+0.000408692*dem.m+0.956139) # NB don't worry about warnings. Calculation assigns NAs to the sea
  v.adj<-vstr*((-0.000000108025)*dem.m^2+0.000408692*dem.m+0.956139)  # NB don't worry about warnings. Calculation assigns NAs to the sea
  # adjust values to correspond to wind speed 1m above theground
  # rescaling factor first derived by comparing values ot Culdrose wind data using wind_downscale3
  # however, in line wiht what you'd expect from: http://www.fao.org/docrep/x0490e/x0490e07.htm#wind profile relationship
  u.adj<-u.adj*0.373686439
  v.adj<-v.adj*0.373686439
  u.adj<-u.adj*udir
  v.adj<-v.adj*vdir
  u.r<-raster(u.adj,template=dem.block)
  v.r<-raster(v.adj,template=dem.block)
  
  # calculate wind direction
  direction = (180/pi)*(atan2(u.adj,v.adj))  # NB this is direction in which wind blows to
  direction<-ifelse(direction<=180,direction+180,direction-180) # NB this direction from which wind originates (360 deg)
  dir.r<-raster(direction,template=dem.block)
  
  # Some code here for plotting altitude adjusted values, currently commmented out
  #par(mfrow=c(2,2))
  #plot(u.r,main="altitude adjusted wind u")
  #plot(v.r,main="altitude adjusted wind v")
  #plot(dem.block,main="dem")
  #plot(dir.r,main="wind direction - stg3")
  
  #############
  # Stage 4: height adjustments done using shelter coefficient maps based on topography and wind direction
  #############
  
  # Create matrixes from rasters of wind vectors and direction - ORIENTATION???
  u.m<-getValues(u.r,format="matrix")
  v.m<-getValues(v.r,format="matrix")
  dir.m<-getValues(dir.r,format="matrix")
  
  # Uses Rounded wind direction of each cell to select correct shelter coefficient
  dir.m<-round(dir.m/interval)*interval # round to interval used for shelter maps
  dir.m<-ifelse(dir.m==0,360,dir.m) # converts 0 to 360 degree direction
  
  # Create windstrength matrix for storing values - original: matrix(NA,nrow=1800,ncol=3000)
  str.m<-matrix(NA,nrow=(nrow(dem.m)),ncol=(ncol(dem.m)) )# matrix for storing all values
  
  # Calculates wind strength from u and v components
  str.m<-sqrt(u.m^2+v.m^2)
  #set str.m to NA if sea (on basis of dem)
  str.m[which(is.na(dem.m))]<-NA
  
  # Applies shelter coefficient to calculate wind strength if land cell (else NA)
  mxrws<-nrow(dem.m)
  mxcls<-ncol(dem.m)
  for (rws in 1:mxrws) {
    for (cls in 1:mxcls) {    
      str.m[rws,cls]<-str.m[rws,cls]*shelter[(dir.m[rws,cls]/interval),rws,cls] 
    }
  }
  
  # Calculate inverse wind strength
  invstr.m<-1/(sqrt(str.m+1))
  
  # Convert to raster and save tif files - QUESTION - SHOULD RASTERS BE MASKED WITH DEM.BLOCK ???
  wstr.r<-raster(str.m,template=dem.block)
  wdir.r<-raster(dir.m,template=dem.block)
  invwstr.r<-raster(invstr.m,template=dem.block)
  
  if (write.files==TRUE){
    fileout.1<-paste(dir_windstrength,"strength_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    fileout.2<-paste(dir_winddirection,"direction_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    fileout.3<-paste(dir_windinvstr,"invstr_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    #print(fileout.1); print(fileout.2); print(fileout.3)
    writeRaster(wstr.r,file=fileout.1,overwrite=TRUE)
    writeRaster(wdir.r,file=fileout.2,overwrite=TRUE)
    writeRaster(invwstr.r,file=fileout.3,overwrite=TRUE)
  }
  
  #Create raster stack and print if requested
  if (print.results==TRUE){
    shelterblock<-raster(shelter[(mean(dir.m,na.rm=TRUE)/interval),,],template=dem.block)
    result.stack<-stack(wstr.r,wdir.r,invwstr.r,dem.block,shelterblock)
    names(result.stack)<-c("wind strength","wind direction","inv wind str","dem block","shelter block")
    par(mfrow=c(2,3))
    plot(result.stack)  
  } # end if
  
  wind.results<-c(wdir.r,wstr.r,invwstr.r)
  return(wind.results)
  
} # end function wind.downscale