
# TWO stage altitude correction to 5km then 100m using tpswith altitude
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
  # Stage 2: convert to 100m resolution raster OSGB grid reference using tps resampling
  #############
  # Convert to raster (original lat long format and resolution - CHECK LAT/LON MAX/MIN
  uwind.r<-raster(uwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
  vwind.r<-raster(vwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
  # Reproject in OSGB projection
  crs(uwind.r)<-latlong
  crs(vwind.r)<-latlong
  u_osgb<-projectRaster(uwind.r,crs=ukgrid)
  v_osgb<-projectRaster(vwind.r,crs=ukgrid)
  
  # Coarsen DEM to 5km resolution
  dem.5km<-aggregate(dem,50)
  plot(dem.5km,main="Coarsened DEM")
  # tps resample wind to 5km then apply altitude correction
  u_5km<-tps.resample(u_osgb,dem.5km)
  v_5km<-tps.resample(v_osgb,dem.5km)

  # correct for altitude
  ustr<-sqrt(u_5km^2) # wind strength
  vstr<-sqrt(v_5km^2) # wind strength
  # adjust to altitude
  u.adj<-overlay(ustr, dem.5km, fun=function(x,y){x*((-0.000000108025)*y^2+0.000408692*y+0.956139)})
  v.adj<-overlay(vstr, dem.5km, fun=function(x,y){x*((-0.000000108025)*y^2+0.000408692*y+0.956139)})
  # adjust to 1m above ground level
  u.adj<-u.adj*0.373686439
  v.adj<-v.adj*0.373686439
  # correct for direction
  u.adj<-calc(u.adj,fun=function(x){ifelse(x>0,x*1,x*-1)})
  v.adj<-calc(v.adj,fun=function(x){ifelse(x>0,x*1,x*-1)})
  
  # Calculate Wind Direction
  dir.5km<-overlay(u.adj,v.adj,fun=function(x,y){(180/pi)*(atan2(x,y))}) # NB this is direction in which wind blows to
  dir.5km<-calc(dir.5km,fun=function(x){ifelse(x<=180,x+180,x-180)} ) # NB this direction from which wind originates (360 deg)
  # Calculates Wind Strength from u and v components
  wstr.5km<-overlay(u.adj,v.adj, fun=function(x,y){return(sqrt(x^2+y^2))} )
  plot(dir.5km,main="dir 5km")
  plot(wstr.5km,main="str 5km")
   
  # Strength tps to 100m including altitude as variable
  wstr.100m<-tps.xyz(dem.5km,wstr.5km,dem) # xy and altitude
  wdir.100m<-tps.resample(dir.5km,dem) # xy only
  plot(wstr.100m,main="2-stage str")
  plot(wdir.100m,main="2-stage dir")
  #plot(crop(wstr.100m,dem.block))
  
  #############
  # Stage 5: Format outputs
  #############
  # Calculate inverse wind strength
  invstr.m<-1/(sqrt(str.m+1))
  # Convert to raster and save tif files 
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
  
  
  
  
  # FUNCTION - tps model and apply
  tps.xyz<-function(coarse.r,zvalues.r,fine.r){
    xy <- data.frame(xyFromCell(coarse.r, 1:ncell(coarse.r))) # extracts xy coordinates of every cell
    xy$z<-getValues(coarse.r) # extracts elevetation data at 5km res for same area
    v <- getValues(zvalues.r) # gets rain values at 5km
    tps <- Tps(xy, v)# fit a model
    # use model to predict values for 100m cells
    output.r<-interpolate(fine.r, tps,xyOnly=F)
    #plot(output.r,main="100m tps wind str")
    return(output.r)
  } # end function

tps.resample<-function(input.r,output.r){
  xy <- data.frame(xyFromCell(input.r,1:ncell(input.r)))
  v <- getValues(input.r)
  tps <- Tps(xy, v) # fit tps model (Don't worry about warning)
  
  # use model to predict values for all locations
  result<- raster(output.r) # create blank raster at resolution of output.r
  result<- interpolate(result,tps)
  result<-mask(result,output.r)
  #plot(result,main="Thin-plate spline")
  
  return(result)
}# end function
