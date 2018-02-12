#NEW CODE FOR RESAMPLING RASTER using thin plate spline
#library(fields) # required for thin plate spline
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





# Possible use for wind speed as alternative to simple altitude adjustment
# Attempt for different speeds - calaulcate corrections factor for each 100m cell
rain<-raster(month.rain[,,mth],template=r)
rs<-raster(array(0,dim=c(32,53)),template=r) # a black raster dataset at 5 km res of the right size
xy <- data.frame(xyFromCell(rs, 1:ncell(rs))) # extracts xy coordinates of every cell
xy$z<-getValues(sac) # extracts elevetation data at 5km res for same area
v <- getValues(rain) # gets rain values at 5km
# fit a model
tps <- Tps(xy, v)
# use model to predict values for all locations
rain.fine<-interpolate(studyarea, tps,xyOnly=F)
rain.fine<-crop(rain.fine,e2)
plot(rain.fine,main=mths[mth])
# write out file



#Code for resampling raster
library(fields) # required for thin plate spline
par(mfrow=c(2,2))
# get example DEM
dem<-raster("C:/Jonathanmodel/DEM/demsw.asc")
plot(dem,main="High-res DEM")
# Coarsen DEM to illustrate how bilinear interpolation versus thin-plate spline works
dem.coarse<-aggregate(dem,100)
plot(dem.coarse,main="Coarsened DEM")
# Conventional bilinar resample
dem.res<-resample(dem.coarse,dem)
dem.res<-mask(dem.res,dem)
plot(dem.res,main="Bilinear-interpolation")
# Thin plate spline
xy <- data.frame(xyFromCell(dem.coarse,1:ncell(dem.coarse)))
v <- getValues(dem.coarse)
tps <- Tps(xy, v) # fit tps model (Don't worry about warning)
dem.tps<- raster(dem) #create blank raster
# use model to predict values for all locations
dem.tps<- interpolate(dem.tps,tps)
dem.tps<- mask(dem.tps,dem)
plot(dem.tps,main="Thin-plate spline")
