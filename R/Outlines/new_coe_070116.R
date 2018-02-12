#  NEW INTERPOLATION METHOD FOR CAL
# input:
# dfi = dataset with missing values (effective Cloud albedo with missing night time values)
# output: dataset with missing values imputed (effective Cloud albedo with no missing night time values)
imputation<-function(dfi)
{
  # Ensure values lie between zero and one
  sel<-which(dfi<0.001); dfi[sel]<-0.001
  sel<-which(dfi>0.999); dfi[sel]<-0.999
  dfi<-log(dfi/(1-dfi))
  x<-c(1:length(dfi))
  # use spline to interpolate
  sp<-spline(x,dfi,n=length(x))
  dfi<-1/(1+exp(-1*sp$y))
  dfi
}

#NEW CODE FOR RESAMPLING RASTER
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
dem.tps<-mask(dem.tps,dem)
plot(dem.tps,main="Thin-plate spline")
