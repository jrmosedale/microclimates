§
#NEW CODE FOR RESAMPLING RASTER using thin plate spline
#library(fields) # required for thin plate spline
tps.resample<-function(input.r,output.r,msk=TRUE){
  xy <- data.frame(xyFromCell(input.r,1:ncell(input.r)))
  v <- getValues(input.r)
  tps <- Tps(xy, v) # fit tps model (Don't worry about warning)
  result<- raster(output.r) # create blank raster at resolution of output.r
  
  # use model to predict values for all locations
  result<- interpolate(result,tps)
  if (msk==TRUE) result<-mask(result,output.r)
  #plot(result,main="Thin-plate spline - min")
  
  return(result)
}# end function
  
splint.resample<-function(input.r,output.r,msk=TRUE){
  input.r<-sis.r
  output.r<-dem.buffer
  xy <- data.frame(xyFromCell(input.r,1:ncell(input.r)))
  v <- getValues(input.r)
  z<-sreg(xy$x, xy$y, v)
  
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




#### TPS for temperature at 5km
# 5km elevation data 
elev5km.r<-raster(extent(grid5km.r),crs=crs(grid5km.r),res=res(grid5km.r))  # empty raster of grid5km.r (historic data cells)
elev5km.r<-aggregate(dem,fact=50,fun=sum) # = number of 100m land cells in each 5km cell
plot(elev5km.r,main="5km elevation")
#########################
# TPS for whole dem area
ptm <- proc.time()
# uses t5km.r
xy <- data.frame(xyFromCell(land5km.r, 1:ncell(land5km.r))) # extracts xy coordinates of every cell
#xy$z<-getValues(elev5km.r)# extracts elevetation data at 5km res for same area
t <- getValues(t5km.r) # gets t values at 5km
tps <- Tps(xy, t)
#plot(tps) # info on model
print("finished model")
print(proc.time()-ptm)

# apply model 
t100.r<-interpolate(dem, tps,xyOnly=T)
print("Finished application")
print(proc.time()-ptm)
#########################
# TPS for one block
ptm <- proc.time()
# uses t5km.r
xy <- data.frame(xyFromCell(land5km.r, 1:ncell(land5km.r))) # extracts xy coordinates of every cell
#xy$z<-getValues(elev5km.r)# extracts elevetation data at 5km res for same area
t <- getValues(t5km.r) # gets t values at 5km
tps <- Tps(xy, t)
#plot(tps) # info on model
print("finished model")
print(proc.time()-ptm)

# apply model 
t100.block<-interpolate(dem.block, tps,xyOnly=T)
print("Finished application")
print(proc.time()-ptm)
plot(mask(t100.block,dem.block),main="tps block")
#########################
# Bilinear resampling - grid effect evident
ptm <- proc.time()
t100.block<-resample(t5km.r,dem.block)
print("Finished application")
print(proc.time()-ptm)
plot(mask(t100.block,dem.block),main="bilinear")


# Compare across several gridcells
ptm <- proc.time()
ukcpcell<-c(916:919,924:932,934:935)
map.r<-raster()
xy <- data.frame(xyFromCell(land5km.r, 1:ncell(land5km.r))) # extracts xy coordinates of every cell
t <- getValues(t5km.r) # gets t values at 5km
tps <- Tps(xy, t)
for (n in(1:length(ukcpcell)))
{ 
  x<-landcells[ukcpcell[n],1]
  y<-landcells[ukcpcell[n],2]
  e.block<-extent(x-2500,x+2500,y-2500,y+2500)
  dem.block<-crop(demuk,e.block)# create day stack of hourly data for whole of area

  t.r<-interpolate(dem.block, tps,xyOnly=T)
  if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  
}
print(proc.time()-ptm)

plot(mask(map.r,crop(dem,map.r)),main="tps")


#################
# Compare tps models from different hours of same day
for (n in c(1,7,13,19,24)){
  t5km.r<-raster(t5km.s,layer=n)
  xy <- data.frame(xyFromCell(land5km.r, 1:ncell(land5km.r))) # extracts xy coordinates of every cell
  t <- getValues(t5km.r) # gets t values at 5km
  tps <- Tps(xy, t)
  print(n)
  print(tps)
}



### Code for resampling raster
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
dem.tps<-mask(dem.tps,dem)
plot(dem.tps,main="Thin-plate spline")
