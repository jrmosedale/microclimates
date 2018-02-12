x<-c(-10:5)
x<-c(10:20)

for (y in 145:215){
if (y<210 & y>150) finalt<-x+((x/3.5)*1/(1+exp(6-y)) ) else finalt<-x
print(finalt)
}  
  
for (y in 329:360){
  if (y<30 | y>330) finalt<-x-((x/5)*1/(1+exp(6-y)) ) else finalt<-x
  print(finalt)
}    
  
  

for (x in 0:20){
 #print(1/(1+exp(6-y)))
  #print(x-((x/5)*1/(1+exp(6-y)) ))
  z<-x-(4*1/(1+exp(6-y)) )
  print(z)
     }

ifelse((y<30 | y>330),x-((x/5)*1/(1+exp(6-y)) ),x)



# Define projections
latlong = "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"
new.crs<-"+init=epsg:3857"
dir_terrain<-paste(in.root,"Terrain/",sep="")
dir.basinmap<-paste(root,"basins/",sep="")

dir_in<-paste(root,"Outputs/",sep="")
dir_rasters<-paste(root,"shinydata/rasters/",sep="")
dem<-raster(paste(dir_dem,"dem.tif",sep=""))
slope.sw<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem)
aspect.sw<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem)

t<-raster(paste(dir_in,"fl_tmean_2010.tif",sep=""))

elev<-raster(paste(dir_rasters,"elevation_cornwall.tif",sep=""))

slope<-raster(paste(dir_rasters,"slope_cornwall.tif",sep=""))
aspect<-raster(paste(dir_rasters,"aspect_cornwall.tif",sep=""))

cornwall.r<-raster(paste(dir_counties,"cornwallmask.tif",sep=""))
crs(cornwall.r)<-new.crs
twi<- raster(paste(dir.basinmap,"topidx.tif",sep=""))
altdif<-raster(paste(dir.basinmap,"altdif.tif",sep="")) # for cold air drainage
elevdif<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # elevation effect
crs(twi)<-ukgrid 
crs(altdif)<-ukgrid 
crs(elevdif)<-ukgrid 
twi<-projectRaster(twi,to=cornwall.r)
altdif<-projectRaster(altdif,to=cornwall.r)
elevdif<-projectRaster(elevdif,to=cornwall.r)
lapserate<- 7/1000 # C per 1000m altitude diference 
valleyef<-(altdif*lapserate)*twi ; # plot(valleyef)

t<-raster(paste(dir_rasters,"tmean_gs_2000_cornwall.tif",sep=""))
compareRaster(t,slope,aspect)
tmax<-raster(paste(dir_rasters,"tmax_year_min_1996-2011_cornwall.tif",sep=""))
tmin<-raster(paste(dir_rasters,"tmin_year_mean_1996-2011_cornwall.tif",sep=""))
crs(tmax)<-new.crs; crs(tmin)<-new.crs
plot(tmax)
plot(tmin)
compareRaster(tmax,tmin,slope,aspect)

cloud<-runif(1)
print(cloud)
# correct for aspect & slope
final.tmax<-overlay(tmax,aspect,slope,fun=function(x,y,z) {x+(cloud*( cos(y)*z*(1/(1+exp(0.7)) )))  }) 
plot(final.tmax)
# Increase tmax on south facing slopes
#final.tmax<-overlay(tmax,aspect,slope,fun=function(x,y,z) {ifelse((y<210 & y>150),x+((x/3)*1/(1+exp(6-z)) ),x)  }) # south facing slopes
#final.tmax<-overlay(final.tmax,aspect,slope,fun=function(x,y,z) {ifelse((y<30 | y>330),x-(3*1/(1+exp(6-y)) ),x)  }) # north facing slopes
plot(final.tmax-tmax,main="MAx dif s & n face")
plot(crop(final.tmax,extent(-560000,-520000,6520000,6560000)) )

#NEW
inversion<-runif(1)
print(inversion)
final.tmax<-overlay(tmax,aspect,slope,fun=function(x,y,z) {x+(cloud*(cos(y)*z*(1/(1+exp(0.7)))))  }) # south facing slopes


# correct for elevation difference
final.tmax<-final.tmax+(elevdif*lapserate)  # ; plot(day.tmax,main="day.max")
plot(final.tmax); plot(final.tmax-tmax,main="DIF")

final.tmin<-overlay(tmin,valleyef,fun=function(x,y) { x-(inversion*y)   })
plot(final.tmin-tmin,main="after basin")
plot(crop(final.tmin,extent(-560000,-520000,6520000,6560000)) )

final.tmin<-final.tmin+(elevdif*lapserate)
plot(final.tmin)


test<-as.str

sfacing<-calc(aspect,fun=function(x) {ifelse((x<210 & x>150),1,0)  })
rasterNA(sfacing)+rasternil(sfacing)-ncell(sfacing)


rasterNA<-function(r){
  vals<-getValues(r)
  length(which(is.na(vals)))
}

rasternil<-function(r){
  vals<-getValues(r)
  length(which(vals==0))
}


for (n in 1:10){
  print(runif(1))
}