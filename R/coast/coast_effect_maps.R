# Calculate and write maps of % land within a set radius of cell
# Input: radius, 5km grid, 100m dem for uk, extend of interest
# Outputs: basic and complete rasters of % land at 5km and 100m res
# Complete rasters estimate % land for cells without historic data from nearest valid cell
# Prev Programs: 


library(raster)
library(rgdal)
dir_percland<-"~/Documents/Exeter/Data2015/CoastEffect/percland/"
dir_grids<-"~/Documents/Exeter/Data2015/Templates/"
dir_lsratio<-"~/Documents/Exeter/Data2015/CoastEffect/lsratio/"
#dir_grids<-"C:/Data2015/Templates/"
#dir_percland<-"C:/Data2015/CoastEffect/percland/"

#####################################################################
# Calculate number of 100m landcells in 5km reference grid cells
landcells<-function(inv.lsratio,gridmask.r)
{ 
  land5km.r<-aggregate(inv.lsratio,50,fun=function(x,...)length(na.omit(x))) # number of 100m landcells in each 5km
  #land5km.r<-mask(land5km.r,gridmask.r)
  #plot(lref5km.r)
  return(land5km.r) # x,y,Lref
} # end function

# Calculate % land within radius (m) of each cell (res in m)
# Output: raster of  % land
percent_land<-function(map.r,radius)
{
  # Create raster where Sea=0, Land=1 no NA
  ls.vals<-getValues(map.r)
  ls.vals<-ifelse(is.na(ls.vals),0,1)
  landsea.r<-setValues(map.r,ls.vals)
  
  # Create filter & Calculate % of cells = land
  f<-focalWeight(map.r,radius,type=c("circle"))
  result.r<-focal(landsea.r,f,pad=TRUE)
  #plot(land_perc)
  return(result.r)
}

# Find nearest coordinates in refx/y to coordinates x/y and return correspoinding ref value
# calculates mean value if several points at equal distance
nearestVal<-function(x,y,refx,refy,refvals)
{ 
  x<-rep(x,length(refx))
  y<-rep(y,length(refy))
  dist<-sqrt( abs(x-refx)^2 + abs(y-refy)^2 ) 
  sel<-which(dist==min(dist))
  if (length(sel)==0){warning(paste("Warning - no min found in nearLref: x=",x," y=",y," sel=",sel," min(dist)=",min(dist),sep=""))}
  if (length(sel)>1){Lref<-mean(refvals[sel],na.rm=TRUE)
  } else {Lref<-refvals[sel] }
  #print(Lref)
  return(Lref)
}

#####################################################################
latlong <- "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"

#demuk<-raster("C:/Data2015/DEM100/demoriginal.tif", crs=(ukgrid))
demuk<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif", crs=(ukgrid))

# Define sw dem of interest (no buffer)
#e.dem<-extent(c(70000,420000,0,180000 )) # includes scilly isles
e.dem<-extent(c( 130000,400000,10000,180000 )) # excludes scilly isles
dem<-crop(demuk,e.dem)

e.dem<-extent(dem)

in.file<-paste(dir_grids,"ukhistmask.grd",sep="")
print(in.file)
gridmask.r<-raster(in.file) #  1=valid cell, NA = sea or not data
gridmask.r<-crop(gridmask.r,e.dem) # IMPORTANT: crop to SAME geographical extent of DEM and inv.lsratio raster
#####################################################################

# Calculate % land within 10 and 30km of each 100m cell
# Output: raster of  % land
radius<-10000
xmn<-xmin(dem)-radius
xmx<-xmax(dem)+radius
ymn<-ymin(dem)-radius
ymx<-ymax(dem)+radius
e<-extent(xmn,xmx,ymn,ymx)
dem_buf<-crop(demuk,e)
# plot(dem_buf)

coast.r<-percent_land(dem_buf,radius)
coast.r<-crop(coast.r,dem)

#  1. Calc 5km map from central 100m cells (before masking sea cells)
vals<-rasterToPoints(gridmask.r)
coast5km<-raster::extract(coast.r,vals[,1:2]) 
coast5km.r<-rasterize(vals[,1:2],gridmask.r,coast5km,fun=mean,background=NA)
crop(coast5km.r,dem)
plot(coast5km.r)
out.file<-paste(dir_percland,"landin_",radius/1000,"km_5km_histgrid.r",sep="")
save(coast5km.r,file=out.file)

# 2. Set sea to NA and output 100m version
coast100m.r<-mask(coast.r,dem)
plot(coast100m.r,main=paste("% land in a ",radius/1000,"km radius", sep=""))
out.file<-paste(dir_percland,"landin_",radius/1000,"km_100mgrid.r",sep="")
save(coast100m.r,file=out.file)

# Convert 5km grid back to 100m cell raster
coast5km100m.r<-disaggregate(coast5km.r,50)
coast5km100m.r<-mask(coast5km100m.r,dem) 
# ...but is missing some 100m land cells where no overlapping 5km cell
#plot(coast5km100m.r,main="5km coast effect at 100m cells")

#  5km cells without historic data but containing 100m land cells set to 0 %land
landcells<-landcells(dem,gridmask.r)
plot(landcells,main="Number of 100m landcells in 5km cells")
missing.r<-overlay(coast5km.r,landcells, fun=function(x,y){ifelse(is.na(x) & y>0,0,x)})
plot(missing.r, main="Cells without historic data=0")

# For each missing cell set %land to 'nearest' 5km cell for which historic data available
# and which will therefore provide temperature data
# Select xy and val of  cells without historic cover
xyvals<-rasterToPoints(missing.r)
sel<-which(xyvals[,3]==0)
missingxyv<-xyvals[sel,1:3]

# Assign value from nearest cell with Lref
sel<-which(!is.na(xyvals[,3]) & xyvals[,3]!=0) # selects land cells with data
refvals<-xyvals[sel,1:3]

# set vectors
x<-missingxyv[,1]
y<-missingxyv[,2]
refx<-refvals[,1]
refy<-refvals[,2]
refvals<-refvals[,3] 

for (i in 1: length(x)){
  missingxyv[i,3]<-nearestVal(x[i],y[i],refx,refy,refvals)
}

# Update missing values (0) in lref5km to missingxyv values  
complete5km.r<-rasterize(missingxyv[,1:2], coast5km.r, missingxyv[,3], fun=max, update=TRUE)
plot(complete5km.r,main="% Land 5km cells - cells missing historic data set to nearest values")

out.file<-paste(dir_percland,"landin_",radius/1000,"km_complete_5km.r",sep="")
save(complete100m.r,file=out.file)

# Convert to 100m cells and mask with sea
complete100m.r<-disaggregate(complete5km.r,50)
complete100m.r<-mask(complete100m.r,dem)
plot(complete100m.r)

out.file<-paste(dir_percland,"landin_",radius/1000,"km_complete_100m.r",sep="")
save(complete100m.r,file=out.file)



