# Calculate and write maps of % land within a set radius of cell
# Input: radius, 5km grid, 100m dem for uk, extend of interest
# Outputs: basic and complete rasters of % land at 5km and 100m res
# Complete rasters estimate % land for cells without historic data from nearest valid cell

#####################################################################
# Minor Functions
#####################################################################

# Calculates number of 100m landcells in each 5km
landcells<-function(dem100m,ref5kmgrid.r,gridonly=TRUE)
{ 
  factor<-res(ref5kmgrid.r)[1] / res(dem100m)[1]
  land5km.r<-aggregate(dem100m,factor,fun=function(x,...)length(na.omit(x))) # all 5km cells
  #plot(land5km.r)
  if (gridonly){land5km.r<-mask(land5km.r,ref5kmgrid.r)} # only those 5km cells in ref5kmgrid.r
  plot(land5km.r)
  return(land5km.r) 
} # end function


# Find nearest coordinates in refx/y to coordinates x/y and return correspoinding ref value
# calculates mean value if several points at equal distance
nearestVal<-function(x,y,refx,refy,refvals)
{ 
  x<-rep(x,length(refx))
  y<-rep(y,length(refy))
  dist<-sqrt( abs(x-refx)^2 + abs(y-refy)^2 ) 
  sel<-which(dist==min(dist))
  if (length(sel)==0){warning(paste("Warning - no min found in nearLref: x=",x," y=",y," sel=",sel," min(dist)=",min(dist),sep=""))}
  if (length(sel)>1){refval<-mean(refvals[sel],na.rm=TRUE)
  } else {refval<-refvals[sel] }
  #print(refval)
  return(refval)
}


#Function: fill.5km.map
# USES FUNCTIONS landcells, nearestVal
# Input:
# ref5kmgrid.r = 5km frame covering same extent as 
# refdata5km.r = existing(incomplete) 5km data raster
# dem100m = 100m resolution land/sea raster (sea = NA, land >0)
# Output:
#   ref5km.filled.r  = filled ref 5km data taking values from nearest 5km cell to fill missing cells that contain 100m land cells

fill.5km.map<-function(refdata5km.r,dem100m,ref5kmgrid.r) 
{
  if(compareRaster(refdata5km.r,ref5kmgrid.r)!=TRUE){warning("!!! 5km ref data and 5km grid do not match !!!")}
  # A. identify missing 5km cells (containing 100m land cells but without values)
  missing.r<-overlay(refdata5km.r,landcells(dem100m,ref5kmgrid.r,FALSE), fun=function(x,y){ifelse(is.na(x) & y>0,0,x)})
  plot(missing.r)
  
  # Select xy and val of  cells without historic cover
  xyvals<-rasterToPoints(missing.r)
  sel<-which(xyvals[,3]==0)
  missingxyv<-xyvals[sel,1:3]
  
  # Find value from nearest cell with historic temp data
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
  # B. Assign missing %land values (0)  to missingxyv values  
  ref5km.filled.r<-rasterize(missingxyv[,1:2], refdata5km.r, missingxyv[,3], fun=max, update=TRUE)
  plot(ref5km.filled.r,main="% land 5km cells - cells missing historic data set to nearest values")
  
  return (ref5km.filled.r)
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

#####################################################################
# Main function for writing % land maps
# Calculate % land within radius of each 100m cell
# Output: i. %land within radius of 5km historic gridcells
#         ii. %land within radius of 100m cells
#         iii.% land for all land 5km cells - cells without historic T data assigned %land of nearest neighbour
#         iv.  % land for all land 100m cells - cells without historic T data assigned %land of nearest neighbour (=Lref)
#####################################################################

percent_land_maps<-function(dembuf,grid5km.r,e.dem,radius,dir_percland) {
    # If rasters not already suitably cropped
    #radius<-10000
    #xmn<-xmin(e.dem)-radius
    #xmx<-xmax(e.dem)+radius
    #ymn<-ymin(e.dem)-radius
    #ymx<-ymax(e.dem)+radius
    #e<-extent(xmn,xmx,ymn,ymx)
    #dem.buf<-crop(demuk,e)
    #grid5km.r<-crop(grid5km.r,e.dem)
    
    # Calculate % land for each 100m cell in dem.buf
    coast.r<-percent_land(dembuf,radius) 
    coast.r<-crop(coast.r,e.dem) # includes values for many 100m cells that are sea
    
    #  i. Write 5km map from central 100m cells (before masking sea cells)
    vals<-rasterToPoints(grid5km.r)
    # extract %land from 100m cell matching 5km xy coords - mean calculated if several cells extracted
    coast5km<-raster::extract(coast.r,vals[,1:2]) 
    coast5km.r<-rasterize(vals[,1:2],grid5km.r,coast5km,fun=mean,background=NA)
    plot(coast5km.r,main="% land 5km grid")
    #out.file<-paste(dir_percland,"percent_land",radius/1000,"km_5km_histgrid.tif",sep="")
    #writeRaster(coast5km.r,file=out.file,overwrite=T)
    
    #  Set sea to NA and output 100m version
    #coast100m.r<-mask(coast.r,crop(dembuf,e.dem))
    #plot(coast100m.r,main=paste("% land in a ",radius/1000,"km radius of all 100m cells", sep=""))
    #out.file<-paste(dir_percland,"landpercent_",radius/1000,"km_100mcells.tif",sep="")
    #writeRaster(coast100m.r,file=out.file,overwrite=T)
    
    # ii. Convert 5km grid back to 100m cell raster
    coast5km100m.r<-disaggregate(coast5km.r,50)
    coast5km100m.r<-mask(coast5km100m.r,crop(dembuf,e.dem)) 
    # ...but is missing some 100m land cells where no overlapping 5km cell
    plot(coast5km100m.r,main="% land reference values for 100m cells - missing values")
    
   # iii. Fill 5km cells without values but containing 100m land cells
    plandref5km.r<-fill.5km.map(coast5km.r,crop(dembuf,e.dem),grid5km.r) 
    plot(plandref5km.r,main="% land 5km cells - cells missing historic data set to nearest values")
    out.file<-paste(dir_percland,"percent_land_",radius/1000,"km_lref_5kmgrid.tif",sep="")
    writeRaster(plandref5km.r,file=out.file,overwrite=T)
    
    # iv. Convert to 100m cells and mask with sea
    plandref100m.r<-disaggregate(plandref5km.r,50)
    plandref100m.r<-mask(plandref100m.r,crop(dembuf,e.dem))
    plot(plandref100m.r,main="% land ref values")
    
    out.file<-paste(dir_percland,"percent_land_",radius/1000,"km_lref_100mgrid.tif",sep="")
    writeRaster(plandref100m.r,file=out.file,overwrite=T)

} # end function



