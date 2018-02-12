# Calculates number of 100m landcells in each 5km
landcells<-function(dem100m,ref5kmgrid.r,gridonly=TRUE)
{ 
  factor<-res(ref5kmgrid.r)[1] / res(dem100m)[1]
  land5km.r<-aggregate(dem100m,factor,fun=function(x,...)length(na.omit(x))) # all 5km cells
  #plot(land5km.r)
  if (gridonly){land5km.r<-mask(land5km.r,ref5kmgrid.r)} # only those 5km cells in ref5kmgrid.r
  #plot(land5km.r,main="Number of 100m cells")
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
  #plot(missing.r,main="Missing cells (grey)")
  
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
  #plot(ref5km.filled.r,main="% land 5km cells - cells missing historic data set to nearest values")
  
  return (ref5km.filled.r)
} # end function


# Write elevation difference maps Eregf (5km cell) - Ecell (100m land cell)
# USES  fill.5km.map
# Compare percentland_map_function

elevation.dif.map<-function(dem,uk5kmgrid.r,dir_terrain){
    refgrid.r<-crop(uk5kmgrid.r,dem)
    
    # A. Not possible to calculate 5km map from central 100m cells as many of these are sea cells for coastal regions
    elev5km.r<-mask(aggregate(dem,fact=c(50,50),fun=mean),refgrid.r) # calculate mean elevation and remove sea cells according to refgrid.r
    plot(elev5km.r,main="Elev mean 5km")
    elev100m.r<-disaggregate(elev5km.r,c(50,50))
    plot(elev100m.r,main="Elev 100m from 5km")
    
    # B. Fill 5km cells without values but containing 100m land cells
    complete5km.r<-fill.5km.map(elev5km.r,dem,refgrid.r) 
    
    plot(complete5km.r,main="% land 5km cells - cells missing historic data set to nearest values")
    
    # C. Write files filled 5km maps
    out.file<-paste(dir_terrain,"ref_elevation_5km.tif",sep="")
    print(out.file)
    writeRaster(complete5km.r,file=out.file,overwrite=T)
    
    # i.v Convert to 100m cells, mask and write file 
    complete100m.r<-disaggregate(complete5km.r,50)
    complete100m.r<-mask(complete100m.r,dem)
    plot(complete100m.r)
    out.file<-paste(dir_terrain,"ref_elevation_100m.tif",sep="")
    print(out.file)
    writeRaster(complete100m.r,file=out.file,overwrite=T)
    
    # Calculate elevation difference and write file
    ediff.r<-complete100m.r-dem
    plot(ediff.r)
    out.file<-paste(dir_terrain,"eref-edem_100m.tif",sep="")
    print(out.file)
    writeRaster(ediff.r,file=out.file,overwrite=T)
} # end function




