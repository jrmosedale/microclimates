##########################################################################################
# PROCESS INPUT DATA
# Calculate all DEM derived maps for use in modelling
# Outputs constant across time scale
# New Version 3 - uses inv.dist maps calculated using circle sector not line
##########################################################################################
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

##########################################################################################
# MAP FUNCTIONS USED BY ABOVE
##########################################################################################
# WIND SHELTERMAPS
# Works out the angle to the horizon in a specified direction (used to calculate the shelter coefficient)
# Inputs:
# dtm = a digital eleveation model stored as a matrix
# Outputs:
# 360 files - map of shelter coefs for each wind direction - Save as RASTERS .tif format

# NB the rotation of the digital elevetation data is important. This is designed to be used for a matrix
# extracted from a raster (see raster package) as follows: my.matrix<-getValues(my.raster,format="matrix")
horizonangle <- function(dtm,azimuth,res=100)
{
  dtm<-(dtm*5)/res
  azimuth<-azimuth-90
  azi <- azimuth * (pi/180)
  horizon <- array(0,dim(dtm))
  dtm3 <- array(0,dim(dtm)+200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)] <- dtm
  for (step in 1:10) {
    horizon[1:x,1:y] <- pmax(horizon[1:x,1:y], (dtm3[(101+sin(azi)*step^2):(x+100+sin(azi)*step^2),(101+cos(azi)*step^2):(y+100+cos(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(5*step^2))
  }
  horizon
}

windindex <- function(dtm,direction)
{
  index <- 1 - atan(0.17 * 100 * horizonangle(dtm,direction))/1.65
  index
}

########################################################################################
# Main FUNCTION - wind_sheltermaps
########################################################################################

# Go through in 10 km blocks and adjust by shelter coefficient
# note, however that actually it actually selects 30km x 30 km area to allow for sheltering effects
# that operate outside the area of each block. The 10km x 10km centre of the block is then selected
# Programme could probably be speeded up without major loss of accuracy by setting buffer to ~5km instead of 10km
# NB chopping into 10km blocks is necessary, as the function for calculating
# the shelter coefficient assumes a single wind direction, a fairly safe assumption over 10km, but not over entire study region

wind_sheltermaps<-function(demuk,e.dem,interval,dir_shelter){
  # 1. Define buffered dem (10km buffer) and derived matrix cells
  
  # ASSUMES dem already cropped to match 10km cells 
  # 10km buffer required to produce maps for e.dem
  buffer<-10000
  e.map<-extent(xmin(e.dem)-buffer,xmax(e.dem)+buffer,ymin(e.dem)-buffer,ymax(e.dem)+buffer)
  dem.map<-crop(demuk,e.map)
  plot(dem.map,main="DEM-map")
  # Convert dem matrix
  dem.m<-getValues(dem.map,format="matrix")
  # Define 10km blocks -2 excludes outside cells for which buffer cannot be calculated
  mxrws<-nrow(dem.m)/100-2
  mxcls<-ncol(dem.m)/100-2
  b.cells<-buffer/100
  
  # 2. For every wind direction (1-360 degrees in steps of 5) calculate a shelter coefficient map
  #interval<-10
  for (d in seq(0,360,interval))
  {
    # creates 2D matrix for storing wind coeff output 
    wcoef<-matrix(NA,nrow=(nrow(dem.m)),ncol=(ncol(dem.m)) )# matrix for storing output values
    sea<-FALSE
    for (rws in 1:mxrws) 
    {
      for (cls in 1:mxcls)
      { 
        xmn<-rws*100+1-100
        ymn<-cls*100+1-100
        xmx=xmn+100+(2*b.cells)-1 
        ymx=ymn+100+(2*b.cells)-1 
        b.dem<-dem.m[xmn:xmx,ymn:ymx]
        # Check if dem block is all sea
        sel<-which(is.na(b.dem)==T)
        if (length(sel)==length(b.dem)) sea==TRUE
        b.dem[sel]<-0
        # create block matrix for wind direction and wind coef
        b.wcoef<-matrix(NA,nrow=100+(2*b.cells),ncol=100+(2*b.cells))
        # Calculates shelter coefficient if wind direction not NA.
        # Wind direction would be NA if all values within 10km block are NA, which happens if the entire 10km block is sea
        if (sea==FALSE) b.wcoef<-windindex(b.dem,d)
        # selects data for just the 10km x 10km centre of each 30km x 30km block
        wcoef[(xmn+b.cells):(xmx-b.cells),(ymn+b.cells):(ymx-b.cells)]<-b.wcoef[(b.cells+1):(b.cells+100),(b.cells+1):(b.cells+100)]
        #print(paste("Row ",rws,", Col ",cls))
      } #rws
    }#cls
    
    # set sea cells to NA on basis of dem.m
    wcoef[which(is.na(dem.m))]<-NA
    
    # Convert and write as raster tif files
    wcoef.r<-raster(wcoef,template=dem.map)
    wcoef.r<-crop(wcoef.r,e.dem)
    par(mfrow=c(1,1))
    plot (wcoef.r,main=paste("Shelter index map: ",d,sep=""))
    out.file<-paste(dir_shelter,"Shelter_",sprintf("%03d",d,sep=""),"_deg.tif",sep="")
    print (out.file)
    writeRaster(wcoef.r,file=out.file,overwrite=T)         
  } # end direction loop
} # end function

########################################################################################

# PERCENT LAND MAPS
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

percent_land_maps<-function(dembuf,grid5km.r,dem,radius,dir_percland) {
  # If rasters not already suitably cropped
  #radius<-10000
  #xmn<-xmin(e.dem)-radius
  #xmx<-xmax(e.dem)+radius
  #ymn<-ymin(e.dem)-radius
  #ymx<-ymax(e.dem)+radius
  #e<-extent(xmn,xmx,ymn,ymx)
  #dem.buf<-crop(demuk,e)
  #grid5km.r<-crop(grid5km.r,e.dem)
  print("dembuf")
  print(dembuf)
    print("dem")
  print(dem)
	print("grid5km.r")
  print(grid5km.r)
  # Calculate % land for each 100m cell in dem.buf
  coast.r<-percent_land(dembuf,radius) 
  coast.r<-crop(coast.r,dem) # includes values for many 100m cells that are sea
  print(coast.r)
  
  #  i. Write 5km map from central 100m cells (before masking sea cells)
  vals<-rasterToPoints(grid5km.r)
  # extract %land from 100m cell matching 5km xy coords - mean calculated if several cells extracted
  coast5km<-raster::extract(coast.r,vals[,1:2]) 
  coast5km.r<-rasterize(vals[,1:2],grid5km.r,coast5km,fun=mean,background=NA)
  plot(coast5km.r,main="% land 5km grid")
  #out.file<-paste(dir_percland,"percent_land",radius/1000,"km_5km_histgrid.tif",sep="")
  #writeRaster(coast5km.r,file=out.file,overwrite=T)
  
  #  Set sea to NA and output 100m version
  coast100m.r<-mask(coast.r,dem)
  #plot(coast100m.r,main=paste("% land in a ",radius/1000,"km radius of all 100m cells", sep=""))
  out.file<-paste(dir_percland,"percent_land_",radius/1000,"km_100mcells.tif",sep="")
  writeRaster(coast100m.r,file=out.file,overwrite=T)
  
  # ii. Convert 5km grid back to 100m cell raster
  coast5km100m.r<-disaggregate(coast5km.r,50)
  coast5km100m.r<-mask(coast5km100m.r,dem) 
  # ...but is missing some 100m land cells where no overlapping 5km cell
  plot(coast5km100m.r,main="% land reference values for 100m cells - missing values")
  
  # iii. Fill 5km cells without values but containing 100m land cells
  plandref5km.r<-fill.5km.map(coast5km.r,dem,grid5km.r) 
  plot(plandref5km.r,main="% land 5km cells - cells missing historic data set to nearest values")
  out.file<-paste(dir_percland,"percent_land_",radius/1000,"km_lref_5kmgrid.tif",sep="")
  writeRaster(plandref5km.r,file=out.file,overwrite=T)
  
  # iv. Convert to 100m cells and mask with sea
  plandref100m.r<-disaggregate(plandref5km.r,50)
  print(plandref100m.r)
  plandref100m.r<-mask(plandref100m.r,dem)
  plot(plandref100m.r,main="% land ref values")
  
  out.file<-paste(dir_percland,"percent_land_",radius/1000,"km_lref_100mgrid.tif",sep="")
  writeRaster(plandref100m.r,file=out.file,overwrite=T)
  
} # end function

#####################################################################
# INVERSE LAND:SEA RATIO
# Calculates inv_land:sea ratio for different wind directions
# Adapted from Ilya's program
# Input:wind direction, sst for time t - sst must include 10kmbuffer around area of interest
# Output: rasters of ls ratio at 100m for each wind direction  

#####################################################################
# FUNCTION -  calculate end x/y coordinates for vector 
# INPUT: angle (degrees) and radius of vector
# OUTPUT: xy coordinates 
vectorxy <- function(angle,radius) #
{ angle<-angle*(pi/180)
x<-sin(angle)*radius
y<-cos(angle)*radius
xy<-cbind(x=x,y=y)
return(xy)
}

# FUNCTION - to calculate if vector v2 clockwise from vector v1 defined using x/y coordinates
# INPUT: v1 x & y coord, v2 x & y 
# OUTPUT: TRUE=clockwise
# ASSUMPTION: assumes vectors within 180 degres - false results if obtuse angle
# Ref: http://stackoverflow.com/questions/13652518/efficiently-find-points-inside-a-circle-sector
areClockwise<-function(v1.x,v1.y,v2.x,v2.y){ # is v2  clockwise of  v1 
  clockwise<- (-v1.x*v2.y )+(v1.y*v2.x) >0 
  return(clockwise)
}

# FUNCTION - to calculate if vector v2 clockwise from vector v1 defined using x/y coordinates
# INPUT: v1 x & y coord, v2 x & y 
# OUTPUT: TRUE=clockwise  Identical vectors=FALSE
# ASSUMPTION: assumes vectors within 180 degres - false results if obtuse angle
# EXPANDED version of above
areClockwise2<-function(v1.x,v1.y,v2.x,v2.y){ # is v2 clockwise of  v1 
  nv.x<--v1.y
  nv.y<-v1.x
  projx<-v2.x*nv.x  ;#print(projx); plot(projx)
  projy<-v2.y*nv.y ;#print(projy); plot(projy)
  clockwise<- (projx+projy) <0 ; #print(projx+projy);plot(projx+projy)
  return(clockwise)
}

# FUNCTION - to calculate filter for coastal effect
# INPUT: blank matrix of results, radius (distance of coastal effect in m), start/end vectors
# Calculates inverse distance measure
# ALL cells for which any of the 4 corners fall within sector are included in filter
# OUTPUT: filled matrix with weighted index that sums to 1.0 for all m
# Ref: https://scrogster.wordpress.com/2012/10/05/applying-a-circular-moving-window-filter-to-raster-data-in-r/

invdist.sector.filter<-function(m, radius, res=100, v1, v2) {  
  # calculate end xy for vectors defining sector
  v1xy<-vectorxy(v1,radius) # start of sector vector
  v2xy<-vectorxy(v2,radius) #  end of sector (clockwise of start) vector
  #print(paste("sector = ",v1,v2,sep="  "))
  #print(paste("v1xy=",v1xy,sep=""))
  #print(paste("v2xy=",v2xy,sep=""))
  
  for (row in 1:nrow(m)){
    for (col in 1:ncol(m)){
      y<-as.numeric(dimnames(m)[[1]])[row] 
      x<-as.numeric(dimnames(m)[[2]])[col] 
      ymax<-y+(res/2) ; ymin<-y-(res/2) 
      xmax<-x+(res/2) ; xmin<-x-(res/2) 
      
      dist<-sqrt(y^2 + x^2)
      #cwv2<-areClockwise(v2xy[,"x"],v2xy[,"y"],x,y)
      #ccwv1<-areClockwise(x,y,v1xy[,"x"],v1xy[,"y"])
      # Check to see if angle with either v1 or v2 is >90degrees
      if (  c(v1xy[,"x"],v1xy[,"y"]) %*% c(x,y) <=0 |  c(v2xy[,"x"],v2xy[,"y"]) %*% c(x,y) <=0)  obtuse=TRUE else obtuse<-FALSE
      # if (floor(x/1000)==ceiling(x/1000) & floor(y/1000)==ceiling(y/1000)) print(paste(x," ", y,"  ", dist,"    ", obtuse))
      # Check if grid cell falls into sector
      if ( (dist<=radius) & (obtuse==FALSE) &
           ( areClockwise2(v2xy[,"x"],v2xy[,"y"],xmax,ymax ) |
             areClockwise2(v2xy[,"x"],v2xy[,"y"],xmax,ymin ) |
             areClockwise2(v2xy[,"x"],v2xy[,"y"],xmin,ymax ) |
             areClockwise2(v2xy[,"x"],v2xy[,"y"],xmin,ymin ) )  &
           ( !areClockwise2(v1xy[,"x"],v1xy[,"y"],xmax,ymax ) |
             !areClockwise2(v1xy[,"x"],v1xy[,"y"],xmax,ymin ) |
             !areClockwise2(v1xy[,"x"],v1xy[,"y"],xmin,ymax ) |
             !areClockwise2(v1xy[,"x"],v1xy[,"y"],xmin,ymin ) )
      )
      { m[row,col]<-1/(dist/radius)  }     
    }
  }
  #sect.cells<-length(which(m>0))
  #index.m<-m/sum(m) # sum of index.m=1.0 (if all sector cells qualify ie are sea)
  index.m<-m/sum(m) ; print(sum(m))
  return(index.m)
}

#####################################################################
# ELEVATION DIF MAPS
# Calculates mean elevation for 5km grid cells then interpolates back to 100m dem cells
# Write elevation difference maps E5km-E100m
# Compare percentland_map_function

interpolated.elevation.dif.map<-function(dem,uk5kmgrid.r,dir_terrain){
  refgrid.r<-crop(uk5kmgrid.r,dem)
  
  # A. Not possible to calculate 5km map from central 100m cells as many of these are sea cells for coastal regions
  elev5km.r<-mask(aggregate(dem,fact=c(50,50),fun=mean),refgrid.r) # calculate mean elevation and remove sea cells according to refgrid.r
  plot(elev5km.r,main="Elev mean 5km")
  
  # B Interpolate using tps back to 100m cells for whole dem
  elev100m.r<-tps.resample(elev5km.r,dem,maskoutput=TRUE) 
  plot(elev100m.r,main=" 5km mean elevations interpolated to 100m")
  out.file<-paste(dir_terrain,"ref_elevation_100m.tif",sep="")
  print(out.file)
  writeRaster(elev100m.r,file=out.file,overwrite=T)
  
  # C Calculate elevation difference and write file
  ediff.r<-elev100m.r-dem
  plot(ediff.r,main="Eref-DEM")
  out.file<-paste(dir_terrain,"eref-edem_100m.tif",sep="")
  print(out.file)
  writeRaster(ediff.r,file=out.file,overwrite=T)
} # end function

#####################################################################
# LDIF MAP 
# Calculate Ldif maps for different wind directions
# Input:  Lref rasters of % land of 5km referencecells within radius at resolution of 100m
#         rasters of inverse land-sea ratio (L) for dif wind directions at 100m
# Output: raster of Ldif (%land index - inv_ls index) for each wind direction
# Prev Programs: perc_land_maps_function, inv_lsratio_maps_function - ASSUMES maps of same dimensions and resolution

ldif.maps<-function(interval,radius,dir_percland,dir_lsratio,dir_ldif){
  # Load % land maps
  #radius<-10000
  # Use %land maps calaculated for each 100m cell
  in.file<-paste(dir_percland,"percent_land_",radius/1000,"km_100mcells.tif",sep="") # change file here for different lref 
  print(in.file)
  pland100m.r<-raster(in.file)
  
  # For each wind direction and coastal effect map load inv_ls ratio maps
  for (direction in seq(0,360,interval))
  {
    # Load coast index maps and extract values for centre of 5km cells
    in.file<-paste(dir_lsratio,"Coast_Index_Sector_invdist_",sprintf("%03d",direction,sep=""),"_",(radius/res),"km_10deg.tif",sep="")
    print(in.file)
    inv.lsratio<-raster(in.file) # load raster of lsratio
    # Calculate and save Ldif
    Ldif.r<- pland100m.r-inv.lsratio
    plot(Ldif.r,main=paste("Ldif for wind direction ",direction,sep=""))
    out.file<-paste(dir_ldif,"ldif_v2_",direction,"deg_from_percland_in_",radius/1000,"km.tif",sep="")
    print(out.file)
    writeRaster(Ldif.r,file=out.file,overwrite=TRUE)
  }# end for
} # end function
#####################################################################

#####################################################################
# CODE calling above functions
#####################################################################


# Terrain maps - used by radprog_blocks
# slope<-terrain(dembuf, opt='slope', unit='degrees')
# aspect<-terrain(dembuf, opt='aspect', unit='degrees')
# plot(aspect,main="Aspect")
# plot(slope,main="Slope")
# writeRaster(slope,file=paste(dir_terrain,"slope.tif",sep=""),overwrite=TRUE)
# writeRaster(aspect,file=paste(dir_terrain,"aspect.tif",sep=""),overwrite=TRUE)

# elevation dif maps using mean of 5km reference cells (central xy can be missing) - Prog: elevdif_map_function
# For each 100m grid cell interpolate vs nearest cells and as function of elevation difference between 100m elevation and each 5km elevation
# interpolated.elevation.dif.map(dem,grid5km.r,dir_terrain)

# Create wind shelter maps - Prog: Wind_sheltermaps_functions
# Input: requires 20km buffer area. Output tif files saved to dir_shelter
# interval<-1 # = division of wind direction in degrees
# wind_sheltermaps(demuk,e.dem,interval,dir_shelter)

# Calculate % land in radius of each cell used in historic temp calculations - Prog: percent_land_maps_function (NB: will differ in UKCP09WG analysis)
# Input: requires buffer = radius. Output 5km and 100m grid tif files saved to dir_percland
# radius<-10000 #10km or20km
 percent_land_maps(dembuf,grid5km.r,dem,radius,dir_percland) 

# Calculate inv.land:sea ratio using circle sector 
 landsea.buffer<-calc(dembuf,function(x) ifelse(is.na(x),0,1))
 plot(landsea.buffer)
 res<-100
 for (angle in seq(240,360,10))  {
   #  angle<-50 # e.g direction from which wind blowing or... 
   
   # Define sector width = 1+ degrees and resulting start/end vectors v1/v2
   sect.width<-10
   sect.rad<-sect.width/2
   if (angle>5) {v2<-angle-sect.rad
   } else {v2<-angle-sect.rad}
   if (angle<356) {v1<-angle+sect.rad
   } else {v1<-(angle-(360-sect.rad))}
   # Define results matrix - ATTENTION - ORIENTATION!!!
   m<-matrix(0,nrow=1+(2*radius/res),ncol=1+(2*radius/res))
   dimnames(m)[[1]]<-rev(seq(-radius,radius,by=res))
   dimnames(m)[[2]]<-seq(-radius,radius,by=res)
   # Calculate filter/index
   #filter<-sector.filter(m,radius,v1,v2) 
   filter<-invdist.sector.filter(m,radius,res=100,v1,v2) 
   filter.r<-raster(filter)
   plot(filter.r,main=angle)
   # Apply to full raster and confirm sea cells as NA via mask
   new.r<-focal(landsea.buffer,filter)
   # new.r<-focal(land.buffer,filter,fun=function(x,y){ifelse(y!=0,grid)}
   #plot(new.r); 
   coast.r<-mask(new.r, dembuf)
   coast.r<-crop(coast.r,dembuf)
   plot(coast.r,main=paste("Coast Effect where wind dir= ",angle))
   # Write coastal index map as file
   out.file<-paste(dir_coast,"Coast_Index_Sector_invdist_",sprintf("%03d",angle,sep=""),"_",(radius/res),"km_",sect.width,"deg.tif",sep="")
   print(out.file)
   writeRaster(coast.r,file=out.file,overwrite=TRUE)
 } # end for angle loop

# %landref-inv.land.sea by wind.dir (OR Lref-Lcell) - Prog: ldif_maps_function
 ldif.maps(interval,radius,dir_percland,dir_coast,dir_ldif)


