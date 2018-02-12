##########################################################################################
# PROCESS INPUT DATA
# Calculate all DEM derived maps for use in modelling
# Outputs constant across time scale
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
inv.dist<-function(x)
{
  d<-c(1:length(x))
  d<-1/d
  id<-d*x
  ido<-sum(id)
  ido<-ido/5.187378 # Source of value???
  ido
}

inv.ls<-function(landsea,direction)
{
  store<-array(0,dim=c(100,100,101))
  store[,,1]<-landsea[101:200,101:200]
  for (i in 1:100)
  {
    xshift<-round(i*sin(direction*(pi/180)),0)
    yshift<-round(i*cos(direction*(pi/180)),0)
    yshift<-yshift*(-1)
    store[,,(i+1)]<-landsea[(101+yshift):(200+yshift),(101+xshift):(200+xshift)]
  }
  storev<-array(store,dim=c(100*100,101))
  distance<-apply(storev,1,inv.dist)
  id<-array(distance,dim=c(100,100))
  id
}

#####################################################################
# Main Function called
# Uses min 10km buffer around area of interest
#####################################################################
inv.lsmaps<-function (dembuf,dem,interval,dir_lsratio){
  
  for (direction in seq(0,360,interval))
  {
    print(direction)  # wind direction (origin)
    # create a template raster for storing inverse land sea files of same dimensions as dem
    inv.lsratio<-raster(extent(dem),res=res(dem),nrows=nrow(dem),ncols=ncol(dem))
    inv.lsratio<-setValues(inv.lsratio,rep(0,ncell(inv.lsratio)))
    
    # NB dataset too big to work out inverse land -sea ratio in entirity, so does for each 10km square in turn 
    numcols<-( ((xmax(inv.lsratio)-xmin(inv.lsratio)) /10000 )-1) 
    numrows<-( ((ymax(inv.lsratio)-ymin(inv.lsratio)) /10000 )-1)
    par(mfrow=c(3,3)) # delete
    for (b.col in 0:numcols ){
      for (b.row in 0:numrows ){
        # define 30km sq block
        xmn<-xmin(inv.lsratio)+(b.col*10000)-10000
        xmx<-xmn+30000
        ymn<-(b.row*10000)+ymin(inv.lsratio)-10000
        ymx<-ymn+30000
        e.block<-extent(c(xmn,xmx,ymn,ymx))
        # test for single block
        #xmn<-250000;xmx<-280000
        #ymn<-30000; ymx<-60000
        #e.block<-extent(c(xmn,xmx,ymn,ymx))
        # end test code
        r.block<-crop(dembuf,e.block)
        
        # convert to land = 1, sea = 0
        v<-getValues(r.block,format="matrix")
        v<-ifelse(is.na(v),0,1)
        
        id<-inv.ls(v,direction) # map of index - sea cells=0
        r2<-raster(id,xmn+10000,xmx-10000,ymn+10000,ymx-10000)
        plot(r2)
        if ( !all(v==0)){ 
          inv.lsratio<-raster::mosaic(inv.lsratio,r2,fun=max) 
        } # end if       
      } # end b.row 
      plot(inv.lsratio)
      print(paste("Col: ",b.col," Row ",b.row,sep=""))
    }# end b.col
    
    #crop to dem and convert sea to NA
    #plot(inv.lsratio)
    inv.lsratio<-mask(crop(inv.lsratio,dem),dem)
    plot(inv.lsratio)
    
    # Write lsratio raster file
    fileout<-paste(dir_lsratio,"invratio_",direction,"deg.tif",sep="")
    print(fileout)
    writeRaster(inv.lsratio,file=fileout,overwrite=T)
  } # end dd loop
  
} # end inv_lsmaps function

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
    in.file<-paste(dir_lsratio,"invratio_",direction,"deg.tif",sep="")
    print(in.file)
    inv.lsratio<-raster(in.file) # load raster of lsratio
    # Calculate and save Ldif
    Ldif.r<- pland100m.r-inv.lsratio
    plot(Ldif.r,main=paste("Ldif for wind direction ",direction,sep=""))
    out.file<-paste(dir_ldif,"ldif_",direction,"deg_from_percland_in_",radius/1000,"km.tif",sep="")
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

# inv land:sea ratio by wind dir - Prog: inv_lsratio_maps_function
# Input: requires 10km buffer around area of interest - assumes radius=10000
# inv.lsmaps(dembuf,dem,interval,dir_lsratio)

# %landref-inv.land.sea by wind.dir (OR Lref-Lcell) - Prog: ldif_maps_function
 ldif.maps(interval,radius,dir_percland,dir_lsratio,dir_ldif)


