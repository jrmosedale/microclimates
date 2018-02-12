##########################################################################################
# PROCESS INPUT DATA
# Progs:  Wind_sheltermaps_functions
#         percent_land_maps_function
#         inv_lsratio_maps_function
#         ldif_maps_function
#         elevdif_map_function

#         wind_1_readdata
#         t5km_to_matrix
#         rel.hum.v2
#         calc.cal.hrly
#         longwav_grids
#         sst_downscale_function 

##########################################################################################
# 1. Calculate and write dem derived MAP files
# Outputs constant across time scale

# Terrain maps - used by radprog_blocks
slope<-terrain(dembuf, opt='slope', unit='degrees')
aspect<-terrain(dembuf, opt='aspect', unit='degrees')
plot(aspect,main="Aspect")
plot(slope,main="Slope")
writeRaster(slope,file=paste(dir_terrain,"slope.tif",sep=""),overwrite=TRUE)
writeRaster(aspect,file=paste(dir_terrain,"aspect.tif",sep=""),overwrite=TRUE)

# Create wind shelter maps - Prog: Wind_sheltermaps_functions
# Input: requires 20km buffer area. Output tif files saved to dir_shelter
interval<-10 # = division of wind direction in degrees
wind_sheltermaps(demuk,e.dem,interval,dir_shelter)

# Calculate % land in radius of each cell used in historic temp calculations - Prog: percent_land_maps_function (NB: will differ in UKCP09WG analysis)
# Input: requires buffer = radius. Output 5km and 100m grid tif files saved to dir_percland
radius<-10000 #10km or20km
percent_land_maps(dembuf,grid5km.r,e.dem,radius,dir_percland) 

# inv land:sea ratio by wind dir - Prog: inv_lsratio_maps_function
# Input: requires 10km buffer around area of interest - assumes radius=10000
inv.lsmaps(dembuf,dem,interval,dir_lsratio)

# %landref-inv.land.sea by wind.dir (OR Lref-Lcell) - Prog: ldif_maps_function
ldif.maps(interval,radius,dir_percland,dir_lsratio,dir_ldif)

# elevation dif maps using mean of 5km reference cells (central xy can be missing) - Prog: elevdif_map_function
#For each 100m grid cell interpolate vs nearest cells and as function of elevation difference  between 100m elevation and each 5km elevation
elevation.dif.map(dem,grid5km.r,dir_terrain)

##########################################################################################
# 2. Unzip data files 
# Various programs to be run only once (eg unzip, download, basic provessing of data sources)


# SEA SURFACE PRESSURE data - daily 0.25 deg - load and crop raster bands to area of interest
gzfile<-paste(dir_pressure025,"pp_0.25deg_reg_v11.0.nc.gz",sep="")
ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")
#gunzip(filename=gzfile, destname=ncfile, overwrite=TRUE)

# WIND DATA - Write single files for wind u and v (all times) - Prog: wind_1_readdata

# Extract Radiation datafiles. Prog: extract_tar_ncdf_function
# WARNING - currently extracts ALL tar files - produces 8750 files per year per factor
extract.tar.to.ncdf(dir_dnitar,dir_dnigz,dir_dni)
extract.tar.to.ncdf(dir_sistar,dir_sisgz,dir_sis)
extract.tar.to.ncdf(dir_caltar,dir_calncgz,dir_calnc)


##########################################################################################

# 3. Run processes specific for chosen time period 
# Progs:t5km_to_hourly_blocks
#       sst_downscale_function

# Set time range for which data will be analysed - uses JDdmy function
start.year<-1992
start.month<-7
start.day<-11
hr<-0
end.year<-1992
end.month<-7
end.day<-20
start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)

### Downscale SST to daily 5km grid - Prog: sst_downscale_function ###

# IMPORTANT: NEEDS TO BE FOR 20KM BUFFER REGION TO ALLOW CALC OF UPWIND SST
in.file<-paste(dir_sst,"HadISST_sst.nc",sep="") # in.file required for sst.spdownsc
ncfile<-nc_open(paste(dir_sst,"HadISST_sst.nc",sep="")) # summary of file variables,  dimensions, attributes

# 1. WRITE monthly files for buffered uk region - normally 10km buffer - NEED DATA MONTH BEFORE START!!!
# Calculate jd for start of PREVIOUS month to start.jd *
# ?? Uses land mask of grid5km.r
if (start.month==1) { 
  year1<-start.year-1
  month1<-12 } else {
    year1<-start.year
    month1<-start.month-1 }
if (end.month==12) { 
  year2<-end.year+1
  month2<-1 } else {
    year2<-end.year
    month2<-end.month+1 }

sst.spdownsc(year1,month1,year2,month2,grid5kmbuf.r,dir_sstm)  

# 2. WRITE daily files 
for (jd in start.jd:end.jd) {
  sst.time.int(jd, dir_sstm, dir_ssth)
} # end day

# 3. delete monthly files
for(y in year1:year2){
  for(m in month1:month2){
    filename<-paste(dir_sstm,"sst_",y,"_",m,".tif",sep="")
    print(filename)
    file.remove(filename)
  }}
nc_close(ncfile)
###   END  ###

## Downscale daily TEMPERATURE data to hourly  Prog: t5km_to_hrmatrix ###
# WRITES t5km.day dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),".r"
# REQUIRES: elevdif_map FUNCTIONS to fill in values for missing 5km cells with LAND cells at 100m
### Unzip TEMPERATURE files required - start.jd -1?? 
unzip_yearfiles(start.jd,end.jd,dir_zip,dir_temp)
# includes data for previous day to start.day
hourly_temperatures(start.jd-1,end.jd+1,dir_temp,dir_hrtemp,grid5km.r) # +1 end day so can be used for RH

#Downscale REL HUMIDITY to hourly 5km Prog: rel.hum.v2
# Writes rh.day as dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R"
rh.hourly(start.jd-1,end.jd,dir_rh,dir_rh5km,dir_hrtemp,grid5km.r)

# Downscale and impute CLOUD ALBEDO - WARNING - LONG TIME - RUN SEPERATELY??  Prog: calc.cal.hrly
# WRITES calimp.day as dir_calimp,"CALimp_5km_",year,"_",month,"_",day,".R
cal.5km.impute(dir_cal,dir_calimp,start.jd,end.jd,grid5km)

# Calculate LONG WAVE RADITAION at 5km hourly res ffrom CAL, RH T  Prog: longwav_grids
write_lwr_dayfiles(start.jd,end.jd,grid5km.r)

# Pressure - assign file - already unpacked
p.ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")





##########################################################################################
# MAP FUNCTIONS USED BY ABOVE

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
        #plot(r2)
        if ( !all(v==0)){ 
          inv.lsratio<-raster::mosaic(inv.lsratio,r2,fun=max) 
        } # end if       
      } # end b.row 
      #plot(inv.lsratio)
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
  in.file<-paste(dir_percland,"percent_land_",radius/1000,"km_lref_100mgrid.tif",sep="")
  print(in.file)
  pland100m.r<-raster(in.file)
  
  # For each wind direction and coastal effect map load inv_ls ratio maps
  for (direction in seq(0,350,interval))
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
# DATA PROCESSING FUNCTIONS
#####################################################################

# GENERATE HOURLY TEMPERATURE DATA AT %KM RESOLUTION

# Calculates sunrise and sunset (set Timezone and DST to a vector of 0)
# uses julian day input rather than doy and year
# CHECK - how jd varies from start to end of day
# Input / output as vectors of same length
suntimes.grid<-function(JD,Lat,Long,Timezone,DST){
  J<-JD
  lw<-Long*-1
  n<-J-2451545-0.0009-(lw/360)
  n<-floor(n)+0.5
  sn<-2451545+0.0009+(lw/360)+n
  msa<-(357.5291+0.98560028*(sn-2451545))%%360
  eoc<-1.9148*sin(msa*pi/180)+0.02*sin(2*msa*pi/180)+0.0003*sin(3*msa*pi/180)
  ecl<-(msa+102.9372+eoc+180); ecl<-ecl%%360
  st<-sn+(0.0053*sin(msa*pi/180))-(0.0069*sin(2*ecl*pi/180))
  d<-asin(sin(ecl*pi/180)*sin(23.45*pi/180))
  cos.has<-((sin(-0.83*pi/180)-sin(Lat*pi/180)*sin(d))/(cos(Lat*pi/180)*cos(d)))
  h.set<-vector(length=length(JD)); h.rise<-h.set; dl<-h.set
  # next three lines may provoke warnings in some cases?
  has<-acos(cos.has)
  J.set<-2451545+0.0009+(((has*180/pi+lw)/360)+n+0.0053*sin(msa*pi/180))-0.0069*sin(2*ecl*pi/180)
  J.rise<-st-(J.set-st)
  
  ifelse(cos.has^2<1,
         h.set<-(J.set%%1)*24+Timezone+DST,
         ifelse(cos.has>1,
                h.set<-12,
                ifelse(cos.has<(-1),
                       h.set<-0,
                       warning("Cos.has case not found") ) ) )
  
  ifelse(cos.has^2<1,
         h.rise<-(J.rise%%1)*24+Timezone+DST,
         ifelse(cos.has>1,
                h.rise<-12,
                ifelse(cos.has<(-1),
                       h.rise<-0,
                       warning("Cos.has case not found") ) ))
  
  ifelse(cos.has^2<1, dl<-(J.set-J.rise)*24,
         ifelse(cos.has>1,
                dl<-0,
                ifelse(cos.has<(-1),
                       dl<-24,
                       print("Cos.has case not found") ) ) )    
  
  if(any(dl==0)) {warning("sun below horizon for 24 hours")}   
  if(any(dl==24)) {warning("sun above horizon for 24 hours")}    
  
  sun.vars<-data.frame(sunrise=h.rise,sunset=h.set,daylight=dl)
  sun.vars
}
# calculates sunrise
sunrise.grid<-function(JD,Lat,Long,Timezone,DST){
  sun.rise<-suntimes.grid(JD-1,Lat,Long,Timezone,DST)[,1]
  sun.rise
}
# calculates sunset
sunset.grid<-function(JD,Lat,Long,Timezone,DST){
  sun.set<-suntimes.grid(JD,Lat,Long,Timezone,DST)[,2]
  sun.set
}
# generates hourly values for use between dawn and hotest part of day (~13:35)
generate.hrtemps.grid.day1<-function(min.temp,max.temp,day.length)
{
  x<-matrix(rep(0:23,NROW(min.temp)),nrow=NROW(min.temp), ncol=24, byrow=TRUE)
  A=(max.temp-min.temp)/2
  fr=1/(day.length*1.5)
  phase=13.58989
  of=A+min.temp
  y<-A*cos(2*pi*fr*(x-phase))+of
  y
}
# generates hourly values for use between hotest part of day (~13:35) and midnight
generate.hrtemps.grid.day2<-function(min.temp,max.temp,sun.rise)
{
  x<-matrix(rep(0:23,NROW(min.temp)),nrow=NROW(min.temp), ncol=24, byrow=TRUE)
  A=(max.temp-min.temp)/2
  lengt<-(24-13.58989)+sun.rise
  fr=1/(lengt*1.5)
  phase=13.58989
  of=A+min.temp
  y<-A*cos(2*pi*fr*(x-phase))+of
  y
}
# generates hourly values between midnight and dawn. Coolest part of day ~12 mins before dawn
generate.hrtemps.grid.night<-function(min.temp,max.temp,sun.rise,sun.set)
{
  x<-matrix(rep(0:23,NROW(min.temp)),nrow=NROW(min.temp), ncol=24, byrow=TRUE)
  A=(max.temp-min.temp)/2
  day.length<-sun.set-sun.rise
  lengt<-(24-13.58989)+sun.rise
  fr=1/(lengt*1.5)
  phase=sun.rise-0.2044346
  of=A+min.temp
  y<-A*sin(2*pi*fr*(x-phase)-(2*pi/4))+of
  y
}
# generates hourly values
#inputs:
# min.temp = minimum daily temperature
# max.temp = maximum daily temperature
# next.min = minimum daily temperature the next day
# prev.max = maximum daily temperature the previous day
# sun.rise = sunrise time expressed a decimal hour (24hrs) (see functions above)
# sun.set  = sunset time expressed a decimal hour (24hrs)  (see functions above)
# output:
# matrix of 24 values(cols) corresponding to estimated temperature in each hour for each grid cell (row)  
#(first value=midnight, last value= 23:00 hrs)
generate.hrtemps.grid<-function(min.temp,max.temp,next.min,prev.max,sun.rise,sun.set)
{
  x<-matrix(rep(0:23,NROW(min.temp)),nrow=NROW(min.temp), ncol=24, byrow=TRUE)
  day.length<-sun.set-sun.rise
  day1<-generate.hrtemps.grid.day1(min.temp,max.temp,day.length)
  day2<-generate.hrtemps.grid.day2(next.min,max.temp,sun.rise)
  night<-generate.hrtemps.grid.night(min.temp,prev.max,sun.rise,sun.set)
  x<-day1
  x[,15:24]<-day2[,15:24]
  # Replace x when column<= sun.rise and replace with night
  x<-ifelse(col(x)<=sun.rise,night,x) 
  return(x)
}

# Unzip yearly temperature files in dir_zip required for time period and save daily files in dir_temp
unzip_yearfiles<-function(start.jd,end.jd,dir_zip,dir_temp)
{
  # Unzip year files of temperature data here
  zip.file<-paste(dir_zip,"MinTemp_", DMYjd(start.jd)$year,".zip", sep="")
  print (paste("Unzipping ",zip.file,sep=""))
  unzip(zip.file, exdir=dir_temp)
  zip.file<-paste(dir_zip,"MaxTemp_", DMYjd(start.jd)$year,".zip", sep="")
  print (paste("Unzipping ",zip.file,sep=""))
  unzip(zip.file, exdir=dir_temp)
  
  # Check if start or end dates are 1st or last of a year and unzips data for these years as well
  if ( DMYjd(start.jd)$day==1 & DMYjd(start.jd)$month==1) { # then extract previous year of data
    zip.file<-paste(dir_zip,"MinTemp_", DMYjd(start.jd)$year-1,".zip", sep="")
    print (paste("Unzipping ",zip.file,sep=""))
    unzip(zip.file, exdir=dir_temp) 
    zip.file<-paste(dir_zip,"MaxTemp_", DMYjd(start.jd)$year-1,".zip", sep="")
    print (paste("Unzipping ",zip.file,sep=""))
    unzip(zip.file, exdir=dir_temp)
  }
  if ( DMYjd(end.jd)$day==31 & DMYjd(end.jd)$month==12) { # then extract next year of data
    zip.file<-paste(dir_zip,"MinTemp_", DMYjd(end.jd)$year+1,".zip", sep="")
    print (paste("Unzipping ",zip.file,sep=""))
    unzip(zip.file, exdir=dir_temp) 
    zip.file<-paste(dir_zip,"MaxTemp_", DMYjd(end.jd)$year+1,".zip", sep="")
    print (paste("Unzipping ",zip.file,sep=""))
    unzip(zip.file, exdir=dir_temp)
  }
} # end function

#######################################################################################
# Main Function
# Writes daily 5km temperature files containing hourly temperature data for whole 5km area
# USES: fill.5km.map and relted FUNCTIONS from elevdif_map_functions
# i.e. records from nearest 5km historic data grid cell
#######################################################################################
hourly_temperatures<-function(start.jd,end.jd,dir_temp,dir_hrtemp,grid5km.r){ 
  plothrs=FALSE
  e.dem<-extent(grid5km.r)  
  # Import temperature file available for day before start date (start.jd-1) 
  infile<-paste(dir_temp,"MaxTemp_", DMYjd(start.jd-1)$year, "-",sprintf("%02d",DMYjd(start.jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(start.jd-1)$day,sep=""),"_ACTUAL.txt", sep="")
  tmpdata.r<-raster(infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
  tmpdata.r<-crop(x=tmpdata.r,y=e.dem) # crop to geographical extent of DEM raster
  
  #Create lat/lon grid for daylength etc calculations from cropped temp data
  osgrid<-SpatialPoints(coordinates(tmpdata.r), proj4string=CRS("+init=epsg:27700"), bbox = NULL)
  os.m<-coordinates(osgrid)
  llgrid<-spTransform(osgrid,CRS("+init=epsg:4326"))
  ll.m<-coordinates(llgrid) # creates vector of coordinates from top left  grid cell to bottom right by rows
  
  for (jd in start.jd:end.jd) {  # Loop to calculate hourly files for each day
    # Define output matrix
    t5km.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
    # Read day temperature data
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),"_ACTUAL.txt", sep="")
    day.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    day.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Read previous day temperature data
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_ACTUAL.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_ACTUAL.txt", sep="")
    prev.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    prev.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Read NEXT day file 
    max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="")
    min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="")
    next.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    next.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
    
    # Crop all daily rasters  using dem as model - USE BRICK?? - COMBINE WITH BELOW?
    day.tmax<-crop(x=day.tmax,y=e.dem)
    day.tmin<-crop(x=day.tmin,y=e.dem)
    next.tmax<-crop(x=next.tmax,y=e.dem)
    next.tmin<-crop(x=next.tmin,y=e.dem)
    prev.tmax<-crop(x=prev.tmax,y=e.dem)
    prev.tmin<-crop(x=prev.tmin,y=e.dem)
    
    # OPTION - Fill missing 5km cells without temperature values but containing land cells from nearest neighbour
    day.tmax<-fill.5km.map(day.tmax,dem,grid5km.r) 
    day.tmin<-fill.5km.map(day.tmin,dem,grid5km.r) 
    next.tmax<-fill.5km.map(next.tmax,dem,grid5km.r) 
    next.tmin<-fill.5km.map(next.tmin,dem,grid5km.r) 
    prev.tmax<-fill.5km.map(prev.tmax,dem,grid5km.r) 
    prev.tmin<-fill.5km.map(prev.tmin,dem,grid5km.r) 
    
    # 1. Calculate sunrise, sunset and day length for each grid cell
    
    # Define vectors to calculate sunrise/set and daylength from lat/lon 
    lat<-ll.m[,"y"]
    long<-ll.m[,"x"]
    numcells<-length(long) # number of grid cells in map = length of vector
    
    sunup<-rep(0,numcells)
    sundown<-rep(0,numcells)
    daylength<-rep(0,numcells)
    jd.v<-rep(jd,numcells)
    #year.v<-rep(DMYjd(jd)$year,numcells)
    
    # CHECK JD calculations !!
    sunup<-sunrise.grid(jd.v,lat,long,rep(0,numcells),rep(0,numcells))
    sundown<-sunset.grid(jd.v,lat,long,rep(0,numcells),rep(0,numcells))
    daylength<-sundown-sunup
    
    ##### 2. Generate hourly temperatures using matrices
    # Create vectors of grid cell daily min/max temperature values
    tmax<-getValues(crop(x=day.tmax,y=dem)) # from top left to bottom right by row
    tmin<-getValues(crop(x=day.tmin,y=dem))
    p.tmax<-getValues(crop(x=prev.tmax,y=dem))
    n.tmin<-getValues(crop(x=next.tmin,y=dem))
    
    # Create hourly temperatures for each grid cell
    temp.hr<-generate.hrtemps.grid(tmin,tmax,n.tmin,p.tmax,sunup,sundown) 
    print(paste("Date: ", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),sep=""))
    
    # Write one hourly temperatures as single R matrix cellxhour for each day at 5km resolution 
    # file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),".r", sep="") # define file name from year,month,day,hr
    # save(temp.5km.hr, file=file.out)
    
    # Plot hourly values
    
    for (hr in 0:23) {
      t5km.day[,,hr+1]<-matrix(temp.hr[,hr+1],nrow=NROW(day.tmax),ncol=NCOL(day.tmax),byrow=TRUE) # convert vector to matrix
      # plot raster option
      if (plothrs==TRUE){
        hr.r<-raster(nrow=nrow(day.tmax),ncol=ncol(day.tmax))
        hr.r<-setValues(hr.r,temp.hr[,hr+1]) # convert matrix to raster
        temp.plot(hr.r,hr,DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year) }
    }
    
    # Write day of hrly data as 3D matrix
    file.out<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year[1], "-",sprintf("%02d",DMYjd(jd)$month[1],sep=""),"-", sprintf("%02d",DMYjd(jd)$day[1],sep=""),".r", sep="") # define file name from year,month,day,hr
    save(t5km.day, file=file.out)
    
  } # end for day loop
  
}# end function

#####################################################################
# WRITE DAILY FILES OF RELATIVE HUMIDITY - rel.hum.v2 
#####################################################################
#dir_rh5km<-"~/Documents/Exeter/Data2015/RelHumidity/rh5km/"
add.zero<-function(x)
{
  y<-x
  if (y<9) y<-paste("0",x,sep="")
  y
}

rel.to.abs<-function(rh,t)
{
  e0.1<-0.6108*exp(17.27*t/(t+237.3))
  e<-e0.1*(rh/100)
  abso<-(2165*e)/(t+273.16) # grams per metre cubed
  abso
}

abs.to.rel<-function(ab,t)
{
  e<-ab*(t+273.16)/2165
  e0.1<-0.6108*exp(17.27*t/(t+237.3))
  rh<-100*(e/e0.1)
  rh
}

# Centres raster on GML 
reslice.raster<-function(r)
{
  arr<-array(0,dim=c(73,144))
  e1<-extent(1.25,181.25,-91.25,91.25)
  e2<-extent(181.25,358.75,-91.25,91.25)
  e3<-extent(-1.25,1.25,-91.25,91.25)
  arr[,1:71]<-getValues(crop(r,e2),format="matrix")
  arr[,72]<-getValues(crop(r,e3),format="matrix")
  arr[,73:144]<-getValues(crop(r,e1),format="matrix")
  r2<-raster(arr,xmn=-178.25,xmx=181.25,ymn=-91.25,ymx=91.25)
  res(r2)<-2.5
  r2
}

reproject.r<-function(a,e,grid5km.r)
{
  r<-raster(a,xmn=xmin(e),xmx=xmax(e),ymn=ymin(e),ymx=ymax(e))
  projection(r)<-"+init=epsg:4326"
  r2<-projectRaster(r,crs="+init=epsg:27700")
  r3<-resample(r2,grid5km.r)
  #r4<-crop(r3,e2)
  r3
}

rh.plot<-function(rh,hr,day,month,year) 
{
  t<-paste("RH on ",day,"-",month,"-",year," ",hr,":00",sep="")
  brk<-c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,160,170,180,190,200)
  col<-rev(heat.colors(21))
  plot(rh,main=t,col=col,breaks=brk)
}
###########################################################################
# Write daily files of hourly 5km Rel Humidity 
# RH data in = 4xdaily global data

rh.hourly<-function(start.jd,end.jd,dir_rh,dir_rh5km,dir_hrtemp, grid5km.r,hourplot=FALSE)
{
  e<-extent(-5.75,-0.75,48.75,53.75)  # lon/lat extent of cropped rasters
  #e<-extent(353.75,359.25,48.75,53.75)
  
  for (jd in start.jd:end.jd)
  {
    # Define output matrix to hold one day of hourly relative humidity values at 5km resolution for area covered by grid5km.r
    rh.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
    
    # Reads in Relative Humidity data from yearly file if new year
    if (jd==start.jd | (DMYjd(jd)$day==1 & DMYjd(jd)$month==1)){
      infile.rh<-paste(dir_rh,"rhum.sig995.",DMYjd(jd)$year,".nc",sep="") # reads yearly data file
      #print(in.file)
      rh.brick<-brick(infile.rh) # all times
      # get times
      netRH<-nc_open(infile.rh)
      tm<-ncvar_get(netRH,"time")
    }
    
    # Read in daily Temperature file and resample to low res lat/lon
    filein1<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd)$year, "-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-", sprintf("%02d",DMYjd(jd)$day,sep=""),".r", sep="") # define file name from year,month,day,hr
    print(filein1)
    load(filein1)
    hr.temp1<-t5km.day
    
    filein2<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),".r", sep="") # define file name from year,month,day,hr
    print(filein2)
    load(filein2)
    hr.temp2<-t5km.day
    
    for (period in 0:3)
    {
      hr<-(period*6)
      # Identify rh layers and crop
      Jul.base<-JDdmy(1,1,1800)
      #Jul.actual<-JD(ISOdate(year,month,day))
      hr.val<-(jd-Jul.base)*24+hr
      sel<-which(tm==hr.val)
      rh.all0<-subset(rh.brick,sel)
      rh.all6<-subset(rh.brick,(sel+1))
      rh.all0<-reslice.raster(rh.all0)
      rh.all6<-reslice.raster(rh.all6)
      projection(rh.all0)<-"+init=epsg:4326"
      projection(rh.all6)<-"+init=epsg:4326"
      rh0<-crop(rh.all0,e)
      rh6<-crop(rh.all6,e)     
      
      # read in temperature rasters
      t0.5km<-raster(hr.temp1[,,hr+1],template=grid5km.r)
      if (period<3) t6.5km<-raster(hr.temp1[,,hr+7],template=grid5km.r)
      if (period==4) t6.5km<-raster(hr.temp2[,,1],template=grid5km.r) # 0hr00 for following day
      
      # reproject and calculate aggregate temperatures to match humidity data resolution
      projection(t0.5km)<-"+init=epsg:27700"
      projection(t6.5km)<-"+init=epsg:27700"
      t0.ll<-projectRaster(t0.5km,crs="+init=epsg:4326")
      t6.ll<-projectRaster(t6.5km,crs="+init=epsg:4326")
      t0.ll<-raster::extend(t0.ll,e) # required to allow aggregate values of rh cells only partly covered by temp data
      t6.ll<-raster::extend(t6.ll,e)
      tem<-raster(array(0,dim=c(2,2)),xmn=xmin(e),xmx=xmax(e),ymn=ymin(e),ymx=ymax(e))
      t0<-raster::resample(t0.ll,rh6)
      t6<-raster::resample(t6.ll,rh6)
      
      #########################################
      
      # Convert to absolute humidity
      abs.hum0<-rel.to.abs(getValues(rh0,format="matrix"),getValues(t0,format="matrix"))
      abs.hum6<-rel.to.abs(getValues(rh6,format="matrix"),getValues(t6,format="matrix"))
      
      # interpolate for missing hours
      abs.hum1<-(abs.hum0*5+abs.hum6*1)/6
      abs.hum2<-(abs.hum0*4+abs.hum6*2)/6
      abs.hum3<-(abs.hum0*3+abs.hum6*3)/6
      abs.hum4<-(abs.hum0*2+abs.hum6*4)/6
      abs.hum5<-(abs.hum0*1+abs.hum6*5)/6
      
      # Convert absolute humidity to 5km grid cell resolution 
      e2<-raster::extent(grid5km.r)
      abs.h1<-reproject.r(abs.hum1,e,grid5km.r)
      abs.h2<-reproject.r(abs.hum2,e,grid5km.r)
      abs.h3<-reproject.r(abs.hum3,e,grid5km.r)
      abs.h4<-reproject.r(abs.hum4,e,grid5km.r)
      abs.h5<-reproject.r(abs.hum5,e,grid5km.r)
      abs.h6<-reproject.r(abs.hum6,e,grid5km.r)
      
      # Load remaining hrs of 5km temperature data 
      t1.5km<-raster(hr.temp1[,,hr+2],template=grid5km.r)
      t2.5km<-raster(hr.temp1[,,hr+3],template=grid5km.r)
      t3.5km<-raster(hr.temp1[,,hr+4],template=grid5km.r)
      t4.5km<-raster(hr.temp1[,,hr+5],template=grid5km.r)
      t5.5km<-raster(hr.temp1[,,hr+6],template=grid5km.r)
      
      # Convert to 5km relative humidity
      hr<-period*6
      rh.day[,,hr+1]<-abs.to.rel(getValues(abs.h1,format="matrix"),getValues(t1.5km,format="matrix"))
      rh.day[,,hr+2]<-abs.to.rel(getValues(abs.h2,format="matrix"),getValues(t2.5km,format="matrix"))
      rh.day[,,hr+3]<-abs.to.rel(getValues(abs.h3,format="matrix"),getValues(t3.5km,format="matrix"))
      rh.day[,,hr+4]<-abs.to.rel(getValues(abs.h4,format="matrix"),getValues(t4.5km,format="matrix"))
      rh.day[,,hr+5]<-abs.to.rel(getValues(abs.h5,format="matrix"),getValues(t5.5km,format="matrix"))
      rh.day[,,hr+6]<-abs.to.rel(getValues(abs.h6,format="matrix"),getValues(t6.5km,format="matrix"))
      
      
    } # end period
    if (hourplot==TRUE){
      for (t in 1:24) {
        rh.plot(raster(rh.day[,,t],template=grid5km.r),t-1,DMYjd(jd)$day,DMYjd(jd)$month,DMYjd(jd)$year) 
      } }
    
    # WRITE DAILY files of 5km RELATIVE humidity
    fileout<-paste(dir_rh5km,"RH_5km_",DMYjd(jd)$year,"_",DMYjd(jd)$month,"_",DMYjd(jd)$day,".R",sep="")
    print(fileout)
    save(rh.day,file=fileout)
  } # end jd
  
} # end function

#####################################################################
# Imputes and write daily files of hourly data of CAL
#####################################################################

# FUNCTION - imputation
# input: dfi = dataset with missing values (effective Cloud albedo with missing night time values)
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

# FUNCTION - cal.5km.impute - writes 5km matrix with imputed nightime values
# Writes files of imputed CAL values from files missing nightime values
# Ensure CAL input files already cropped to correct extent

cal.5km.impute<-function(dir_cal,dir_calimp,start.jd,end.jd,grid5km,plotcal=FALSE){
  files<-list.files(dir_cal) 
  
  # Define initial cropping in lat/lon and resulting number of cols/rows to be stored
  e<-extent(c(-7,-2,49,52)) # extent for initial cropping 
  e.5km<-extent(grid5km.r)
  numrows<-nrow(grid5km.r)
  numcols<-ncol(grid5km.r)
  max.i<-(1+end.jd-start.jd)*24 # = total number of hrs to be stored
  CAL.vals<-array(NA,c(numrows,numcols,max.i))
  
  # load and store CAL files 
  i<-1
  for (t in start.jd:end.jd){
    for (hr in 0:23){ 
      year<-DMYjd(t)$year; month<-DMYjd(t)$month ; day<-DMYjd(t)$day 
      # check if file exists and print warning and record missing if not
      infile.cal<-paste(dir_cal,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),"00002UD1000101UD.nc",sep="")
      print(infile.cal)
      #if (infile.cal not in files) print("File not found") 
      
      # Resample CAL data to 5km grid
      cal.r<-raster(infile.cal)
      projection(cal.r)<-"+init=epsg:4326"
      calosgb.r<-projectRaster(cal.r,crs="+init=epsg:27700")
      #plot(calosgb.r)
      cal5km.r<-raster::resample(crop(calosgb.r,e.5km),grid5km.r)
      #plot(cal5km.r)
      #print(dim(m))
      #print (i)
      m<-getValues(cal5km.r,format="matrix")
      CAL.vals[,,i]<-m
      i<-i+1
      # else{warning(paste("Files not found: ",infile.cal,sep=""))}
    }# end hr loop
  }# end day loop
  
  # Impute missing values for each grid cell in turn
  CAL.imp<-array(NA,c(numrows,numcols,max.i))
  for (r in 1:numrows){
    for (c in 1:numcols){
      v<-CAL.vals[r,c,]
      CAL.imp[r,c,]<-imputation(v)
      print(paste("row: ",r," col: ",c,sep=""))
    } # rows
  } # cols
    
  # Write DAILY files of imputed CAL values at 5km resolution 
  for (t in start.jd:end.jd){
    CALimp.day<-array(NA,c(numrows,numcols,24))
    first<-(t-start.jd)*24+1
    last<-first + 23
    print(paste("Day: ",t," first=",first," last=",last,sep=""))
    CALimp.day[,,1:24]<-CAL.imp[,,first:last]
    
    # plot each hour of each day ?
    if (plotcal==TRUE){
      par(mfrow=c(4,3))
      for (i in 1:24){
        calsw.r<-raster(CALimp.day[,,i],template=grid5km.r)
        titletext<-paste(DMYjd(t)$day,"/",DMYjd(t)$month,"/",DMYjd(t)$year," ",i,"hr00",sep="")
        plot(calsw.r,main=titletext)
      } }
    
    fileout<-paste(dir_calimp,"CALimp_5km_",DMYjd(t)$year,"_",DMYjd(t)$month,"_",DMYjd(t)$day,".R",sep="")
    print (fileout)
    save(CALimp.day,file=fileout) 
  }# end for write files
  
} # end function

#####################################################################
# CALCULATE LW RADIATION from CAL , Temp, RH
# Units of output = MJ/m2
# Historic data: run rel.hum.v2, cal.cal.hrly & t5km_to_hrmatrix
# Requires jd functions

# FUNCTION - calculate long wave radiation from Temp, RH & CAL
lwr<-function(Temp,RH,CAL)
{
  e0<-0.6108*exp(17.27*Temp/(Temp+237.3)) # saturated vapour pressure
  ea<-e0*(RH/100) # actual vapour pressure
  moisture.absorb<-0.34-0.14*sqrt(ea)
  cloud.absorb<-1.35*(1-CAL)-0.35
  rnl<-2.043*10^-10*(Temp+273.16)^4*moisture.absorb*cloud.absorb
  rnl
}

write_lwr_dayfiles<-function(start.jd,end.jd,grid5km.r,plotlwr=FALSE)
{
  # Define output file - one day of hourly 5km data
  lwr.day<-array(0,dim=c(nrow(grid5km.r),ncol(grid5km.r),24))
  for (t in start.jd:end.jd)
  {
    year<-DMYjd(t)$year
    month<-DMYjd(t)$month
    day<-DMYjd(t)$day
    
    # Reads in daily files - all assumed to be 5km OSGB gridcells with identical extents of grid5km.r
    rh.filein<-paste(dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R",sep="")
    print(rh.filein)
    load(rh.filein) # rh.day
    
    t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),".r", sep="")
    print(t.filein)
    load(t.filein) # t5km.day
    
    cal.filein<-paste(dir_calimp,"CALimp_5km_",year,"_",month,"_",day,".R",sep="")
    print(cal.filein)
    load(cal.filein) # CALimp.day
    
    for (hr in 1:24)
    {
      # rel hum
      m.rh<-rh.day[,,hr]
      r.rh<-raster(m.rh,template=grid5km.r)
      # temperature
      m.temp<-t5km.day[,,hr]
      r.temp<-raster(m.temp,template=grid5km.r)
      # CAL 
      m.cal<-CALimp.day[,,hr]
      r.cal<-raster(m.cal,template=grid5km.r)
      
      # Calculate longwave radiation
      lwr.day[,,hr]<-lwr(getValues(r.temp,format="matrix"),
                         getValues(r.rh,format="matrix"),
                         getValues(r.cal,format="matrix"))
      
      if (plotlwr==TRUE) {
        par(mfrow=c(2,2))
        dayhr<-paste(day,"/",month,"/",year," ",hr-1,"h00",sep="")
        plot(r.cal,main=paste("Effective cloud albedo ",dayhr,sep=""))
        plot(r.rh,main=paste("Relative humdidity ",dayhr,sep=""))
        plot(r.temp,main=paste("Temperature ",dayhr,sep=""))
        plot(raster(lwr.day[,,hr],template=r.rh),main=paste("Long-wave radiation ",dayhr,sep=""))
      }
    } # for hr
    # WRITE DAILY lwr files
    file.out<-paste(dir_lwr,"lwr_",year,"_",month,"_",day,".R",sep="")
    print(file.out)
    save(lwr.day,file=file.out)
    
  } # for day
} # end function

#####################################################################
# DOWNSCALE  SST
#####################################################################
# Input: single ncdf file of Met Office Hadley Centre HadISST1 data 
# Output: Downscaled SST to hourly 100m data
# Issues/Problems: interpolation/resampling - in effect across N/S coasts
# set buffer eg to 10km
# Input data
#dir_sst<-"~/Documents/Exeter/Data2015/sst/"
#dir_sst<-"C:/Data2015/SST/"  
#gunzip(paste(dir_sst,"HadISST_sst.nc.gz",sep=""))
#dir_sstm<-paste(dir_sst,"monthly/",sep="")
#dir_ssth<-paste(dir_sst,"hourly/",sep="")
####################################################
# Functions to return year/month from sst t value
ttoyr<-function(t){
  yr<-ceiling(t/12)-1+1960
  return(yr)
}

ttomonth<-function(t){
  month<-t-((ttoyr(t)-1960)*12)
  return(month)
}

# Function to spatially downscale sst data for month t to dsgrid  
sst.spdownsc<-function(start.year,start.month,end.year,end.month,dsgrid.r,dir_sstm) {
  start.n<-(start.year-1960)*12+start.month
  end.n<-((end.year-1960)*12)+end.month
  for (n in start.n:end.n){
    # Read ncdf direct to raster
    sst.world.r<-raster(in.file,band=n) # read all of level t of ncdf file
    #plot(sst.world.r)
    # Extract relevant 1 deg cells 353-360E, 48:53N, 1960-2015 
    long.mn<--7;long.mx<--1 # -7 to -1
    lat.mn<-49; lat.mx<-52 # 49 to 52
    e.sst<-extent(c(long.mn,long.mx,lat.mn,lat.mx))
    sst.r<-crop(sst.world.r,e.sst)
    #plot(sst.r,main=paste("SST long: ",long.mn," to ",long.mx,", lat: ",lat.mn," to ",lat.mx,sep=""))
    
    # Set land cells to have same temperature as nearby sea
    #sst.r<-focal(sst.r,w=matrix(1,3,3),fun=mean,na.rm=TRUE,NAonly=TRUE)
    sst.m<-getValues(sst.r,format="matrix")
    sst.m[1,5]<-(sst.m[1,4]) # Bristol channel - count as sea = to adjacent cell (NOT mean including CHannel!)
    #sst.m[1,6]<-(sst.m[2,6]+sst.m[1,5])/2 # keep as NA all land cells
    sst.r<-raster(sst.m,template=sst.r)
    #plot(sst.r)
    
    # Convert projection to OSGB
    projection(sst.r)<-"+init=epsg:4326"
    sst.os<-projectRaster(sst.r,crs="+init=epsg:27700")
    
    # resample to 5km resolution - cf with ngb method
    sst5km.r<-tps.resample(sst.os,dsgrid.r,msk=FALSE)
    #sst5km.r<-raster::resample(sst.os,dsgrid.r,method="bilinear",na.rm=TRUE)
    # sst100m.r<-resample(sst.os,dem.buffer,method="bilinear",na.rm=TRUE) # alternative to 100m
    
    #sst5km.r<-mask(sst5km.r,dsgrid.r,inverse=TRUE) # No masking to allow future downscaling to higher coastal resolutions
    # sst100m.r<-mask(sst5km.r,dsgrid.r,inverse=TRUE)  ; plot(sst100m.r)
    
    # write raster file at 5km resolution for every month
    tl<-paste("Year: ",ttoyr(n)," Month: ",ttomonth(n),sep="")
    plot(sst5km.r,main=tl)
    fileout<-paste(dir_sstm,"sst_",ttoyr(n),"_",ttomonth(n),".tif",sep="")
    print(fileout)
    writeRaster(sst5km.r,file=fileout,overwrite=T)
  }      
}# end of function

##########################################################################################
# 2. TEMPORAL DOWNSCALING FUNCTIONS to DAILY  time period (could also be used for hourly period if necessary)

nday.in.month <- function(date)
{
  m<-format(date, format="%m")
  while(format(date, format="%m") == m) date <- date + 1
  return(as.integer(format(date - 1, format="%d")))
}

nday.in.month.jd <- function(jd)
{
  m<-DMYjd(jd)$month
  jd1<-jd-DMYjd(jd)$day+1 # jd for 1st of the month
  while(DMYjd(jd1)$month == m) jd1 <- jd1 + 1
  return(DMYjd(jd1-1)$day)
}

# Uses julian date functions and function nday.in.month.jd
sst.time.int<-function(jd,dir_sstm,dir_ssth,hr=12)
{
  # convert julian date to day/month/year
  day<-DMYjd(jd)$day ; month<-DMYjd(jd)$month ; year<-DMYjd(jd)$year
  
  # How many days in month
  dimth<-nday.in.month.jd(jd)
  
  # Read in SSTs for the two months that lie either side of the date in question
  # reads in month before if date is in first half of month - using jd allows for month being in different year
  mth1<-DMYjd(jd-(dimth/2))$month 
  mth2<-DMYjd(jd+(dimth/2))$month 
  
  filein1<-paste(dir_sstm,"sst_",year,"_",mth1,".tif",sep="")
  filein2<-paste(dir_sstm,"sst_",year,"_",mth2,".tif",sep="")
  r1<-raster(filein1)
  r2<-raster(filein2)
  v1<-getValues(r1,format="matrix")
  v2<-getValues(r2,format="matrix")
  
  # work out weighting to attach
  
  # How many days in both months
  dimth1<-nday.in.month.jd(jd-(dimth/2))
  dimth2<-nday.in.month.jd(jd+(dimth/2))
  #print(paste(day,"/",month, ".  dimth1= ",dimth1," dimth2=", dimth2,sep=""))
  
  # time after dimth1
  tt<-day+hr/24+dimth1/2
  if (day+hr/24>dimth/2) tt<-day+hr/24-dimth1/2
  wgt<-tt/((dimth1+dimth2)/2)
  v<-v2*wgt+v1*(1-wgt)
  ssthr.r<-raster(v,template=r1)
  
  # OPTION - downscale to 100m dembuf and mask
  #ssthr.r<-resample(ssthr.r,dembuf)
  #ssthr.r<-mask(sst.test,dembuf,inverse=TRUE)
  
  # write daily file
  tl<-paste("Year: ",year," Month: ",month," Day: ",day," Hour: ",hr,sep="")
  #plot(ssthr.r,main=tl)
  fileout<-paste(dir_ssth,"sst_",year,"_",month,"_",day,"_",hr,"h.tif",sep="")
  print(fileout)
  writeRaster(ssthr.r,file=fileout,overwrite=T)
  
} # end function

##########################################################################################

# CODE FOR READING WIND DATA AND RECORDING AS FILES

######################################################################################
# # This bit of the programme converts netcdf files to a single array and stores array
######################################################################################
# define file time groupings
vars<-c(1960,1970,1980,1990,2000,2010)
filenum<-6

# u wind
# This bit works out what size the array should be
store<-0
for (i in 1:filenum){
  infile.u<-paste(dir_wind,"u_",vars[i],".nc",sep="")
  ncdf_windu<-nc_open(infile.u) # summary of file variables,  dimensions, attributes
  wind.u<-ncvar_get(ncdf_windu) # read data
  d<-dim(wind.u) # dimension lengths of wind.u
  store[i]<-d[3] 
}
wind_u<-array(0,dim=c(d[1],d[2],sum(store))) # defines array size using sum of time and d1 and d2 of last file
# This bit stores the data in an array
mn<-1
for (i in 1:filenum)
{
  infile.u<-paste(dir_wind,"u_",vars[i],".nc",sep="")
  ncdf_windu<-nc_open(infile.u)
  wind.u<-ncvar_get(ncdf_windu)
  mx<-mn+store[i]-1
  tp<-paste("min=",mn," max=",mx,sep="")
  print(tp)
  wind_u[,,mn:mx]<-wind.u
  mn<-mn+store[i]
}
save(wind_u,file=paste(dir_wind,"wind_u.r",sep=""))

# v wind
# This bit works out what size the array should be
store<-0
for (i in 1:filenum){
  infile.v<-paste(dir_wind,"v_",vars[i],".nc",sep="")
  ncdf_windv<-nc_open(infile.v)
  wind.v<-ncvar_get(ncdf_windv)
  d<-dim(wind.v)
  store[i]<-d[3]
}
wind_v<-array(0,dim=c(d[1],d[2],sum(store)))
# This bit stores the data in an array
mn<-1
for (i in 1:filenum)
{
  infile.v<-paste(dir_wind,"v_",vars[i],".nc",sep="")
  ncdf_windv<-nc_open(infile.v)
  wind.v<-ncvar_get(ncdf_windv)
  mx<-mn+store[i]-1
  tp<-paste("min=",mn," max=",mx,sep="")
  print(tp)
  wind_v[,,mn:mx]<-wind.v
  mn<-mn+store[i]
}
save(wind_v,file=paste(dir_wind,"wind_v.r",sep=""))
