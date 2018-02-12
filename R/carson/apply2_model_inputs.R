### INPUTS required for model run - apply2

# block.width<-5000 = ASSUMPTION
buffer<-20000 # already set in setup??

### Load or set MODEL PARAMETERS ###
parameters<-"/home/ISAD/jm622/rscripts/inputs/lizardparams.csv"
#parameters<-"~/Documents/Exeter/Data2015/parameters/testparams.csv"
#parameters<-"F:/Data2015/parameters/lizardparams.csv"
params<-read.csv(parameters)
modify_params<-TRUE

if (modify_params==TRUE){   # change parameter values from those in file 
  params$estimates[1]<- 0.528 #intercept
  params$estimates[2]<- 0.137 #wc lapse
  params$estimates[3]<- 0.8  # shortwave rad^2.5
  params$estimates[4]<- 4.3   # longwave
  params$estimates[5]<- -0.43 # albedo
  params$estimates[6]<- -0.015 #windspeed
  params$estimates[7]<- -1.96 #inv wind speed
  params$estimates[8]<- -0.2 #ldif -0.2, -1.5
  params$estimates[9]<- 0.148 #sst-tref
  params$estimates[10]<- -15 #evapdif
  params$estimates[11]<- 0.06 # condendif
  params$estimates[12]<- -0.0012 # flow.acc
  params$estimates[13]<- -0.805 #tic
  params$estimates[14]<- 3.11 # albedo*rad
  params$estimates[15]<- -0.4 #rad*windspeed
  params$estimates[16]<- -0.59 # longwave*windspeed
  params$estimates[17]<- -0.233  #inv windspeed*ldif -2.2
  params$estimates[18]<-0.012 # wind x sst-tref 0.1
  params$estimates[19]<-0.33 # ldif*sst-tref
  params$estimates[20]<-0.0033 #flowacc*tic
}     
print("Parameters: ")
print(params$estimates)
##########################################################################################
# Prepare data required for ANY cell
##########################################################################################
### Load WIND data into memory
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))

# Link to sea surface pressure file 
p.ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")
nc_open(p.ncfile)

print("Get Cell Coordinates")
#in.file<-paste(dir_grids,"ukcpmask.grd",sep="")
#gridmask.r<-raster(in.file) #  1=valid cell, NA = sea or not data
#gridmask.r<-crop(gridmask.r,dem) 
gridmask.r<-land5km.r # land5km.r defined in setup
vals<-getValues(gridmask.r)
xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
sel<-which(vals==1)
landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell
print(dim(landcells))

##########################################################################################
# Prepare data required for UKCP cell fixed for all timeperiods
##########################################################################################
# Crop dem to fit cell and buffered cell
x<-landcells[ukcpcell,1]
y<-landcells[ukcpcell,2]
e.block<-extent(x-2500,x+2500,y-2500,y+2500)
dem.block<-crop(demuk,e.block) ; print (dem.block)
e.buffer<-extent(x-2500-buffer,x+2500+buffer,y-2500-buffer,y+2500+buffer)
dem.buffer<-crop(demuk,e.buffer)
# plot(dem.buffer,main=ukcpcell);plot(dem.block,main=ukcpcell)
print(paste("Plotting location of cell ",ukcpcell,sep=""))
cell_location(gridmask.r,ukcpcell) # plot map of cell location

# Perform operations for each block but fixed for all timeperiods - ie CROPPING Terrain datasets
# Load shelter matrix for block
shelter.block<-block.sheltermap(dem.block,dir_shelter)

# Slope & aspect cropped to dem.buffer for use in radiation_downscale
slope.r<-raster(paste(dir_terrain,"slope.tif",sep=""))
aspect.r<-raster(paste(dir_terrain,"aspect.tif",sep=""))
projection(slope.r)<-CRS("+init=epsg:27700")
projection(aspect.r)<-CRS("+init=epsg:27700")
slope.buffer<-crop(slope.r,dem.buffer)
aspect.buffer<-crop(aspect.r,dem.buffer)

# Calculate elevation difference from 5km mean
elevdif.r<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # created by elevation.dif.map function
projection(elevdif.r)<-CRS("+init=epsg:27700")
elevdif.block<-crop(elevdif.r,dem.block)
#elevdif.block<-dem.block-cellStats(dem.block, stat='mean', na.rm=TRUE)

# Load and extract ldif for each wind direction for block - used with hourly wind direction to calc coastal effect - write it as temp file?
for (direction in seq(0,350,interval)) {
  ldif.r<-raster(paste(dir_ldif,"ldif_",direction,"deg_from_percland_in_",radius/1000,"km.tif",sep=""))
  projection(ldif.r)<-CRS("+init=epsg:27700")
  ldif.layer<-crop(ldif.r,dem.block)
  if (direction==0) ldif.stack<-stack(ldif.layer) else ldif.stack<-stack(ldif.stack,ldif.layer)
} # end for direction

# Extract ALBEDO for area and set NA to default value
albedomap.r<-raster(paste(dir_albedo,"albedo_mean_dembuf.tif",sep=""))
projection(albedomap.r)<-CRS("+init=epsg:27700")
albedo.block<-crop(albedomap.r,dem.block)

# Extract FLOW ACC for area
flowacc<-raster(paste(dir_flowacc,"flowacc_multi.tif",sep=""),res=100)
projection(flowacc)<-CRS("+init=epsg:27700")
minval<-cellStats(flowacc,min)  
flowacc<-flowacc/minval  # divide by min val (area of single cell?)
flow.block<-crop(flowacc,dem.block)
