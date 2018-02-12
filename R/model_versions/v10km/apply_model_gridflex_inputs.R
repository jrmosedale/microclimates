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

##########################################################################################
# Prepare data required for grid cell fixed for all timeperiods
##########################################################################################
gridcell.file<-paste(dir_dem,"gridcellxy",sep="")
cells.df<-read.csv(file=gridcell.file)

print(paste("Grid cell= ", cellnum,sep=""))
print(paste("Coordinates of grid cell =",cells.df[cellnum,"xmin"],cells.df[cellnum,"xmax"],cells.df[cellnum,"ymin"],cells.df[cellnum,"ymax"]))
# Crop dem to fit cell and buffered cell
e.block<-extent(cells.df[cellnum,"xmin"],cells.df[cellnum,"xmax"],
                cells.df[cellnum,"ymin"],cells.df[cellnum,"ymax"])
dem.block<-crop(dem,e.block) 

# size of buffer region remains same using larger 10km cells
e.buffer<-extent(cells.df[cellnum,"xmin"]-buffer,cells.df[cellnum,"xmax"]+buffer,
                 cells.df[cellnum,"ymin"]-buffer,cells.df[cellnum,"ymax"]+buffer)
dem.buffer<-crop(dembuf,e.buffer)

plot(dem.buffer,main=paste("DEM.BUFFER ",cellnum))
plot(dem.block,main=paste("DEM.BLOCK ",cellnum))

# Perform operations for each block but fixed for all timeperiods - ie CROPPING Terrain datasets
# Load shelter matrix for block
shelter.block<-block.sheltermap(dem.block,dir_shelter,interval)

# Slope & aspect cropped to dem.buffer for use in radiation_downscale
slope.r<-raster(paste(dir_terrain,"slope.tif",sep=""))
aspect.r<-raster(paste(dir_terrain,"aspect.tif",sep=""))
projection(slope.r)<-CRS("+init=epsg:27700")
projection(aspect.r)<-CRS("+init=epsg:27700")
slope.buffer<-crop(slope.r,dem.buffer)   # Are these needed??
aspect.buffer<-crop(aspect.r,dem.buffer)

# Calculate elevation difference from 5km mean
elevdif.r<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # created by elevation.dif.map function
projection(elevdif.r)<-CRS("+init=epsg:27700")
elevdif.block<-crop(elevdif.r,dem.block)

# Load and extract ldif for each wind direction for block - used with hourly wind direction to calc coastal effect - write it as temp file?
print("Calculating ldif stack...")
for (direction in seq(0,360,interval)) {
  # ldif.r<-raster(paste(dir_ldif,"ldif_",direction,"deg_from_percland_in_",radius/1000,"km.tif",sep=""))
  ldif.r<-raster(paste(dir_ldif,"ldif_v2_",direction,"deg_from_percland_in_",radius/1000,"km.tif",sep=""))
  projection(ldif.r)<-CRS("+init=epsg:27700")
  ldif.layer<-crop(ldif.r,dem.block)
  if (direction==0) ldif.stack<-stack(ldif.layer) else ldif.stack<-stack(ldif.stack,ldif.layer)
} # end for direction

# Extract ALBEDO for area and set NA to default value
albedomap.r<-raster(paste(dir_albedo,"albedo_mean_dembuf.tif",sep=""))
projection(albedomap.r)<-CRS("+init=epsg:27700")
albedo.block<-crop(albedomap.r,dem.block)

# Extract TWI and Alt difference for area
twi.block<-crop(raster(paste(dir.basinmap,"topidx.tif",sep="")),dem.block)
altdif.block<-crop(raster(paste(dir.basinmap,"altdif.tif",sep="")),dem.block)

# Remove larger raster not needed
remove(albedomap.r,ldif.r,ldif.layer,elevdif.r,slope.r,aspect.r)
