### INPUTS required for model run - apply2



##########################################################################################
# Prepare data required for UKCP cell fixed for all timeperiods
##########################################################################################
print(paste("UKCP cell= ", ukcpcell,sep=""))
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
shelter.block<-block.sheltermap(dem.block,dir_shelter,interval)

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
print("Calculating ldif stack...")
for (direction in seq(0,360,interval)) {
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
