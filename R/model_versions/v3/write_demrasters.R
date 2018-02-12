# Create basic DEM files used
dir_dem<-paste(root,"DEM/",sep="")

# DEM rasters used by programs
dem.infile<-paste(dir_dem,"demoriginal.tif",sep="")

# Names for dem rasters to be created
dem.file<-paste(dir_dem,"dem.tif",sep="")
dembuf.file<- paste(dir_dem,"dembuf.tif",sep="")
grid5km.file<-paste(dir_dem,"grid5km.tif",sep="")
grid5kmbuf.file<-paste(dir_dem,"grid5kmbuf.tif",sep="")
land5km.file<-paste(dir_dem,"land5km.tif",sep="")
demland.file<-paste(dir_dem,"demland.tif",sep="")

demuk<-raster(dem.infile)

latlong <- "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"

# Define area of interest (no buffer)
#e.dem<-extent(c(70000,420000,10000,180000 )) # includes scilly isles
e.dem<-extent(c( 120000,420000,10000,180000 )) # excludes scilly isles

# Define 100m dem rasters
demuk<-raster(dem.infile)
projection(demuk)<-"+init=epsg:27700"
e.ukexp<-c(0,7e+05,-10000,1200000) # expand to allow 20km buffer to south of area of interest - set to sea (NA)
demuk<-extend(demuk,e.ukexp,values=NA)
print(demuk)
dem<-crop(demuk,e.dem)
plot(dem,main="dem")
writeRaster(dem,file=dem.file)

# define  buffered area around region of interest
buffer<-20000
e.buf20km<-extent(xmin(dem)-buffer,xmax(dem)+buffer,ymin(dem)-buffer,ymax(dem)+buffer)# Run setup programs for creating constant raster maps etc
dembuf<-crop(demuk,e.buf20km)
plot(dembuf,main="dembuf")
writeRaster(dembuf,file=dembuf.file)

# Define 5km grid rasters for which there is historical Met Off data
grid5kmuk.r<-raster(paste(dir_grids,"ukhistmask.grd",sep="")) #  1=valid cell, NA = sea or not data
projection(grid5kmuk.r)<-"+init=epsg:27700"
grid5km.r<-crop(grid5kmuk.r,e.dem)
grid5kmbuf.r<-crop(grid5kmuk.r,e.buf20km)
plot(grid5km.r,main="grid5km.r")
plot(grid5kmbuf.r,main="grid5kmbuf.r")
writeRaster(grid5km.r,file=grid5km.file,overwrite=TRUE)
writeRaster(grid5kmbuf.r,file=grid5kmbuf.file,overwrite=TRUE)


# Define 5km grid cells containing 100m land cells
land100m.r<-calc(dem,function(x) ifelse(is.na(x),NA,1))
land5km.r<-raster(extent(grid5km.r),crs=crs(grid5km.r),res=res(grid5km.r))  # empty raster of grid5km.r (historic data cells)
land5km.r<-aggregate(land100m.r,fact=50,fun=sum) # = number of 100m land cells in each 5km cell
land5km.r<-calc(land5km.r,function(x) ifelse(is.na(x),NA,1)) # convert to either 1 or NA 
plot(land5km.r,main="land5km.r")
writeRaster(land5km.r,file=land5km.file,overwrite=TRUE)

# dem raster where land=0, sea=NA
demland.r<-calc(dem,function(x) ifelse(is.na(x),NA,0))
plot(demland.r)
writeRaster(dem.land,file=demland.file,overwrite=TRUE)

