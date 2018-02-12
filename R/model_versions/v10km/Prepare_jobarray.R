# SET-UP for job-array of model runs 
# Performs calculations and writes files used for all model runs within a defined area
# For example dem and other maps

# Set root directory according to whether run on Mac, PC , Carson
root<-"/home/ISAD/jm622/Data2015/"   # Source data and output data
in.root<-"/data/jm622/ModelData/" # model input data 
#root<-"~/Documents/Exeter/Data2015/"; in.root<-"~/Documents/Exeter/Data2015/"
#root<-"C:/Data2015/"
print(paste("Input root= ",in.root,sep=""))
print(paste("Output root= ",root,sep=""))
print(paste("Working directory= ",getwd(),sep=""))

# Define DIRECTORIES required
dir_dem<-paste(in.root,"DEM/",sep="")
dir_counties<-paste(in.root,"counties/",sep="")

# Define LIBRARIES
print("Defining libraries")
library(RNetCDF)
library(ncdf4)
library(raster)
library(rgdal)
library(sp)
library(insol) # required for julian day functions
library(mgcv) # require package for inputation (of CAL)
library(fields) # required for thin plate spline

# Write DEMs required from dem (100m) of UK
#  DEM file locations
ukdem.file<-paste(dir_dem,"demoriginal.tif",sep="")
dem.file<-paste(dir_dem,"dem.tif",sep="")
dembuf.file<- paste(dir_dem,"dembuf.tif",sep="")
demland.file<-paste(dir_dem,"demland.tif",sep="")

# Define AREA OF INTEREST - divisible by 5, 10 & 20km
e.dem<-extent(c(60000,420000,0,180000 )) # includes scilly isles#
#e.dem<-extent(c( 120000,420000,0,180000 )) # excludes scilly isles
#e.ios<-extent(c(80000,100000,0,20000 ))

latlong <- "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"

# Define DEM of interest and DEM buffered from UK dem
demuk<-raster(ukdem.file)
projection(demuk)<-"+init=epsg:27700"

dem<-crop(demuk,e.dem)
plot(dem,main="dem")
#writeRaster(dem,file=dem.file)

# Define  buffered area around region of interest
buffer<-20000
e.bufkm<-extent(xmin(dem)-buffer,xmax(dem)+buffer,ymin(dem)-buffer,ymax(dem)+buffer)# Run setup programs for creating constant raster maps etc
dembuf<-crop(demuk,e.bufkm)
plot(dembuf,main="dembuf")
#writeRaster(dembuf,file=dembuf.file)

### Calculate larger Job_ARRAY grid cells containing 100m land cells from area of interest
# Set GridCell Size eg 5, 10 20 km in metres
gridcell<-20000

# Calculate number of 100m land cells in gridcell 
land100m.r<-calc(dem,function(x) ifelse(is.na(x),NA,1))
landgridcell.r<-aggregate(land100m.r,fact=(gridcell/100),fun=sum,expand=TRUE)
landgridcell.r<-calc(landgridcell.r,function(x) ifelse(is.na(x),NA,1))
plot(landgridcell.r)

vals<-getValues(landgridcell.r)
xy<-xyFromCell(landgridcell.r,1:ncell(landgridcell.r))
sel<-which(vals==1)

cells.df<-as.data.frame(xy[sel,1:2])   # = coordinates for middle of each cell
cells.df$xmin<-cells.df$x-(gridcell/2)
cells.df$xmax<-cells.df$x+(gridcell/2)
cells.df$ymin<-cells.df$y-(gridcell/2)
cells.df$ymax<-cells.df$y+(gridcell/2)

# Sort so 1st gridcell in SW corner
cells.df<-cells.df[order(cells.df$x, cells.df$y),]

# Assign cell number 
cells.df$cellnum<-c(1:nrow(cells.df))

# Print dems of eaach gridcell
for (n in 1:nrow(cells.df)){
  plot(crop(dem,extent(cells.df$xmin[n],cells.df$xmax[n],cells.df$ymin[n],cells.df$ymax[n])),main=paste("Grid cell ",n))
}


#### Option - reduce area to only part of total ###
# Assign whole sw dem to demsw
demsw<-crop(demuk,sw.dem)
projection(demsw)<-"+init=epsg:27700"


# Choose area - divisible by 20km - to align with grid same as if whole area used
area.e<-extent(c(60000,260000,0,120000 )) # cornwall & IoS
# area.e<-extent(c(260000,320000,0,160000 ))  # W Devon
# area.e<-extent(c(320000,420000,0,180000 )) # Rest
dem<-crop(demuk,area.e)
plot(dem,main="dem") 
e.bufkm<-extent(xmin(dem)-buffer,xmax(dem)+buffer,ymin(dem)-buffer,ymax(dem)+buffer)# Run setup programs for creating constant raster maps etc
dembuf<-crop(demuk,e.bufkm)
plot(dembuf,main="dembuf") 

# Keep only grid cells that fall within area
in_area<-extract(dem, cells.df[,1:2], cellnumbers=TRUE) 
sel<-which(!is.na(in_area[,1]))
cells.df<-cells.df[sel,]  # keep only those cells in map

### Write grid cell coordinates file ###
gridcell.file<-paste(dir_dem,"gridcellxy",sep="")
write.csv(cells.df,file=gridcell.file)

### 





# Define areas covered by Met Office 5km grids ??????

# Define 5km grid rasters for which there is historical Met Off data
grid5kmuk.r<-raster(paste(dir_grids,"ukhistmask.grd",sep="")) #  1=valid cell, NA = sea or not data
projection(grid5kmuk.r)<-"+init=epsg:27700"
grid5km.r<-crop(grid5kmuk.r,e.dem)
grid5kmbuf.r<-crop(grid5kmuk.r,dembuf)
plot(grid5km.r,main="grid5km.r")
plot(grid5kmbuf.r,main="grid5kmbuf.r")
#writeRaster(grid5km.r,file=grid5km.file,overwrite=TRUE)
#writeRaster(grid5kmbuf.r,file=grid5kmbuf.file,overwrite=TRUE)


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


demland.r<-raster(demland.file)


remove(demuk)
