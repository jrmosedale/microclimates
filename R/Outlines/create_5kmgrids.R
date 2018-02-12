
# Set t5km.r to a raster of 5km historic temperature data cropped to correct dem extent
# Could use UKCP09 gridsw.r - created below 
# Assumes already rounded to 5km resolution
t5km.r<-hr.r

# METHOD 1 - produces too many NA 5km cells - NO!
# Create x y vectors of central point (Easting,Northing) for each grid cell
# for which there is historic temperature data
xy<-coordinates(t5km.r) # extracts coordinates for centre point of each cell
dem5km<-extract(dem,xy,na.rm=FALSE, method='simple') # creates grid cell vector of cental dem values
# dem5km<-ifelse(is.na(dem5km),NA,1) # To create outline raster only
dem5km.m<-matrix(dem5km,nrow=NROW(t5km.r),ncol=NCOL(t5km.r),byrow=TRUE)
dem5km.r<-raster(dem5km.m,template=t5km.r)
plot(dem5km.r)

# METHOD 1b - bilinear interpolation - NO!
dem5km<-extract(dem,xy,na.rm=FALSE, buffer=500, fun=mean) # creates grid cell vector of cental dem values
# dem5km<-ifelse(is.na(dem5km),NA,1) # To create outline raster only
dem5km.m<-matrix(dem5km,nrow=NROW(t5km.r),ncol=NCOL(t5km.r),byrow=TRUE)
dem5km.r<-raster(dem5km.m,template=t5km.r)
plot(dem5km.r)

# METHOD 2. # Calculate mean elevation for 5km grid cells matching temperature 5km grid
dem5km<- aggregate(dem,fact=50,fun=mean)
# Set sea cells to NA if NA in temperature 5km 
t5km<-getValues(t5km.r)
dem5km<-getValues(dem5km)
dem5km<-ifelse( is.na(t5km), NA , dem5km)
dem5km.m<-matrix(dem5km,nrow=NROW(t5km.r),ncol=NCOL(t5km.r),byrow=TRUE)
dem5km.r<-raster(dem5km.m,template=t5km.r)
plot(dem5km.r)
compareRaster(dem5km.r,t5km.r,values=FALSE)

###############################################################################
# Create 5km matrix containing Weather Generator grid cell ID values
# NB: historic 5km temp files exclude certain cells (where mid point is not land?)
# Input file from: http://ukclimateprojections-ui.metoffice.gov.uk/ui/docs/grids/wg_5km/index.php
###############################################################################
dir_grids<-"~/Documents/Exeter/Data2015/Templates/"
#dir_grids<-"C:/Data2015/Templates/"

file.in<-paste(dir_grids,"grid_box_ids_5km.csv",sep="")
#file.in<-paste(dir_grids,"grid_box_ids_with_mask.csv",sep="")

print(file.in)
grid.m<-as.matrix(read.csv(file=file.in,header=FALSE))

# Convert -9999 values to NA if required
sel<-which(grid.m==-9999)
grid.m[sel]<-NA

# Create raster using easting/northings for bottom left of grid cells (centre=+2500)
grid.5km<-raster(nrows=290,ncols=180,xmn=-200000, ymn=-200000, xmx=700000,ymx=1250000,res=c(5000,5000), crs="+init=epsg:27700")
#grid.5km<-raster(nrows=290,ncols=180,xmn=-197500, ymn=-197500, xmx=697500,ymx=1247500,res=c(5000,5000), crs="+init=epsg:27700")
gridid.r<-raster(grid.m,template=grid.5km)
file.out<-paste(dir_grids,"ukcp_gridcells.r",sep="")
ukcp_gridcells<-getValues(gridid.r)
write(ukcp_gridcells,file.out)

sel<-which(!is.na(ukcp_gridcells))
ukcp_gridmask<-ukcp_gridcells
ukcp_gridmask[sel]<-1
ukcp_gridmask.r<-setValues(gridid.r,ukcp_gridmask)
plot(ukcp_gridmask.r)

file.out<-paste(dir_grids,"ukcpmask.grd",sep="")
writeRaster(ukcp_gridmask.r,file.out,overwrite=TRUE)
         
#gridsw.r<-crop(x=gridid.r,y=dem) # crop to geographical extent of DEM raster
#centrexy<-xyFromCell(gridsw.r,1:ncell(gridsw.r)) # centre xy coordinates for 5km grid cells

###############################################################################
# Create historic grid cell raster mask where cells for which no data=NA
###############################################################################
dir_temp<-"~/Documents/Exeter/Data2015/Temp5km/extract/"
#dir_temp<-"C:/Data2015/Temp5km/"

max.infile<-paste(dir_temp,"MaxTemp_", 2011, "-",sprintf("%02d",6,sep=""),"-", sprintf("%02d",10,sep=""),"_ACTUAL.txt", sep="")
ukhist.r<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
ukhist_gridmask<-getValues(ukhist.r)
sel<-which(!is.na(ukhist_gridmask))
ukhist_gridmask[sel]<-1
ukhist.r<-setValues(ukhist.r,ukhist_gridmask)
plot(ukhist.r)

file.out<-paste(dir_grids,"ukhistmask.grd",sep="")
writeRaster(ukhist.r,file.out)

###############################################################################
# Downscale t5km.r to 100m accounting for difference in altitude
# INPUT: 5km dem and 100m dem raster
# Convert 5km dem back to 100m to store 5km grid cell averages at 100m resolution

# Treat 5km temperatures and 5km dem as point data x,y,z???

# define 100m temperature raster
t100.r<-raster(nrows=NROW(dem),ncols=NCOL(dem),res=res(dem), 
               xmn=xmin(dem), xmx=xmax(dem), ymn=ymin(dem),ymx=ymax(dem), 
                vals=NULL)

# LOOP here for each time period
# METHOD 1 convert 5km t to sea level temperatures, interpolate, reconvert for 100m altitude
# Temperature declines by 9.8C per 100m ??? CHECK ???

# Calculate temperature at sea level for each 5km cell
t5km.sl.r<-raster(nrows=NROW(dem5km.r),ncols=NCOL(dem5km.r),res=res(dem5km.r), 
                  xmn=xmin(dem5km.r), xmx=xmax(dem5km.r), ymn=ymin(dem5km.r),ymx=ymax(dem5km.r), vals=NULL)
vals<-getValues(t5km.r+((dem5km/1000)*9.8))
t5km.sl.r<-setValues(t5km.sl.r,vals)

#t5km.sl.r<-t5km.r+((dem5km/1000)*9.8)

# Interpolate (RESAMPLE) RASTER to 100m cells
t100.sl.r<-resample(t5km.sl.r,dem, nrow=NROW(dem),ncol=NCOL(dem),method='bilinear')

# Convert 100m sea level temperature to actual from dem
vals<-getValues(t100.sl.r-((dem/1000)*9.8))
t100.r<-setValues(t100.r,vals)
# plot(t100.r)




###############################################################################
#CUT OUTS
###############################################################################

#### CREATE matrices if needed #####
# Create matrix of 5km dem
dem5km.v<-getValues(dem5km)
dem5km.m<-matrix(dem5km.v,nrow=NROW(dem5km),ncol=NCOL(dem5km),byrow=TRUE) # convert array to matrix

# Set temperature data at 5km
#t5km.r<-READ HOURLY TEMPERATURE DATA as matrix
t5km.v<-getValues(t5km.r)
t5km.m<-matrix(t5km.v,nrow=NROW(dem5km),ncol=NCOL(dem5km),byrow=TRUE) # convert array to matrix
###############################################################################
