# Works out the angle to the horizon in a specified direction (used to calculate the shelter coefficient)
# Inputs:
# dtm = a digital eleveation model stored as a matrix
# Outputs:
# 360 files - map of shelter coefs for each wind direction - Save as RASTERS .tif format

library(ncdf4)
library(raster)
library(rgdal)

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
#dem<-raster("C:/Data2015/DEM100/demoriginal.tif")
demuk<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif", crs=("+init=epsg:27700"))

#e.dem<-extent(c(70000,420000,0,180000)) # includes scilly isles
#e.dem<-extent(c(120000,420000,0,180000)) # excludes scilly isles
dem<-crop(demuk,e.dem)

# 20km buffer required to prodiuce maps for e.dem
buffer<-10000
e.map<-extent(xmin(e.dem)-buffer,xmax(e.dem)+buffer,ymin(e.dem)-buffer,ymax(e.dem)+buffer)

dem.map<-crop(demuk,e.map)
plot(dem.map,main="DEM-map")


########################################################################################
#dir_wind<-"C:/Data2015/Wind/"
#dir_shelter<-"C:/Data2015/Wind/Shelter/"

dir_wind<-"~/Documents/Exeter/Data2015/Wind/"
dir_shelter<-"~/Documents/Exeter/Data2015/Wind/Shelter/"

# Go through in 10 km blocks and adjust by shelter coefficient
# note, however that actually it actually selects 30km x 30 km area to allow for sheltering effects
# that operate outside the area of each block. The 10km x 10km centre of the block is then selected
# Programme could probably be speeded up without major loss of accuracy by setting buffer to ~5km instead of 10km
# NB chopping into 10km blocks is necessary, as the function for calculating
# the shelter coefficient assumes a single wind direction, a fairly safe assumption over 10km, but not over entire study region

# 1. Create matrix from dem
# ASSUMES dem already cropped to match 10km cells 
# Convert dem matrix
dem.m<-getValues(dem.map,format="matrix")

# Create 3d array to hold shelter coeeficient values shelter(1:360,nrow(m.dem),ncol(m.dem))
#shelter<-array(NA, 360,nrow(dem.m),ncol(dem.m)))

# Goes through and does each 10km block 
#  -2 excludes outside cells for which buffer cannot be calculated
mxrws<-nrow(dem.m)/100-2
mxcls<-ncol(dem.m)/100-2
b.cells<-buffer/100

# 2. For every wind direction (1-360 degrees in steps of 5) calculate a shelter coefficient map
interval<-10
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
          }
        }
        # set sea cells to NA on basis of dem
        wcoef[which(is.na(dem.m))]<-NA
        
        # save to single 3d array - NOT USED TO SPEED UP
        #shelter[d,,]<-w.coef
        # Convert and write as raster
        wcoef.r<-raster(wcoef,template=dem.map)
        wcoef.r<-crop(wcoef.r,e.dem)
        par(mfrow=c(1,1))
        plot (wcoef.r,main=paste("Shelter index map: ",d,sep=""))
        out.file<-paste(dir_shelter,"Shelter_",sprintf("%03d",d,sep=""),"_deg.tif",sep="")
        print (out.file)
        writeRaster(wcoef.r,file=out.file,overwrite=T)
  
} # end direction loop


# Calculate Lref - shelter index for each wind dir for central 5km pt.
# 1. Calculate using 5km historic grid cells

# 2. Calculate for 5km UKCP gridcells


# CUT OUTS or TESTS

# For reading and plotting raster of saved file
d<-0
w.coef<-matrix(NA,nrow=(nrow(m.dem)),ncol=(ncol(m.dem)) )
in.file<-paste(dir_out,"Shelter_",sprintf("%03d",d,sep=""),"_deg.r",sep="")
load(file=in.file)
r<-raster(w.coef,template=dem)
#plot (r,main=paste("Shelter index map: ",d,sep=""))

e.test<-extent(c(250000,280000,30000,60000))
test.r<-crop(r,e.test)
plot(test.r,main=paste("Shelter index map: ",d,sep=""))

# test raster 
vals<-

