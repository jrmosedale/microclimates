# Works out the angle to the horizon in a specified direction (used to calculate the shelter coefficient)
# Inputs:
# dtm = a digital eleveation model stored as a matrix
# Outputs:
# 360 files - map of shelter coefs for each wind direction

library(ncdf4)
library(raster)
library(rgdal)

# NB the rotation of the digital elevetation data is important. This is designed to be used for a matrix
# extracted from a raster (see raster package) as follows: my.matrix<-getValues(my.raster,format="matrix")
horizonangle <- function(dtm,azimuth,res=100,steps=40)
{
  azimuth<-azimuth-90
  azi <- azimuth * (pi/180)
  horizon <- array(0,dim(dtm))
  dtm3 <- array(0,dim(dtm)+200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)] <- dtm
  m<-10^2/steps^2
  for (step in 1:steps) {
    horizon[1:x,1:y] <- pmax(horizon[1:x,1:y], (dtm3[(101+sin(azi)*m*step^2):(x+100+sin(azi)*m*step^2),(101+cos(azi)*m*step^2):(y+100+cos(azi)*m*step^2)]-dtm3[101:(x+100),101:(y+100)])/(m*res*step^2))
  }
  horizon
}

windindex <- function(dtm,direction)
{
  index <- 1 - atan(0.17 * 100 * horizonangle(dtm,direction))/1.65
  index
}

########################################################################################
#dem<-raster("C:/Data2015/DEM100/dem_sw_x60-420k_y-10-180k.tif")
dem<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif", crs=("+init=epsg:27700"))
#plot(dem,main="DEM-full")
extent(dem)
#e.dem<-extent(c(70000,420000,0,180000)) # includes scilly isles
e.dem<-extent(c(120000,420000,0,180000)) # excludes scilly isles
dem<-crop(dem,e.dem)
plot(dem,main="DEM-sw")
e.dem <-extent(dem)

########################################################################################
#dir_wind<-"C:/Data2015/Wind/"
dir_wind<-"~/Documents/Exeter/Data2015/Wind/"
dir_out<-"~/Documents/Exeter/Data2015/Wind/Shelter/"

# Go through in 10 km blocks and adjust by shelter coefficient
# note, however that actually it actually selects 30km x 30 km area to allow for sheltering effects
# that operate outside the area of each block. The 10km x 10km centre of the block is then selected
# Programme could probably be speeded up without major loss of accuracy by setting buffer to ~5km instead of 10km
# NB chopping into 10km blocks is necessary, as the function for calculating
# the shelter coefficient assumes a single wind direction, a fairly safe assumption over 10km, but not over entire study region

# 1. Create matrix from dem
# ASSUMES dem already cropped to match 10km cells 
# Convert dem matrix
m.dem<-getValues(dem,format="matrix")

# Create 3d array to hold shelter coeeficient values shelter(1:360,nrow(m.dem),ncol(m.dem))
#shelter<-array(NA, 360,nrow(m.dem),ncol(m.dem)))

# Goes through and does each 10km block (original: rws in 1:16, cls in 1:28) 
#  -2 excludes outside cells for which bufer cannot be calculated
# set buffer.km to 10km
buffer.km<-10
buffer<-buffer.km*10
mxrws<-nrow(m.dem)/100-2
mxcls<-ncol(m.dem)/100-2


# 2. For every wind direction (1-360 degrees in steps of 5) calculate a shelter coefficient map
for (d in 1:360)
        {
        # creates 2D matrix for storing wind coeff output 
        w.coef<-matrix(NA,nrow=(nrow(m.dem)),ncol=(ncol(m.dem)) )# matrix for storing output values
        sea<-FALSE
        for (rws in 1:mxrws) 
        {
          for (cls in 1:mxcls)
          { 
            xmn<-rws*100+1-100
            ymn<-cls*100+1-100
            xmx=xmn+100+(2*buffer)-1 
            ymx=ymn+100+(2*buffer)-1 
            b.dem<-m.dem[xmn:xmx,ymn:ymx]
            # Check if dem block is all sea
            sel<-which(is.na(b.dem)==T)
            if (length(sel)==length(b.dem)) sea==TRUE
            b.dem[sel]<-0
            # create block matrix for wind direction and wind coef
            b.wcoef<-matrix(NA,nrow=100+(2*buffer),ncol=100+(2*buffer))
            # Calculates shelter coefficient if wind direction not NA.
            # Wind direction would be NA if all values within 10km block are NA, which happens if the entire 10km block is sea
            if (sea==FALSE) b.wcoef<-windindex(b.dem,d)
            # selects data for just the 10km x 10km centre of each 30km x 30km block
            w.coef[(xmn+buffer):(xmx-buffer),(ymn+buffer):(ymx-buffer)]<-b.wcoef[(buffer+1):(buffer+100),(buffer+1):(buffer+100)]
            #print(paste("Row ",rws,", Col ",cls))
          }
        }
        # set sea cells to NA on basis of dem
        w.coef[which(is.na(m.dem))]<-NA
        # save to single 3d array - NOT USED TO SPEED UP
        #shelter[d,,]<-w.coef
        # Convert to raster
        #r<-raster(w.coef,template=dem)
        #plot (r,main=paste("Shelter index map: ",d,sep=""))
        
        # Output wind coef file for each direction
        out.file<-paste(dir_out,"Shelter_",sprintf("%03d",d,sep=""),"_deg.r",sep="")
        print (out.file)
        save(w.coef,file=out.file)
  
} # end direction loop

# Output wind coef file for each direction - NOT USED 
#out.file<-paste(dir_out,"shelter_indexes.r",sep="")
#print (out.file)
#save(shelter,file=out.file)

