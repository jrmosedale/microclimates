library(ncdf4)
library(raster)
library(rgdal)

####################################################
# Function for computing Julian data 
# NOTE: JD reconverts to one day later - NOT a problem as only used for ID of correct array location
# In this instance just used to work out number of days after 1st Jan 1960 (1st wind data observation),
# so that correct element of array can be extracted
# Inputs:
  # day: the day of the month (0-31)
  # month: the month of the year (numeric: 0-12)
  # year: any year from 1950 to 2014
# Output: the Julian day (https://en.wikipedia.org/wiki/Julian_day)
  JD<-function(day,month,year){
    a<-(14-month)/12
    y<-year+4800-a
    m<-month+12*a-3
    JDN<-floor((153*m+2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045 + day
    JDN
  }
# Based on hour, day, month and year, extracts the required value from the array of values stored by wind_downscale1.R
# Inputs:
  # hr: the hour (0-23) 0 = midnight
  # day: the day of the month (0-31)
  # month: the month of the year (numeric: 0-12)
  # year: any year from 1960 to 2014
# Output: the element of the array stored by wind_downscale1.R that corresponds to either that hour, or the latest
# period immediatly before that hour (data only available 6-hourly)
array.val<-function(hr,day,month,yr)
{
 jd.base=JD(1,1,1960)
 jd<-JD(day,month,yr)
 dval<-(jd-jd.base)*4
 hval<-floor(hr/6)
 val<-dval+hval+1
 val
}
# Works out the angle to the horizon in a specified direction (used to calculate the shelter coefficient)
# Inputs:
  # dtm = a digital eleveation model stored as a matrix

####################################################
# Read in Digital Eelevation data - DECISION: what area to cover
# **** ENSURE here that dem is rounded to 10km ****

#dem<-raster("C:/Data2015/DEM100/dem_sw_x60-420k_y-10-180k.tif")
dem<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif", crs=("+init=epsg:27700"))
#plot(dem,main="DEM-full")
extent(dem)
#e.dem<-extent(c(70000,420000,0,180000)) # includes scilly isles
e.dem<-extent(c(120000,420000,0,180000)) # excludes scilly isles
dem<-crop(dem,e.dem)
plot(dem,main="DEM-sw")
e.dem <-extent(dem)
# Create dem matrix
m.dem<-getValues(dem,format="matrix")

####################################################
# INPUT:  interval to which wind direction is rounded - eg 5 = every 5 degrees
interval<-5
####################################################
# Load shelter maps for each degree each saved as R name w.coef - takes about 1hr40
# Create 3d array to hold shelter coeeficient values shelter(1:360,nrow(m.dem),ncol(m.dem))
dir_shelter<-"~/Documents/Exeter/Data2015/Wind/Shelter/"
shelter<-array(NA, dim=c((360/interval),nrow(m.dem),ncol(m.dem)))
for (i in 1:(360/interval)) {
  dir<-i*interval
  in.file<-paste(dir_shelter,"Shelter_",sprintf("%03d",dir,sep=""),"_deg.r",sep="")
  load(file=in.file) 
  shelter[i,,]<-w.coef
  print(paste("i= ",i," dir= ",dir))
}

# CHECK SHELTER MAPS
# Plot E/S/W/N direction indices
par(mfrow=c(2,2))
brk<-c(seq(0,1,0.1))
col<-rev(terrain.colors(11))

for (i in c((90/interval),(180/interval),(270/interval),(360/interval))){
r<-raster(shelter[i,,],template=dem)
plot(r,main=paste("Shelter coef - wind from ", (i*interval), " degrees", sep=""), col=col, breaks=brk)
}     
# calculatae mean index for each direction
#apply(shelter,1,mean,na.rm=TRUE)

####################################################
# Load wind data (single file for whole time period - created by wind_downscale1)
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v

#dir_wind<-"C:/Data2015/Wind/"
dir_wind<-"~/Documents/Exeter/Data2015/Wind/"
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))

####################################################
# # # This bit downscales the wind
####################################################
# NB running this for 1 month takes ~20 hours
# # # Stages:
# (1) get wind values for a given hour, day, month and year
# (2) convert to 100m resolution raster OSGB grid reference
# (3) adjust based on altitude
# (4) adjust based on shelter coefficient

# Specify month and year for which data are required
yr=2014
month=6
day<-10


for (day in 10:11){
  ptm <- proc.time()
  
      for (hr in 0:23){
            # Can be quite slow. Allows you to keep tabs on progress by printing hour, day, month & year
            tp<-paste("year=",yr," month=",month," day=",day," hour=",hr,sep="")
            print(tp)
            
            #############
            # Stage 1: get wind values for a given day month and year
            #############
            # As original data are 4x daily, but data are required for each hour,
            # this bit reads in the data for the periods immediatly before after for which there are data and calculates
            # weighted mean
            av1<-array.val(hr,day,month,yr)
            av2<-av1+1
            rem<-hr/6-floor(hr/6)
            uwind1<-wind_u[,,av1]
            uwind2<-wind_u[,,av2]
            vwind1<-wind_v[,,av1]
            vwind2<-wind_v[,,av2]
            uwind<-(1-rem)*uwind1+rem*uwind2
            vwind<-(1-rem)*vwind1+rem*vwind2
            
            #############
            # Stage 2: convert to 100m resolution raster OSGB grid reference
            #############
            # Convert to raster (original lat long format and resolution - CHECK LAT/LON MAX/MIN
            uwind.r<-raster(uwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
            vwind.r<-raster(vwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
            # Reproject in OSGB projection
            crs(uwind.r)<-"+init=epsg:4326"
            crs(vwind.r)<-"+init=epsg:4326"
            u_osgb<-projectRaster(uwind.r,crs="+init=epsg:27700")
            v_osgb<-projectRaster(vwind.r,crs="+init=epsg:27700")
            # Trim to desired area and resample the data at a 100m resolution. Interpolation is set at default simple bilinear
            e<-extent(u_osgb) 
            xmn=round(e@xmin/100)*100
            xmx=round(e@xmax/100)*100
            ymn=round(e@ymin/100)*100
            ymx=round(e@ymax/100)*100
            template<-raster(xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,resolution=100)
            u_100<-resample(u_osgb,template)
            v_100<-resample(v_osgb,template)
            
            # Once resampled, can crop wind data to match DEM raster
            u_100_b<-crop(u_100,e.dem) # define u_100_b to avoid re-running resampling in testing
            v_100_b<-crop(v_100,e.dem)
            
            #############
            # Stage 3: adjust based on altitude of terrain
            # NB Height adjustment based on wind spped values at different pressures downloaded  Earth System Research Lab
            # Typical heights at different pressures calculated from Allen et al 1998 http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric pressure (p)
            # Quadratic function fitted - NB this works well for heights up to ~1800m. IT won't work above ~2000m
            # Function was first derived by comparing values at different pressures (heights) over the course of a year (2014)
            #############
            # adjust wind speeds by height of dem
             # convert to matrices
            uwind.m<-getValues(u_100_b,format="matrix")
            vwind.m<-getValues(v_100_b,format="matrix")
            dem.m<-getValues(dem,format="matrix")
             # adjust wind by height
            ustr<-sqrt(uwind.m^2) # wind strength
            vstr<-sqrt(vwind.m^2) # wind strength
            udir<-ifelse(uwind.m>0,1,-1) # positive or negative
            vdir<-ifelse(vwind.m>0,1,-1) # positive or negative
            u.adj<-ustr*((-0.000000108025)*dem.m^2+0.000408692*dem.m+0.956139) # NB don't worry about warnings. Calculation assigns NAs to the sea
            v.adj<-vstr*((-0.000000108025)*dem.m^2+0.000408692*dem.m+0.956139)  # NB don't worry about warnings. Calculation assigns NAs to the sea
            # adjust values to correspond to wind speed 1m above theground
            # rescaling factor first derived by comparing values ot Culdrose wind data using wind_downscale3
            # however, in line wiht what you'd expect from: http://www.fao.org/docrep/x0490e/x0490e07.htm#wind profile relationship
            u.adj<-u.adj*0.373686439
            v.adj<-v.adj*0.373686439
            u.adj<-u.adj*udir
            v.adj<-v.adj*vdir
            ru<-raster(u.adj,template=dem)
            rv<-raster(v.adj,template=dem)
            
            # Some code here for plotting altitude adjusted values, currently commmented out
            #par(mfrow=c(1,1))
            #plot(ru,main="altitude adjusted wind u")
            #plot(rv,main="altitude adjusted wind v")
            
            # calculate wind direction
            direction = (180/pi)*(atan2(u.adj,v.adj))  # NB this is direction in which wind blows to
            direction<-ifelse(direction<=180,direction+180,direction-180) # NB this direction from which wind originates (360 deg)
            rd<-raster(direction,template=dem)
            #plot(rd,main="wind direction")
            
            #############
            # Stage 4: height adjustments done using shelter coefficient maps based on topography and wind direction
            # **** ASSUMES: shelter maps calculated using identical dem ****
            #############
            
            # Extend rasters to match 10km cells
            # BUT does this produce empty cells - better to reduce size to 10
            #e<-extent(ru) 
            #print(e)
            #xmn=ceiling(e@xmin/10000)*10000  # Decide if crop or extend best - whether NA around border desirable?
            #xmx=ceiling(e@xmax/10000)*10000
            #ymn=ceiling(e@ymin/10000)*10000
            #ymx=ceiling(e@ymax/10000)*10000
            #e<-extent(c(xmn,xmx,ymn,ymx))
            #print(e)
            
            #ru<-extend(ru,e)
            #rv<-extend(rv,e)
            #dem<-extend(dem,e)
            #rd<-extend(rd,e)

            # Create matrixes from rasters of wind vectors and direction - ORIENTATION???
            m.u<-getValues(ru,format="matrix")
            m.v<-getValues(rv,format="matrix")
            m.dem<-getValues(dem,format="matrix")
            m.dir<-getValues(rd,format="matrix")
            
            # Uses Rounded wind direction of each cell to select correct shelter coefficient
            # Question:  use aggregate to mean across 100m cells to match block used for shelter calc??
            
            m.dir<-round(m.dir/interval)*interval # round to interval used for shelter maps
            m.dir<-ifelse(m.dir==0,360,m.dir) # converts 0 to 360 degree direction
            
            # TEST using m.dir of 1 or 2 deg direction
            #m.dir<-(rep(c(1,2),length(m.dir/2)))
            #m.dir<-matrix(rep(1,2),nrow=nrow(m.dem), ncol=ncol(m.dem))
            
            # Create windstrength matrix for storing values - original: matrix(NA,nrow=1800,ncol=3000)
            m.str<-matrix(NA,nrow=(nrow(m.dem)),ncol=(ncol(m.dem)) )# matrix for storing all values
            
            # Calculates wind strength from u and v components
            m.str<-sqrt(m.u^2+m.v^2)
            #set m.str to NA if sea (on basis of dem)
            m.str[which(is.na(m.dem))]<-NA
            
            # Applies shelter coefficient to calculate wind strength if land cell (else NA)
            mxrws<-nrow(m.dem)
            mxcls<-ncol(m.dem)
            for (rws in 1:mxrws) {
                  for (cls in 1:mxcls) {    
                        m.str[rws,cls]<-m.str[rws,cls]*shelter[(m.dir[rws,cls]/interval),rws,cls] 
                  }
            }
            
            
            # converts to raster (and crops) if required
            #r<-raster(m.str,template=dem)
            #plot(r,main="wind speed")
            
            # Saves matrix as R dataset
            # seperate datasets saved for each hour,day,month,year (entire study area)
            
            fileout.1<-paste(dir_wind,"strength/strength_",yr,"_",month,"_",day,"_",hr,".r",sep="")
            fileout.2<-paste(dir_wind,"direction/direction_",yr,"_",month,"_",day,"_",hr,".r",sep="")
            save(m.str,file=fileout.1)
            save(m.dir,file=fileout.2)
            print(fileout.1)
            print(fileout.2)
      } # end of hr loop
  proc.time() - ptm
  
} # end of day loop



# 


# Load and print 4 maps per day

for (hr in seq(0,23,6)) {
wstr.in<-paste(dir_wind,"strength_",yr,"_",month,"_",day,"_",hr,".r",sep="")
print(wstr.in)
load(file=wstr.in) # creates m.str
w.str.r<-raster(m.str,template=dem)
plot(w.str.r,main=paste("strength_",yr,"/",month,"/",day," ",hr,":00",sep=""))
 }

hr<-0


# Runtime values for one day loop
#user   system  elapsed 
#1731.551  531.737 2313.040 
