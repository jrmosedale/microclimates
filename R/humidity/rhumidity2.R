library(ncdf4)
library(raster)
library(rgdal)

dir_rh<-"C:/Data2015/RelHumidity/"
dir_rh_out<-"C:/Data2015/RelHumidity/swfiles"

###################################################
###                FUNCTIONS
###################################################
# FUNCTION - returns number of days in a month for any year between 1980-2015 (accounts for leap years)
monthdays<-function(year,month){
  y=year-1979 # so 1=1980
  feb.d<-c(29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28,29,28,28,28,
           29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28)
  monthdays<-c(31,feb.d[y],31,30,31,30,31,31,30,31,30,31)
  return(monthdays[month])
}

# FUNCTION for computing Julian data
# = number of days after 1st Jan 1950,
JD<-function(day,month,year){
  a<-(14-month)/12
  y<-year+4800-a
  m<-month+12*a-3
  JDN<-floor((153*m+2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045 + day
  JDN
}

# FUNCTION Based on hour, day, month and year, extracts the required value from the array of values 
# REQUIRES: jd.baseto be set to first observation of dataset ###
# Input: day, month, year
# Output: the element of the array stored that corresponds to either that hour, or the latest
# period immediatly before that hour (data only available 6-hourly)
rh.array.val<-function(hr,day,month,yr)
{
  jd.base=JD(1,1,1980)
  jd<-JD(day,month,yr)
  dval<-(jd-jd.base)*4
  hval<-floor(hr/6)
  val<-dval+hval+1
  val
}

#### Read in Digital Eelevation data - DECISION: what area we want to cover ####
dem<-raster("C:/Data2015/DEM100/dem_sw_x60-420k_y-10-180k.tif")
#dem<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif")
plot(dem,main="DEM")
e<-extent(dem)
# Crop areas so that it neatly has whole 100m grid cells around the edge and save extent as e.dem for use later
xmn=floor(e@xmin/100)*100 # CHANGED from Round to Floor - to ensure cropping
xmx=floor(e@xmax/100)*100
ymn=floor(e@ymin/100)*100
ymx=floor(e@ymax/100)*100
e.dem<-extent(c(xmn,xmx,ymn,ymx))
dem<-crop(dem,e.dem)
plot(dem,main="DEM-original")
e.dem <-extent(dem)
#crs(dem)<-"+init=epsg:4326"
#plot(dem,main="DEM-OScrs")

####################################################
# # # This bit downscales rh                     ###
####################################################
load(file=paste(dir_rh,"rh.r",sep=""))

# Specify month and year for which data are required
yr=2014
month=6
day<-30

# set period for which you want to create 100m resolution data
# Set to do all hours in January 2014
for (day in 1:monthdays(year,month)){
  par(mfrow=c(2,2)) # for printing
  for (hr in 0:23){
    # Can be quite slow. Allows you to keep tabs on progress by printing hour, day, month & year
    print(paste("year=",year," month=",month," day=",day," hour=",hr,sep=""))
        
    #############
    # Stage 1: get  values for a given day month and year
    #############
    # As original data are 4x daily, but data are required for each hour, this bit reads in the data for the periods immediatly before after for which there are data and calculates
    # weighted mean
    av1<-rh.array.val(hr,day,month,yr)
    av2<-av1+1
    rem<-hr/6-floor(hr/6)
    rh1<-rh[,,av1]
    rh2<-rh[,,av2]
    rh.hr<-(1-rem)*rh1+rem*rh2 # = lat/lon rh grid for single hour (of day,month,year) 
    
    #############
    # Stage 2: convert to 100m resolution raster OSGB grid reference
    #############
    # Convert to raster (original lat long format and resolution - CHECK LAT/LON MAX/MIN
    rh.r<-raster(rh.hr,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5) # check with ncdf_rh
    # Reproject in OSGB projection
    crs(rh.r)<-"+init=epsg:4326"
    rh_osgb<-projectRaster(rh.r,crs="+init=epsg:27700")
    
    
    # Trim to desired area and resample the data at a 100m resolution. Interpolation is set at default simple bilinear
    # Derive template from dem
    # Q. perform once on 24 hoursly stacks ?? Memory?
    xmn=round(e.dem@xmin/100)*100
    xmx=round(e.dem@xmax/100)*100
    ymn=round(e.dem@ymin/100)*100
    ymx=round(e.dem@ymax/100)*100
    template<-raster(xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,resolution=100)
    
    rh_100<-resample(rh_osgb,template)

    # Set sea cells (NA in dem) to NA i
    rh_100_land <- overlay(x = dem,y = rh_100,fun = function(x,y) {
      y[is.na(x)] <- NA
      return(y)
    } )   
    
    #####################################################################
    # Plot rh by hour using same scale
    brk<-c(0,10,20,30,40,50,60,70,80,90,100)
    col<-(rainbow(11,start=1/6,end=4/6))
    t<-paste("rh ",day,"-",month,"-",year," Hour: ", hr, ":00",sep="")
    plot(rh_100_land, col=col,breaks=brk,main=t)
    
    #####################################################################
    # Write outputs file for every hour
    rh.out<-getValues(rh_100_land,format="matrix")
    fileout.rh<-paste(dir_rh_out,"rh_",year,"_",month,"_",day,"_",hr,".r",sep="")
    save(rh.out,file=fileout.rh)
    
    #####################################################################
  } # end hr in day loop
} # end day in month loop
    