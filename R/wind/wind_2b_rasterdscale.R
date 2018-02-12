library(ncdf4)
library(raster)
library(rgdal)
# FUNCTION to calculate raster of wind strength from wind direction and shelter maps\
# Input: raster of wind direction & wind strength,  interval for which shelter maps available
# Internal parameter = dir of shelter maps and names
# Output: final raster of wind strength corrected using appropriate shelter maps
# CHECK - extent of shelter maps compared with dem etc
windstrength<-function(wdir,wstr,interval=5)
{
  dir_shelter<-"~/Documents/Exeter/Data2015/Wind/Shelter/"
  # create raster of wind direction rounded to 5 deg
  wdir<-round((wdir/interval))*interval
  # create result raster all cells set to 0
  final<-raster(wdir) 
  final<-setValues(final,rep(NA,ncell(final)))
  #define stack layers = length(seq(cellStats(rd5,stat="min"),cellStats(rd5,stat="max"),5))
  for(d in seq(cellStats(wdir,stat="min"),cellStats(wdir,stat="max"),5)) { 
    in.file<-paste(dir_shelter,"Shelter_",sprintf("%03d",d,sep=""),"_deg.r",sep="")
    load(in.file) # loads matrix w.coef
    wcoef<-raster(w.coef,template=dem)
    #compare(wc.r,final)
    dir<-wdir %in% d 
    wstr.d<-mask(wstr,dir,maskvalue=0)
    final<-sum(final,(wstr.d*wcoef),na.rm=TRUE)
    #plot(final)
    remove(in.file)
  }
  final<-mask(final,wstr) # converts sea to NA on basis of wstr raster
  return(final)
}

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
# Load wind data (single file for whole time period - created by wind_downscale1)
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v

#dir_wind<-"C:/Data2015/Wind/"
dir_wind<-"~/Documents/Exeter/Data2015/Wind/"
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))

# Specify month and year for which data are required
yr=2014
month=6
day<-10

for (day in 1:31){
  
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
ustr<-calc(u_100_b,fun=function(x){sqrt(x^2)}) # wind strength
vstr<-calc(v_100_b,fun=function(x){sqrt(x^2)}) # wind strength

udir<-calc(u_100_b,fun=function(x){ifelse(x>0,1,-1)})# positive or negative
vdir<-calc(v_100_b,fun=function(x){ifelse(x>0,1,-1)}) # positive or negative

u.adj<-overlay(ustr,dem,fun=function(x,y){return(x*((-0.000000108025)*y^2+0.000408692*y+0.956139))})
v.adj<-overlay(vstr,dem,fun=function(x,y){return(x*((-0.000000108025)*y^2+0.000408692*y+0.956139))})

# adjust values to correspond to wind speed 1m above theground
# rescaling factor first derived by comparing values ot Culdrose wind data using wind_downscale3
# however, in line wiht what you'd expect from: http://www.fao.org/docrep/x0490e/x0490e07.htm#wind profile relationship
u.adj<-u.adj*0.373686439
v.adj<-v.adj*0.373686439
ru<-u.adj*udir
rv<-v.adj*vdir

# Some code here for plotting altitude adjusted values, currently commmented out
#par(mfrow=c(1,1))
#plot(ru,main="altitude adjusted wind u")
#plot(rv,main="altitude adjusted wind v")

# calculate wind direction
direction<-(180/pi)*atan2(ru,rv)  # overlay(ru,rv,fun=function(x,y){(180/pi)*atan2(u.adj,v.adj)}) # NB this is direction in which wind blows to
rd<-calc(direction,fun=function(x){ifelse(x<=180,x+180,x-180)})  # NB this direction from which wind originates (360 deg)

#plot(rd,main="wind direction")


#############
# Stage 4: height adjustments done using shelter coefficient maps based on topography and wind direction
# **** ASSUMES: shelter maps calculated using identical dem ****
#############

# Calculates wind strength from u and v components
rstr<-overlay(ru,rv,fun=function(x,y){return(sqrt(x^2+y^2))})
rstr<-mask(rstr,dem) #set  to NA if sea (on basis of dem)
final<-windstrength(rd,rstr)
#plot(final,main="Final wind stength")

m.str<-getValues(final)
m.dir<-getValues(rd)

fileout.1<-paste(dir_wind,"strength/strength_",yr,"_",month,"_",day,"_",hr,".r",sep="")
fileout.2<-paste(dir_wind,"direction/direction_",yr,"_",month,"_",day,"_",hr,".r",sep="")
save(m.str,file=fileout.1)
save(m.dir,file=fileout.2)
print(fileout.1)
print(fileout.2)

} # end hr
proc.time() - ptm
}#end  day

#
# one day 24 hrs loop
#user   system  elapsed 
#2513.607  662.597 3201.185 
