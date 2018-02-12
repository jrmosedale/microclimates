library(ncdf4)
library(raster)
library(rgdal)
# Function for computing Julian data
# In this instance just used to work out number of days after 1st Jan 1950,
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
# Based on hour, day, month and year, extracts the required value form the array of values stored by wind_downscale1.R
# Inputs:
  # hr: the hour (0-23) 0 = midnight
  # day: the day of the month (0-31)
  # month: the month of the year (numeric: 0-12)
  # year: any year from 1950 to 2014
# Output: the element of the array stored by wind_downscale1.R that corresponds to either that hour, or the latest
# period immediatly before that hour (data only available 6-hourly)
array.val<-function(hr,day,month,yr)
{
 jd.base=JD(1,1,1950)
 jd<-JD(day,month,yr)
 dval<-(jd-jd.base)*4
 hval<-floor(hr/6)
 val<-dval+hval+1
 val
}
# Works out the angle to the horizon in a specified direction (used to calculate the shelter coefficient)
# Inputs:
  # dtm = a digital eleveation model stored as a matrix

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
####################################################
# # # This bit downscales the wind
####################################################
# NB running this for 1 month takes ~30 hours
# # # Stages:
# (1) get wind values for a given hour, day, month and year
# (2) convert to 100m resolution raster OSGB grid reference
# (3) adjust based on altitude
# (4) adjust based on shelter coefficient
# loads data output by wind_downscale1.R
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v
load(file="C:/Jonathanmodel/wind/newdata/uwind.r")
load(file="C:/Jonathanmodel/wind/newdata/vwind.r")
# Specify month and year for which data are required
yr=2014
month=1
# set period for which you want to create 100m resolution wind data
# Set to do all hours in January 2014
for (day in 1:31){
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
# Convert to raster (original lat long format and resolution
uwind.r<-raster(uwind,xmn=-7.5,xmx=-2.5,ymn=47.5,ymx=52.5)
vwind.r<-raster(vwind,xmn=-7.5,xmx=-2.5,ymn=47.5,ymx=52.5)
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
# read in Digital Eelevation data - I've chopped this to just cover SW Britain
# We probably need to make a definative decision as to what area we want to cover and adjust
# areas accordingly.
dem<-raster("C:/Jonathanmodel/wind/demsw.asc")
e<-extent(dem)
# Crop areas so that it neatly has whole 100m grid cells around the edge
# Dem and wind rasters set to same extent here as well
xmn=round(e@xmin/100)*100
xmx=round(e@xmax/100)*100
ymn=round(e@ymin/100)*100
ymx=round(e@ymax/100)*100
e<-extent(c(xmn,xmx,ymn,ymx))
u_100<-crop(u_100,e)
v_100<-crop(v_100,e)
dem<-crop(dem,e)
#############
# Stage 3: adjust based on altitude of
# NB Height adjustment based on wind spped values at different pressures downloaded  Earth System Research Lab
# Typical heights at different pressures calculated from Allen et al 1998 http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric pressure (p)
# Quadratic function fitted - NB this works well for heights up to ~1800m. IT won't work above ~2000m
# Function was first derived by comparing values at different pressures (heights) over the course of a year (2014)
#############
# adjust wind speeds by height of dem
 # convert to matrices
uwind.m<-getValues(u_100,format="matrix")
vwind.m<-getValues(v_100,format="matrix")
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
ru<-raster(u.adj,template=u_100)
rv<-raster(v.adj,template=v_100)
# Some code here for plotting altitude adjusted values, currently commmented out
#par(mfrow=c(2,2))
#plot(ru,main="altitude adjusted wind u")
#plot(rv,main="altitude adjusted wind v")
# calculate wind direction
direction = (180/pi)*(atan2(u.adj,v.adj))  # NB this is direction in which wind blows to
direction<-ifelse(direction<=180,direction+180,direction-180) # NB this direction from which wind originates
rd<-raster(direction,template=u_100)
#plot(rd,main="wind direction")
#############
# Stage 4: height adjustments done using a shelter coefficient based on topography and wind direction
#############
# go through in 10 km blocks and adjust by shelter coefficient
# note, however that actually it actually selects 30km x 30 km area to allow for sheltering effects
# that operate outside the area of each block. The 10km x 10km centre of the block is then selected
# Programme could probably be speeded up without major loss of accuracy by setting buffer to ~5km instead of 10km
# NB chopping into 10km blocks is necessary, as the function for calculating
# the shelter coefficient assumes a single wind direction, a fairly safe assumption over 10km, but not over entire study region
# first extend rasters to match 10km cells
e<-extent(c(60000,360000,-10000,170000))
ru<-extend(ru,e)
rv<-extend(rv,e)
dem<-extend(dem,e)
rd<-extend(rd,e)
m.u<-getValues(ru,format="matrix")
m.v<-getValues(rv,format="matrix")
m.dem<-getValues(dem,format="matrix")
m.d<-getValues(rd,format="matrix")
# creates matrix for storing values
windstrength<-matrix(NA,nrow=1800,ncol=3000) # matrix for storing all values
# Goes through and does each 10km block
for (rws in 1:16)
{
 for (cls in 1:28)
 {
  xmn<-rws*100+1-100
  ymn<-cls*100+1-100
  xmx=xmn+300-1
  ymx=ymn+300-1
  b.u<-m.u[xmn:xmx,ymn:ymx]
  b.v<-m.v[xmn:xmx,ymn:ymx]
  b.dem<-m.dem[xmn:xmx,ymn:ymx]
  sel<-which(is.na(b.dem)==T)
  b.dem[sel]<-0
  # calculates mean direction
  b.dir<-mean(m.d[xmn:xmx,ymn:ymx],na.rm=T)
  # Calculates wind strength from u and v components
  m.str<-sqrt(b.u^2+b.v^2)
  wcoef<-matrix(0,nrow=300,ncol=300)
  # Applies shelter coefficient if wind direction not NA.
  # Wind direction would be NA if all values within 10km block are NA, which happens if the entire 10km block is sea
  # (for which values nOt required)
  if (is.na(b.dir)==F) wcoef<-windindex(b.dem,b.dir)
  m.str<-m.str*wcoef
  # selects data for just the 10km x 10km centre of each 30km x 30km block
  windstrength[(xmn+100):(xmx-100),(ymn+100):(ymx-100)]<-m.str[101:200,101:200]
 }
}
# converts to raster and crops to desired extent
r<-raster(windstrength,template=dem)
e<-extent(c(79400,343500,0,159300))
r<-crop(r,e)
rd<-crop(rd,e)
#plot(r,main="wind speed")
# converts raster to matrix and saves matrix as R dataset
# seperate datasets saved for each hour (entire study area)
m1.out<-getValues(r,format="matrix")
m2.out<-getValues(rd,format="matrix")
fileout.1<-paste("C:/Jonathanmodel/wind/dataout/strength_",yr,"_",month,"_",day,"_",hr,".r",sep="")
fileout.2<-paste("C:/Jonathanmodel/wind/dataout/direction_",yr,"_",month,"_",day,"_",hr,".r",sep="")
save(m1.out,file=fileout.1)
save(m2.out,file=fileout.2)
}}









