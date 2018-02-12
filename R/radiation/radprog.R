library(raster)
# Julian day from day of year and year
JD<-function(DOY,year)
{
  month<-ifelse(DOY==365,12,floor(DOY*12/365)+1)
  day<-DOY%%(365/12)
  a<-(14-month)/12
  y<-year+4800-a
  m<-month+12*a-3
  JDN<-day+(153*m+2)/5+365*y+y/4-y/100+y/400-32045
  JDN<-JDN-1.5
  JDN
}
# Converts a raster to a matrix for use with solar index functions
use.raster<-function(r)
{
     xr<-dim(r)[1]
     xc=dim(r)[2]
     m<-array(getValues(r),dim=c(xc,xr))
     m<-t(m)
     m
}
# OS GB grid ref to Lat and Long
OSGBtolatlong<-function(x,y)
{
  pt = data.frame(x,y)
  coordinates(pt)=~x+y
  proj4string(pt)=CRS("+init=epsg:27700")
  latlong<-spTransform(pt,CRS("+init=epsg:4326"))
  ll<-as.data.frame(latlong)
  ll
}
# Needed for solar index function
solartime <- function(localtime,Long,Julian,merid=0,dst=0)
{
	Bn <- 2 * 3.141 * (Julian - 81) / 364
	eot <- 9.87 * sin(2 * Bn) - 7.53 * cos(Bn) - 1.5 * sin(Bn)
	solartime <- localtime + (4 / 60) * (2 * 3.141 * (merid - Long) / 360) + (eot / 60) - dst
	solartime
}
solalt <- function(localtime,Lat,Long,Julian,merid=0,dst=0)
{
  stime<-solartime(localtime,Long,Julian,merid,dst)
	tt <- 0.261799 * (stime - 12)
	declin <- (pi * 23.5 / 180) * cos(2 * pi * ((Julian - 171) / 365.25))
	Sinh = sin(declin) * sin(Lat * pi / 180) + cos(declin) * cos(Lat * 3.141 / 180) * cos(tt)
	solalt = (180 * atan(Sinh / sqrt(1 - Sinh * Sinh))) / pi
	solalt
}
solazi <- function(localtime,Lat,Long,Julian,merid=0,dst=0)
{
  stime<-solartime(localtime,Long,Julian,merid,dst)
	tt = 0.261799 * (stime - 12)
	declin = (pi * 23.5 / 180) * cos(2 * pi * ((Julian - 171) / 365.25))
	Sinh = sin(declin) * sin(Lat * pi / 180) + cos(declin) * cos(Lat * pi / 180) * cos(tt)
	hh = (atan(Sinh / sqrt(1 - Sinh * Sinh)))
	Sinazi = cos(declin) * sin(tt) / cos(hh)
	cosazi = (sin(Lat * pi / 180) * cos(declin) * cos(tt) - cos(pi * Lat / 180) * sin(declin)) / sqrt((cos(declin) *
            sin(tt)) ^ 2 + (sin(pi * Lat / 180) * cos(declin) * cos(tt) - cos(pi * Lat / 180) * sin(declin)) ^ 2)
	solazi = 180 + (180 * atan(Sinazi / sqrt(1 - Sinazi * Sinazi))) / pi
	if (cosazi < 0) {
		if (Sinazi < 0) {
		solazi = 180 - solazi
		} else {
		solazi = 540 - solazi
		}
	}
	solazi
}
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
solarindex <- function(slope,aspect,localtime,Lat,Long,Julian,dtm=array(0,dim=c(1,1)),res=100,merid=0,dst=0,shadow=TRUE)
{
	saltitude<-solalt(localtime,Lat,Long,Julian,merid,dst)
  alt <- saltitude * (pi/180)
	zen <- pi/2 - alt
	sazimuth<-solazi(localtime,Lat,Long,Julian,merid,dst)
	azi <- sazimuth * (pi/180)
	sl <- slope * (pi/180)
	asp <- aspect * (pi/180)
	shadowmask <- array(1,dim(dtm))
	horangle<-horizonangle(dtm,sazimuth)
  if(shadow) {
		shadowmask[horizonangle(dtm,sazimuth)>tan(alt)] <- 0
	}
	index <- array(0,dim(dtm))
	index <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - asp)
	index[index<0] <- 0
	index <- index * shadowmask
index
}
skyview <-function(dtm,steps=36)
{
	sky <- array(1,dim(dtm))
	for (s in 1:steps) {
		sky <- sky-atan(horizonangle(dtm,s*360/steps))/((pi/2)*steps)
	}
sky
}
# # # # # # # # # # # # # # #
# Direct normal radiation: to downscaled direct radiation
# # # # # # # # # # # # # # #
library(ncdf4)
library(raster)
# read in DEM and calculate slope and aspect
dem<-raster("C:/Jonathanmodel/wind/demsw.asc")
# extend dem so divisable by 10 km
e<-extent(c(69000,351000,-1000,161000))
dem<-extend(dem,e) # raster should now comprise
projection(dem)<-"+init=epsg:27700"
slope<-terrain(dem, opt='slope', unit='degrees')
aspect<-terrain(dem, opt='aspect', unit='degrees')
# read in list of files in directory
files<-list.files("C:/Jonathanmodel/radiation/DNI/")
# read in 13th file in directory of the nc files of direct normal irradiance (as an example)
#for (fls in 1:48)
#{
fls<-40
infile.dnr<-paste("C:/Jonathanmodel/radiation/DNI/",files[fls],sep="")
# ensure SIS file read in matches DNR file
infile.sis<-paste("C:/Jonathanmodel/radiation/SIS/SIS",
                  substr(infile.dnr,nchar(infile.dnr)-30,nchar(infile.dnr)),
                  sep="")
ncdf_dnr<-nc_open(infile.dnr)
ncdf_sis<-nc_open(infile.sis)
dnr<-ncvar_get(ncdf_dnr)
sis<-ncvar_get(ncdf_sis)
dnr<-t(dnr)
sis<-t(sis)
# convert to raster
r1<-raster(dnr,xmn=-20,xmx=20.05,ymn=30,ymx=65.05)
r2<-raster(sis,xmn=-20,xmx=20.05,ymn=30,ymx=65.05)
e<-extent(c(-7,-2,49,52))
r1<-crop(r1,e)
r2<-crop(r2,e)
projection(r1)<-"+init=epsg:4326"
projection(r2)<-"+init=epsg:4326"
# reproject to OSGB and set extent to same as DEM
rprj1<-projectRaster(r1,crs="+init=epsg:27700")
rprj2<-projectRaster(r2,crs="+init=epsg:27700")
e<-extent(dem)
rprj1<-crop(rprj1,e)
rprj2<-crop(rprj2,e)
# resample to 100m
r100.1<-resample(rprj1,dem)
r100.2<-resample(rprj2,dem)
par(mfrow=c(2,2))
plot(r100.1)
plot(r100.2)
# work out Julian day and time
jul.base<-JD(1,1983)
hrs<-ncvar_get(ncdf_dnr,"time")
days<-floor(hrs/24)
jul.day<-as.numeric(jul.base+days)
hr<-as.numeric(hrs%%24)
# run each block in turn, ignoring blocks that are entirely sea
blocks<-read.csv("C:/Jonathanmodel/blocks.csv")
sel<-which(is.na(blocks$val)==F)
blocks.good<-blocks[sel,]
# creates raster template for storing direct and diffuse radiation values
direct.raster<-dem*0-9
diffuse.raster<-dem*0-9
total.raster<-dem*0-9
for (block in 1:length(blocks.good$val))
{
   # crop 12 km x 12 km block (includes 1 km buffer for terrain shadow effect)
   vals<-blocks.good[block,]
   e<-extent(c(vals$b.xmn,vals$b.xmx,vals$b.ymn,vals$b.ymx))
   dem.crop<-crop(dem,e)
   slope.crop<-crop(slope,e)
   aspect.crop<-crop(aspect,e)
   # convert values to matrices for use with solar index function
   m.dem<-use.raster(dem.crop)
   m.slope<-use.raster(slope.crop)
   m.aspect<-use.raster(aspect.crop)
   # converts NAs to zeros
   sel<-which(is.na(m.dem)==T); m.dem[sel]<-0
   sel<-which(is.na(m.slope)==T); m.slope[sel]<-0
   sel<-which(is.na(m.aspect)==T); m.aspect[sel]<-0
   # Calculate lat and long of centre of grid
   ll<-OSGBtolatlong(vals$c.xmn+5000,vals$c.ymn+5000)
   lat<-as.numeric(ll[2])
   long<-as.numeric(ll[1])
   si<-solarindex(slope=m.slope,aspect=m.aspect,localtime=hr,
                  Lat=lat,Long=long,Julian=jul.day,dtm=m.dem)
   si.flat<-solarindex(slope=0,aspect=0,localtime=hr,
                  Lat=lat,Long=long,Julian=jul.day,shadow=F)[1,1]
   # Direct normal radiation:
   dnr.crop<-crop(r100.1,e)
   dnr.m<-getValues(dnr.crop,format="matrix")
   # Direct radiation: flat
   dir.flat<-dnr.m*si.flat
   # Diffuse radiation flat
   sis.crop<-crop(r100.2,e)
   sis.m<-getValues(sis.crop,format="matrix")
   dif.m<-sis.m-dir.flat
   # downscaled direct radiation
   direct<-dnr.m*si
   direct<-direct[11:110,11:110] # Crop out 10 km
   # downscaled diffuse radiation
   sv<-skyview(m.dem)
   diffuse<-dif.m*sv
   diffuse<-diffuse[11:110,11:110] # Crop out 10 km
   # add mask back in so that SIs are only produced for land
   mask<-getValues(dem.crop,format="matrix")*0
   mask<-mask[11:110,11:110]
   direct<-direct+mask
   diffuse<-diffuse+mask
   total<-direct+diffuse
   # set new extent
   e<-extent(c(vals$c.xmn,vals$c.xmx,vals$c.ymn,vals$c.ymx))
   dnr.crop<-crop(dnr.crop,e)
   # convert to rasters
   dir.r<-raster(direct,template=dnr.crop)
   dif.r<-raster(diffuse,template=dnr.crop)
   tot.r<-raster(total,template=dnr.crop)
   direct.raster<-mosaic(direct.raster,dir.r,fun=max)
   diffuse.raster<-mosaic(diffuse.raster,dif.r,fun=max)
   total.raster<-mosaic(total.raster,tot.r,fun=max)
   par(mfrow=c(2,2))
   plot(direct.raster)
   plot(diffuse.raster)
   plot(total.raster)
}
fileout1<-paste("C:/Jonathanmodel/radiation/rasters/direct/hr_",hrs,"_all.tif",sep="")
fileout2<-paste("C:/Jonathanmodel/radiation/rasters/diffuse/hr_",hrs,"_all.tif",sep="")
fileout3<-paste("C:/Jonathanmodel/radiation/rasters/total/hr_",hrs,"_all.tif",sep="")
writeRaster(direct.raster,file=fileout1,overwrite=T)
writeRaster(diffuse.raster,file=fileout2,overwrite=T)
writeRaster(total.raster,file=fileout3,overwrite=T)