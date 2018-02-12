
test.fun<-function(x,m1,m2) {
  m3<-SZA(x,m1,m2)
  #m3<-(m1*m2) + x
  print(m3)
return(m3)
}

sz.angle<-function(x,m1,m2) {
  
  
d<- 23.45p / 180 * sin(2p * (284 + n) / 365) 
Tsolar<-Tlocal + Eqt / 60 + (Longsm - Longlocal) / 15
#Longsm is the longitude for the standard meridian for the observer's time zone
w = p * (12 - Tsolar) / 12  
  
sza[i] <- 90 - asin(sin(l) * sin(d) + cos(latt) * cos(d) * cos(w))

}

#TEST 
long<-matrix(c(49.5,50,50.5,51),nrow=2,ncol=2)
latt<-matrix(c(0,1,2,3),nrow=2,ncol=2)
long
latt
m3<-(long*latt)+10
x<-datetime
test.fun(x,latt,long)

# Using insol sunvector and sunpos
day<-30;month<-6;year<-2013; hr<-13
lt<-49.5;ln<-0
lt2<-51;ln2<-0

sza.angle<- function(jd,lt,ln,tz=0){
  sv<-sunvector(jd,lt,ln,tz)
  sza<-sunpos(sv)[,2]
  return(sza)
}

sza.angle(jd,lt,ln)


# CUT OUTS FROM HR>TEMPS

# Alternative methods
#newproj<-"+init=epsg:4326"
#r<-projectExtent(day.tmax,newproj)
#grid<-projectRaster(day.tmax,r)  
#llgrid<-projectRaster(day.tmax,crs="+init=epsg:4326")
#long<-xFromCell(grid,c(1:ncell(grid)))
#lat<-yFromCell(grid,c(1:ncell(grid)))


# For TESTING on vectors from p:p+q
#tday<-day.v[p:p+q];tyear<-year.v[p:p+q];tlong<-long[p:p+q];tlat<-lat[p:p+q]
#tsd<-sundown[p:p+q];tsu<-sunup[p:p+q]
#tsu<-sunrise(doy,tyear,tlat,tlong,Timezone=0,DST=0)
#tsd<-sunset(doy,tyear,tlat,tlong,Timezone=0,DST=0)

# create test m of xyz values and matching raster where value = cell number
z<-1:2135
OStest.m<-cbind(os.m,z)
OStest.r<-rasterFromXYZ(test.m, res=c(5000,5000), crs="+init=epsg:27700")  
plot(OStest.r)

LLtest.m<-cbind(ll.m,z)
LLtest.r<-rasterFromXYZ(LLtest.m, crs="+init=epsg:27700") 


# atomic values
lt<-49.5;ln<-0;tz=0
lt2<-51;ln2<-0
jd<-JD(day,month,year)+hr/24
sv1<-sunvector(jd,lt,ln,tz)
sv2<-sunvector(jd,lt2,ln2,tz)

sza1<-sunpos(sv1)[,2]
sza2<-sunpos(sv2)[,2]

#merge sv
sv.a<-as.list(sv1,sv2)
sza.a<-sunpos(sv.a[1])

# argumnets = matrices
lt.m<-matrix(rep(lt,4),nrow=2, ncol=2)
ln.m<-matrix(rep(ln,4),nrow=2, ncol=2)

jd.m<-matrix(rep(jd,4),nrow=2, ncol=2)
tz.m<-matrix(rep(0,4),nrow=2, ncol=2)
sv<-sunvector(jd.m,lt.m,ln.m,tz.m)
sza<-sunpos(sv)[,2]

# replace parts of matrix
your.mat[your.mat == 1] <- 0
matrix(a[cbind(c(row(b)),c(b))],nrow=nrow(a))

night.h<-paste(",1:",as.character(ceiling(sunup)),sep="" )

# REmove all from memory
rm(list=ls())
