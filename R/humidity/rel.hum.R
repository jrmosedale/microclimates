library(raster)
library(insol) # require package to be installed
library(raster)
library(ncdf4)
add.zero<-function(x)
{
 y<-x
 if (y<9) y<-paste("0",x,sep="")
 y
}
rel.to.abs<-function(rh,t)
{
    e0.1<-0.6108*exp(17.27*t/(t+237.3))
    e<-e0.1*(rh/100)
    abso<-(2165*e)/(t+273.16) # grams per metre cubed
    abso
}
abs.to.rel<-function(ab,t)
{
 e<-ab*(t+273.16)/2165
 e0.1<-0.6108*exp(17.27*t/(t+237.3))
 rh<-100*(e/e0.1)
 rh
}
reslice.raster<-function(r)
{
 arr<-array(0,dim=c(73,144))
 e1<-extent(1.25,181.25,-91.25,91.25)
 e2<-extent(181.25,358.75,-91.25,91.25)
 e3<-extent(-1.25,1.25,-91.25,91.25)
 arr[,1:71]<-getValues(crop(r,e2),format="matrix")
 arr[,72]<-getValues(crop(r,e3),format="matrix")
 arr[,73:144]<-getValues(crop(r,e1),format="matrix")
 r2<-raster(arr,xmn=-178.25,181.25,-91.25,91.25)
 r2
}
reproject.r<-function(a,e)
{
 r<-raster(a,xmn=-8.75,xmx=1.25,ymn=46.25,ymx=53.75)
 projection(r)<-"+init=epsg:4326"
 r2<-projectRaster(r,crs="+init=epsg:27700")
 r3<-resample(r2,tr1)
 r4<-crop(r3,e)
 r4
}


# # # # # # # # # # # # # # # #
# Reads in Relative Humidity
# # # # # # # # # # # # # # # #
rh.brick<-brick("C:/Jonathanmodel/relhum/RH/rhum.sig995.1992.nc") # all times
# get times
netRH<-nc_open("C:/Jonathanmodel/relhum/RH/rhum.sig995.1992.nc")
tm<-ncvar_get(netRH,"time")


# template raster
r<-raster("C:/Jonathanmodel/temperature/Temp5km/MaxTemp_1992-06-18_Actual.txt")
year<-1992
month<-add.zero(6)
start.day<-19
end.day<-23
# Store all temperature values in one array
hrtemp.store<-array(0,dim=c(290,180,24*5))
i<-1
for (day in start.day:end.day)
{
   dy<-add.zero(day)
   filein1<-paste("C:/Jonathanmodel/temperature/Temp5km/HRtemp/Temp_",year,"_",month,"_",dy,".R",sep="")
   load(filein1)
   hrtemp.store[,,i:(i+23)]<-Hr.Temps
   i<-i+24
}
i<-1
max.i<-(end.day-start.day+1)*24
for (day in start.day:end.day)
{
   relhum.store<-array(0,dim=c(32,56,24))
   for (period in 0:3)
   {
      hr<-period*6
      # read in temperature raster1, crop and resample
      if (i==0) tr1<-raster(hrtemp.store[,,i+1],template=r)
      if (i>0)  tr1<-raster(hrtemp.store[,,i],template=r)
      # read in temperature raster1, crop and resample
      if (i==max.i-5) tr6<-raster(hrtemp.store[,,i+5],template=r)
      if (i<max.i-5)  tr6<-raster(hrtemp.store[,,i+6],template=r)
      # reproject and crop
      projection(tr1)<-"+init=epsg:27700"
      projection(tr6)<-"+init=epsg:27700"
      tr1.ll<-projectRaster(tr1,crs="+init=epsg:4326")
      tr6.ll<-projectRaster(tr6,crs="+init=epsg:4326")
      e<-extent(c(-8.75,1.25,46.25,53.75))
      tr1.ll<-crop(tr1.ll,e)
      tr6.ll<-crop(tr6.ll,e)
      tem<-raster(array(0,dim=c(3,4)),xmn=-8.75,xmx=1.25,ymn=46.25,ymx=53.75)
      t1<-resample(tr1.ll,tem)
      t6<-resample(tr6.ll,tem)
      # read in correct rh grids and crop
      Jul.base<-JD(ISOdate(1800,1,1))
      Jul.actual<-JD(ISOdate(year,month,day))
      hr.val<-(Jul.actual-Jul.base)*24+hr
      sel<-which(tm==hr.val)
      rh.all1<-subset(rh.brick,sel)
      rh.all2<-subset(rh.brick,(sel+1))
      rh.all1<-reslice.raster(rh.all1)
      rh.all2<-reslice.raster(rh.all2)
      rh1<-crop(rh.all1,e)
      rh2<-crop(rh.all2,e)
      # Convert to absolute humidity
      abs.hum0<-rel.to.abs(getValues(rh1,format="matrix"),getValues(t1,format="matrix"))
      abs.hum6<-rel.to.abs(getValues(rh2,format="matrix"),getValues(t6,format="matrix"))
      # interpolate for missing hours
      abs.hum1<-(abs.hum0*5+abs.hum6*1)/6
      abs.hum2<-(abs.hum0*4+abs.hum6*2)/6
      abs.hum3<-(abs.hum0*3+abs.hum6*3)/6
      abs.hum4<-(abs.hum0*2+abs.hum6*4)/6
      abs.hum5<-(abs.hum0*1+abs.hum6*5)/6
      # Convert to high res absolute humidity
      e<-extent(c(70000,350000,0,160000))
      abs.h1<-reproject.r(abs.hum1,e)
      abs.h2<-reproject.r(abs.hum2,e)
      abs.h3<-reproject.r(abs.hum3,e)
      abs.h4<-reproject.r(abs.hum4,e)
      abs.h5<-reproject.r(abs.hum5,e)
      abs.h6<-reproject.r(abs.hum6,e)
      # crop temperatures to same extent
      tr2<-raster(hrtemp.store[,,i+2],template=r)
      tr3<-raster(hrtemp.store[,,i+3],template=r)
      tr4<-raster(hrtemp.store[,,i+4],template=r)
      tr5<-raster(hrtemp.store[,,i+5],template=r)
      tr1<-crop(tr1,e)
      tr2<-crop(tr2,e)
      tr3<-crop(tr3,e)
      tr4<-crop(tr4,e)
      tr5<-crop(tr5,e)
      tr6<-crop(tr6,e)
      # Convert to high res relative humidity
      j<-period*6
      relhum.store[,,j+1]<-abs.to.rel(getValues(abs.h1,format="matrix"),getValues(tr1,format="matrix"))
      relhum.store[,,j+2]<-abs.to.rel(getValues(abs.h2,format="matrix"),getValues(tr2,format="matrix"))
      relhum.store[,,j+3]<-abs.to.rel(getValues(abs.h3,format="matrix"),getValues(tr3,format="matrix"))
      relhum.store[,,j+4]<-abs.to.rel(getValues(abs.h4,format="matrix"),getValues(tr4,format="matrix"))
      relhum.store[,,j+5]<-abs.to.rel(getValues(abs.h5,format="matrix"),getValues(tr5,format="matrix"))
      relhum.store[,,j+6]<-abs.to.rel(getValues(abs.h6,format="matrix"),getValues(tr6,format="matrix"))
      plot(raster(relhum.store[,,j+1],template=tr1))
   }
   fileout<-paste("C:/Jonathanmodel/relhum/HighRes/RelHum_",year,"_",month,"-",day,".R",sep="")
   save(relhum.store,file=fileout)
}

