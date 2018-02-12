# # # # # # # # # # # # # # # #
# Read in CAL file
# # # # # # # # # # # # # # # #
# template for raster
library(raster)
lwr<-function(Temp,RH,CAL)
{
  e0<-0.6108*exp(17.27*Temp/(Temp+237.3)) # saturated vapour pressure
  ea<-e0*(RH/100) # actual vapour pressure
  moisture.absorb<-0.34-0.14*sqrt(ea)
  cloud.absorb<-1.35*(1-CAL)-0.35
  rnl<-2.043*10^-10*(Temp+273.16)^4*moisture.absorb*cloud.absorb
  rnl
}
add.zero<-function(x)
{
 y<-x
 if (y<9) y<-paste("0",x,sep="")
 y
}

year<-1992
month<-6
mth<-add.zero(month)
# template raster
tem.r<-raster("C:/Jonathanmodel/temperature/Temp5km/MaxTemp_1992-06-18_Actual.txt")
# Effective cloud albedo (CAL)
filein<-paste("C:/Jonathanmodel/longwav/CALimp/CALimp_",year,"_",month,".R",sep="")
load(filein)
i<-1
lwr.store<-array(0,dim=c(32,56,72))
for (day in 19:21)
{
  # # # # # # # # # # # # # # #
  # Reads in Relative Humidity
  # # # # # # # # # # # # # # # #
  filein<-paste("C:/Jonathanmodel/relhum/HighRes/RelHum_",year,"_",mth,"-",day,".R",sep="") 
  load(filein)
  for (hr in 1:24)
  {
     #Keep track fo progress
     tp<-72-i
     print(i)
     # rel hum
     m.rh<-relhum.store[,,hr]
     r.rh<-raster(m.rh,xmn=70000,xmx=350000,ymn=0,ymx=160000)
     # CAL
     m.cal<-CAL.imp[,,i]
     r.cal<-raster(m.cal,xmn=-6.725,xmx=-2.675,ymn=49.775,ymx=51.375)
     projection(r.cal)<-"+init=epsg:4326"
     OSGB.cal<-projectRaster(r.cal,crs="+init=epsg:27700")
     OSGB.cal<-resample(OSGB.cal,r.rh)
     e<-extent(c(70000,350000,0,160000))
     OSGB.cal<-crop(OSGB.cal,e)
     # temperature
     filein<-paste("C:/Jonathanmodel/temperature/Temp5km/HRtemp/Temp_",year,"_",mth,"_",add.zero(day),".R",sep="")
     load(filein)
     m.temp<-Hr.Temps[,,hr]
     r.temp<-raster(m.temp,template=tem.r)
     r.temp<-crop(r.temp,e)
     # Calculate longwave radiation
     lwr.store[,,i]<-lwr(getValues(r.temp,format="matrix"),
                         getValues(r.rh,format="matrix"),
                         getValues(OSGB.cal,format="matrix"))
     
     par(mfrow=c(2,2))
     plot(OSGB.cal,main="Effective cloud albedo")
     plot(r.rh,main="Relative humdidity")
     plot(r.temp,main="Temperature")
     plot(raster(lwr.store[,,i],template=r.rh),main="Long-wave radiation")
     i<-i+1
  }
}
save(lwr.store,file="C:/Jonathanmodel/longwav/LWRout/lwwav.R")































