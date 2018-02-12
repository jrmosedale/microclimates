library(raster)
library(rgdal)
# # # # # # # # # # # # # # # # # # # # # #
# # # Compare model with Culdrose wind data
# # # # # # # # # # # # # # # # # # # # # #
# read in two files:
#(1) Culdrose weather data - some of these data are  missing and the impution is a bit dodgy
#(2) tic - this is a version that has better impututed data, (but the time series isn't quite as long)
weather<-read.csv("C:/Lizardrefugiatidy/Data/Satelittemodeling/allyears.data.csv")
tic<-read.csv("C:/Lizardrefugiatidy/Data/Satelittemodeling/Culdroseall_tic_2014pasted.csv")
# Uses tic imputed values:
weather$pressure[169:332928]<-tic$pressure
weather$temperature[169:332928]<-tic$temperature
weather$rel.hum[169:332928]<-tic$rel.hum
weather$windspeed[169:332928]<-tic$windspeed
weather$winddir[169:332928]<-tic$winddir
weather$cloudcover[169:332928]<-tic$cloudcover
weather$winddir<-ifelse(weather$winddir<355,weather$winddir,0)
# extract data for required period only
sel<-which(weather$year==2014)
weather<-weather[sel,]
sel<-which(weather$decimal.day<32)
weather<-weather[sel,]
# extract wind data
wind.data<-data.frame(c.windspeed=weather$windspeed*0.51444444,c.winddir=weather$winddir)
# converts to wind speed at 1 metre height (Culdrose anometer @ 32m)
# Method based on Allen et al. 1998: http://www.fao.org/docrep/x0490e/x0490e07.htm#wind profile relationship
wind.data$c.windspeed<-wind.data$c.windspeed*0.6341314*0.8487155
# Convert wind to u and v components
# Important: u and v components measure direction wind is blowing too
            #wind directions are usually measured in direction wind is coming from
            #hence reason that wind direction has 180 added to it
wind.data$c.uwind<-wind.data$c.windspeed*sin((wind.data$c.winddir+180)*(pi/180))
wind.data$c.vwind<-wind.data$c.windspeed*cos((wind.data$c.winddir+180)*(pi/180))
# # # # # # # # # # # # # # # # # # # # # # # # #
# Obtain modelled wind data for each day and hour
# In this case, all data for January 2014
# # # # # # # # # # # # # # # # # # # # # # # # #
# create variables in data frame
wind.data$p.windspeed<-0
wind.data$p.winddir<-0
i<-1
for (day in 1:31)
{
 # looping can take quite a long time. This keeps tabs on progress printing the day
 tp<-paste("Day: ",day,sep="")
 print(tp)
 for (hr in 0:23)
 {
  # loads arrays of wind speed and direction for each hour produce by wind_downscale2.R
  # NB - as original names of these were m1.out and m2.out, this is the name they are automatically assigned here
  filein1<-paste("C:/Jonathanmodel/wind/dataout/strength_2014_1_",day,"_",hr,".r",sep="")
  filein2<-paste("C:/Jonathanmodel/wind/dataout/direction_2014_1_",day,"_",hr,".r",sep="")
  load(filein1)
  load(filein2)
  # Convert arrays to rasters
  r.speed<-raster(m1.out,xmn=79400,xmx=343500,ymn=0,ymx=159300)
  r.direction<-raster(m2.out,xmn=79400,xmx=343500,ymn=0,ymx=159300)
  # Extract data for Culdrose (Easting = 167162.8, Northing =  25489.81)
  xy<-data.frame(x=167162.8,y=25489.81)
  # store values in data frame
  wind.data$p.windspeed[i]<-extract(r.speed,xy)
  wind.data$p.winddir[i]<-extract(r.direction,xy)
  i<-i+1
 }
}
# Convert wind to u and v components
# Again, u and v components measure direction wind is blowing too
# but wind directions are usually measured in direction wind is coming from
wind.data$p.uwind<-wind.data$p.windspeed*sin((wind.data$p.winddir+180)*(pi/180))
wind.data$p.vwind<-wind.data$p.windspeed*cos((wind.data$p.winddir+180)*(pi/180))
# Comparison plots + lines through points (abline function)
par(mfrow=c(2,2))
plot(c.windspeed~p.windspeed,data=wind.data,xlab="Modelled wind speed",ylab="Culdrose wind speed")
abline(lm(c.windspeed~p.windspeed,data=wind.data))
# Wind direction a little misleading as not on circular plot, so 359 looks vastly different from 0
plot(c.winddir~p.winddir,data=wind.data,xlab="Modelled wind direction",ylab="Culdrose wind direction")
abline(lm(c.winddir~p.winddir,data=wind.data))
plot(c.uwind~p.uwind,data=wind.data,xlab="Modelled easterly wind component",ylab="Culdrose easterly wind component")
abline(lm(c.uwind~p.uwind,data=wind.data))
plot(c.vwind~p.vwind,data=wind.data,xlab="Modelled northerly wind component",ylab="Culdrose northerly wind component")
abline(lm(c.vwind~p.vwind,data=wind.data))
# Linear regression of modelled versus Culdrose data, with intercept forced at zero
summary(lm(c.windspeed~p.windspeed+0,data=wind.data))




