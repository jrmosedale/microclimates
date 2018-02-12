# INPUTS to correct:
#       Longitude and Lattitude in sunrise/set Function calls (and check timezone)
#       NOBS in daily weather data in hourly generation loop code (currently set as 333096)


# calcuates Julian Day
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
# calculates sunrise and sunset
suntimes<-function(DOY,year,Lat,Long,Timezone=0,DST=0){
  J<-JD(DOY,year)
  lw<-Long*-1
  n<-J-2451545-0.0009-(lw/360)
  n<-floor(n)+0.5
  sn<-2451545+0.0009+(lw/360)+n
  msa<-(357.5291+0.98560028*(sn-2451545))%%360
  eoc<-1.9148*sin(msa*pi/180)+0.02*sin(2*msa*pi/180)+0.0003*sin(3*msa*pi/180)
  ecl<-(msa+102.9372+eoc+180); ecl<-ecl%%360
  st<-sn+(0.0053*sin(msa*pi/180))-(0.0069*sin(2*ecl*pi/180))
  d<-asin(sin(ecl*pi/180)*sin(23.45*pi/180))
  cos.has<-((sin(-0.83*pi/180)-sin(Lat*pi/180)*sin(d))/(cos(Lat*pi/180)*cos(d)))
  h.set<-vector(length=length(DOY)); h.rise<-h.set; dl<-h.set
  for(i in 1:length(DOY)){
    if(cos.has[i]^2<1){
        has<-acos(cos.has[i])
        J.set<-2451545+0.0009+(((has*180/pi+lw)/360)+n[i]+0.0053*sin(msa[i]*pi/180))-0.0069*sin(2*ecl[i]*pi/180)
        J.rise<-st[i]-(J.set-st[i])
        h.set[i]<-(J.set%%1)*24+Timezone+DST
        h.rise[i]<-(J.rise%%1)*24+Timezone+DST
        dl[i]<-(J.set-J.rise)*24
    }
    if(cos.has[i]>1){
      h.set[i]<-12
      h.rise[i]<-12
      dl[i]<-0
      warning("sun below horizon for 24 hours")
    }
    if(cos.has[i]<(-1)){
      h.set[i]<-0
      h.rise[i]<-0
      dl[i]<-24
      warning("sun above horizon for 24 hours")
      }
  }
  sun.vars<-data.frame(sunrise=h.rise,sunset=h.set,daylight=dl)
  sun.vars
}
# calculates sunrise
sunrise<-function(DOY,year,Lat,Long,Timezone=0,DST=0){
  sun.rise<-suntimes(DOY,year,Lat,Long)[,1]
  sun.rise
}
# calculates sunset
sunset<-function(DOY,year,Lat,Long,Timezone=0,DST=0){
  sun.set<-suntimes(DOY,year,Lat,Long)[,2]
  sun.set
}
# generates hourly values for use between dawn and hotest part of day (~13:35)
generate.hrtemps.day1<-function(min.temp,max.temp,day.length)
{
 x<-c(0:23)
 A=(max.temp-min.temp)/2
 fr=1/(day.length*1.5)
 phase=13.58989
 of=A+min.temp
 y<-A*cos(2*pi*fr*(x-phase))+of
 y
}
# generates hourly values for use between hotest part of day (~13:35) and midnight
generate.hrtemps.day2<-function(min.temp,max.temp,sun.rise)
{
 x<-c(0:23)
 A=(max.temp-min.temp)/2
 lengt<-(24-13.58989)+sun.rise
 fr=1/(lengt*1.5)
 phase=13.58989
 of=A+min.temp
 y<-A*cos(2*pi*fr*(x-phase))+of
 y
}
# generates hourly values between midnight and dawn. Coolest part of day ~12 mins before dawn
generate.hrtemps.night<-function(min.temp,max.temp,sun.rise,sun.set)
{
 x<-c(0:23)
 A=(max.temp-min.temp)/2
 day.length<-sun.set-sun.rise
 lengt<-(24-13.58989)+sun.rise
 fr=1/(lengt*1.5)
 phase=sun.rise-0.2044346
 of=A+min.temp
 y<-A*sin(2*pi*fr*(x-phase)-(2*pi/4))+of
 y
}
# generates hourly values
#inputs:
        # min.temp = minimum daily temperature
        # max.temp = maximum daily temperature
        # next.min = minimum daily temperature the next day
        # prev.max = maximum daily temperature the previous day
        # sun.rise = sunrise time expressed a decimal hour (24hrs) (see functions above)
        # sun.set  = sunset time expressed a decimal hour (24hrs)  (see functions above)
# output:
# vector of 24 values corresponding to estimated temperature in each hour (first value=midnight, last value= 23:00 hrs)
generate.hrtemps<-function(min.temp,max.temp,next.min,prev.max,sun.rise,sun.set)
{
 x<-c(0:23)
 day.length<-sun.set-sun.rise
 day1<-generate.hrtemps.day1(min.temp,max.temp,day.length)
 day2<-generate.hrtemps.day2(next.min,max.temp,sun.rise)
 night<-generate.hrtemps.night(min.temp,prev.max,sun.rise,sun.set)
 y<-day1
 y[15:24]<-day2[15:24]
 y[1:ceiling(sun.rise)]<-night[1:ceiling(sun.rise)]
 y
}



# test of model. Weather is a dataset of weather variables from Culdrose
weather<-read.csv("C:/Data/CuldroseData/allyears.data.csv")
# hourly temperature
temp.hr<-data.frame(temperature=weather$temperature,year=weather$year,
                 day=weather$day,decimal.day=weather$decimal.day,imputed=weather$imputed)
# daily min
temp.mn<-aggregate(temp.hr$temperature,by=list(temp.hr$year,temp.hr$day),min)
temp.mn$year<-temp.mn$Group.1; temp.mn$Group.1<-NULL
temp.mn$day<-temp.mn$Group.2; temp.mn$Group.2<-NULL
temp.mn$mintemp<-temp.mn$x; temp.mn$x<-NULL
# daily max
temp.mx<-aggregate(temp.hr$temperature,by=list(temp.hr$year,temp.hr$day),max)
temp.mx$year<-temp.mx$Group.1; temp.mx$Group.1<-NULL
temp.mx$day<-temp.mx$Group.2; temp.mx$Group.2<-NULL
temp.mx$maxtemp<-temp.mx$x; temp.mx$x<-NULL
temp.both<-merge(temp.mn,temp.mx,by=c("year","day"))
# in correct order
o<-order(temp.both$day)
temp.both<-temp.both[o,]
o<-order(temp.both$year)
temp.both<-temp.both[o,]
# next day temperature minimum
temp.both$next.min<-c(temp.both$mintemp[2:length(temp.both$mintemp)],temp.both$mintemp[length(temp.both$mintemp)])
# previous day maximum temperature
temp.both$prev.max<-c(temp.both$maxtemp[1],temp.both$maxtemp[1:(length(temp.both$maxtemp)-1)])
# calculate sunrise, sunset and day length
temp.both$sunrise<-sunrise(temp.both$day,temp.both$year,Lat=50.0838,Long=-5.25609,Timezone=0,DST=0)
temp.both$sunset<-sunset(temp.both$day,temp.both$year,Lat=50.0838,Long=-5.25609,Timezone=0,DST=0)
temp.both$daylength=temp.both$sunset-temp.both$sunrise
# generate predicted temperatures
pred.temps<-data.frame(year=0,day=0,pred.temp=0)
for (i in 1:length(temp.both$daylength))
{
   tp1<-(length(temp.both$daylength)-i)/100
   tp2<-floor(tp1)
   d<-333096-(dim(pred.temps)[1]-1)
   tp<-paste(tp2,": ",d,sep="")
   if (tp1==tp2) print(tp)
   fit1<-generate.hrtemps(temp.both$mintemp[i],temp.both$maxtemp[i],
                          temp.both$next.min[i],temp.both$prev.max[i],
                          temp.both$sunrise[i],temp.both$sunset[i])
   one.temps<-data.frame(year=temp.both$year[i],day=temp.both$day[i],pred.temp=fit1)
   pred.temps<-rbind(pred.temps,one.temps)
}
pred.temps<-pred.temps[2:length(pred.temps$year),]
pred.temps$year<-as.integer(pred.temps$year)
pred.temps$day<-as.integer(pred.temps$day)
# merge datasets to give hourly plus max and min
temp.hr<-merge(temp.hr,temp.both,by=c("year","day"))
# order temp.hr
o<-order(temp.hr$day)
temp.hr<-temp.hr[o,]
o<-order(temp.hr$year)
temp.hr<-temp.hr[o,]
temp.hr$pred.temp<-pred.temps$pred.temp
temp.hr$hour<-round((temp.hr$decimal.day-temp.hr$day)*24,0)
# days with imputed data   (NB some met station readings missing, which were imputed in original data)
impute.day<-aggregate(temp.hr$imputed,by=list(temp.hr$year,temp.hr$day),sum)
impute.day$year<-impute.day$Group.1; impute.day$Group.1<-NULL
impute.day$day<-impute.day$Group.2; impute.day$Group.2<-NULL
impute.day$dayimp<-impute.day$x; impute.day$x<-NULL
temp.hr<-merge(temp.hr,impute.day,by=c("year","day"))
# remove days with imputed data
sel<-which(temp.hr$dayimp==0)
temp.hr.ni<-temp.hr[sel,]
# some accuracy on 2014 data
sel<-which(temp.hr.ni$year==2014)
hr.2014<-temp.hr.ni[sel,]
o<-order(hr.2014$day)
hr.2014<-hr.2014[o,]
o<-order(hr.2014$year)
hr.2014<-hr.2014[o,]
# plots
par(mfrow=c(2,2),mar=c(5,5,5,5))
# Day 2
sel<-which(hr.2014$day==100)
p1.2014<-hr.2014[sel,]
plot(pred.temp~hour,data=p1.2014,type="l",lwd=2,col="grey",
     xlim=c(0,24),ylim=c(0,14),cex.lab=2,cex.axis=2,
     xlab="",ylab="")
par(new=T)
plot(temperature~hour,data=p1.2014,type="l",
     xlim=c(0,24),ylim=c(0,14),cex.lab=2,cex.axis=2,
     xlab="Hour",ylab=expression(paste("Temperature (",~degree~C,")",sep="")))
# Day 76 to 90
sel<-which(hr.2014$day>75 & hr.2014$day<=90)
p2.2014<-hr.2014[sel,]
plot(pred.temp~decimal.day,data=p2.2014,type="l",lwd=2,col="grey",
     xlim=c(76,90),ylim=c(0,14),cex.lab=2,cex.axis=2,
     xlab="",ylab="")
par(new=T)
plot(temperature~decimal.day,data=p2.2014,type="l",
     xlim=c(76,90),ylim=c(0,14),cex.lab=2,cex.axis=2,
     xlab="Day of year",ylab=expression(paste("Temperature (",~degree~C,")",sep="")))
# Day 176 to 205
sel<-which(hr.2014$day>181 & hr.2014$day<=195)
p3.2014<-hr.2014[sel,]
plot(pred.temp~decimal.day,data=p3.2014,type="l",lwd=2,col="grey",
     xlim=c(181,195),ylim=c(10,24),cex.lab=2,cex.axis=2,
     xlab="",ylab="")
par(new=T)
plot(temperature~decimal.day,data=p3.2014,type="l",
     xlim=c(181,195),ylim=c(10,24),cex.lab=2,cex.axis=2,
     xlab="Day of year",ylab=expression(paste("Temperature (",~degree~C,")",sep="")))

plot(pred.temp~temperature,data=hr.2014,pch=3,cex=0.2,
     xlim=c(-5,30),ylim=c(-5,30),cex.lab=2,cex.axis=2,
     xlab="Observed temperature",ylab="Predicted temperature")
abline(a=0,b=1,lwd=2,col="red",lty=2)









