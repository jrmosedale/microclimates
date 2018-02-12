# INPUTS to correct:
#       Longitude and Lattitude in sunrise/set Function calls (and check timezone)
#       NOBS in daily weather data in hourly generation loop code (currently set as 333096)

##########################################################################################
### FUNCTIONS
##########################################################################################

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

##########################################################################################

year<-2011
month<-6
day<-10

# FOR year/month Loops 

For (hr in 0:23) {}
#  read max and min day files and create max/min rasters - READ prev,day,next
# loop with prev<-day, day<-next, READ next
dir_temp<-"C:/Data2015/Temp5km/"
dir_temp<-"~/Documents/Exeter/Data2015/Temp5km/extract/"

max.infile<-paste(dir_temp,"MaxTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_ACTUAL.txt", sep="")
min.infile<-paste(dir_temp,"MinTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_ACTUAL.txt", sep="")
daily.tmax<-raster(max.infile, layer=1,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")
daily.tmin<-raster(min.infile, layer=2,xmn=-200000, ymn=-200000,nrows=180,ncols=290, res=c(5000,5000), crs="+init=epsg:27700")

daily.stk<-stack(daily.tmax,daily.tmin)
names(daily.stk)<-c("tmax","tmin")
par(mfrow=c(1,1)) # for printing
plot(daily.stk)

# Crop to extent of dem
daily.stk<-crop(x=daily.stk,y=dem)




# Calculate next day min and prev day max

# Calculate hourly temperatures - 1. calculate sunrise, sunset and day length
temp.both$sunrise<-sunrise(temp.both$day,temp.both$year,Lat=50.0838,Long=-5.25609,Timezone=0,DST=0)
temp.both$sunset<-sunset(temp.both$day,temp.both$year,Lat=50.0838,Long=-5.25609,Timezone=0,DST=0)
temp.both$daylength=temp.both$sunset-temp.both$sunrise

# 2. generate predicted temperatures
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


# Write hourly temp files ? Orig resolution as downscaling will depend upon other variables??



} # end For hr in day loop

# resample to 100m - do this before or after creation of hourly temps? 
# include elevation correction based on 100m cell difference from 5km mean
# ?? coastal effects etc??
