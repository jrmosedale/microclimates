# Julian date functions for use
library(insol)

# Calculates julian date for 12:00 on the day of parameters
# Used in: setup, t5km_to_hourly_blocks
JDdmy<-function(day,month,year) # correct
{
  options(digits=12) 
  jdate<-insol::JD(ISOdate(year,month,day))
  return(jdate)
}

# Function to convert JD to day, month year as list - allows for jd fractions
# Gives Gregorian dates (even for julian period before 15C) 
# Inputs are vectors; Output as Data.frame
# Check with: http://aa.usno.navy.mil/faq/docs/JD_Formula.php
# Used in: setup, t5km_to_hourly_blocks
DMYjd<-function(jd) # correct 
{
  options(digits=12) 
  newdate<-as.POSIXlt(insol::JD(jd,inverse=TRUE))
  dmy<-data.frame(day=newdate$mday,month=newdate$mon+1,year=newdate$year+1900)
  return (dmy)
}

JDdoy<-function(doy,year)
{
  options(digits=12) 
  newdate<-insol::doyday(year,doy)
  jdate<-insol::JD(ISOdate(newdate$year+1900,newdate$mon+1,newdate$mday))
  return(jdate)
}


# Calculates days in month (from 1980-2015) from day, month, year
# Used to write monthly files of imputed values 
days.in.month<-function(day,month,year){
  y<-year-1979 # so 1=1980
  feb.d<-c(29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28,29,28,28,28,
           29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28) # days of Feb from 1980 to 2015
  monthdays<-c(31,feb.d[y],31,30,31,30,31,31,30,31,30,31)
  return(monthdays[month])
}


############################################################################
# TESTS
############################################################################
JDdmy(1,1,1980)
JDdmy(2,1,1980)
JDdoyy(1,1980)
JDdoyy(2,1980)
test<-DMYjd(2444240)
DMYjd(2444241.4)
test$day

days<-c(1:31,1:29,1:31,1:30,1:31,1:30,1:31,1:31,1:30,1:31,1:30,1:31,1:31,1:28,1:31,1:30,1:31,1:30,1:31,1:31,1:30,1:31,1:30,1:31)
months<-c(rep(1,31),rep(2,29),rep(3,31),rep(4,30),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31),rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))
years<-c(rep(1980,366),rep(1981,365))

res<-JDdmy(days,months,years)
DMYjd(res)
