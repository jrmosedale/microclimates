
# Create arrays of two years including a leap year for testing JD functions etc

days<-c(1:31,1:29,1:31,1:30,1:31,1:30,1:31,1:31,1:30,1:31,1:30,1:31,1:31,1:28,1:31,1:30,1:31,1:30,1:31,1:31,1:30,1:31,1:30,1:31)

months<-c(rep(1,31),rep(2,29),rep(3,31),rep(4,30),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31),rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,31),rep(11,30),rep(12,31))

years<-c(rep(1980,366),rep(1981,365))

length(days)==length(months)
length(months)==length(years)

# Functions to be tested
JDdoy<-function(DOY,year)
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
JD<-function(day,month,year){
  a<-(14-month)/12
  y<-year+4800-a
  m<-month+12*a-3
  JDN<-as.integer(floor((153*m+2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045 + day)
  JDN
}

JDres<-JD(days,months,years)
JDdoyres<-c(JDdoy(1:366,rep(1980,366)),JDdoy(1:365,rep(1981,365)))

# From http://stackoverflow.com/questions/27757994/julian-dates-in-r-chron-versus-us-naval-observatory
a <- floor((14 - months) / 12)
y <- years + 4800 - a
m <- months + 12 * a - 3
julian2 <- days + floor((153*m + 2)/5) + 365*y + floor(y/4) - 32083


# Test chron package - doesn't work for negative origins

library("chron")
options(chron.origin = c(month=1, day=1, year= ))

julres <- julian(months,days,years)



# Corrected functions
# Functions to be tested
JDdoy<-function(DOY,year)
{
  month<-ifelse(DOY>=365,12,floor(DOY*12/365)+1)
  day<-DOY%%(365/12)
  a<-(14-month)/12
  y<-year+4800-a
  m<-month+12*a-3
  JDN<-day+(153*m+2)/5+365*y+y/4-y/100+y/400-32045
  JDN<-JDN-1.5
  JDN
}
