# NB download SST data from here: http://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
# Download all data from 1963 to 2015
# # # # # # # # # # # # # # # # # # # # # # # #
# This part of the programme unzips the files
# # # # # # # # # # # # # # # # # # # # # # # #
library(R.utils)
gunzip("C:/Jonathanmodel/SST/HadISST1_SST_1961-1990.txt.gz")
gunzip("C:/Jonathanmodel/SST/HadISST1_SST_1991-2003.txt.gz")
for (year in 2004:2015)
{
    gz.name<-paste("C:/Jonathanmodel/SST/HadISST1_SST_",year,".txt.gz",sep="")
    gunzip(gz.name)
}
# # # # # # # # # # # # # # # # # # # # # # # #
# This part of the programme tidies up the data
# and writes out seperate files for each month
# # # # # # # # # # # # # # # # # # # # # # # #
# 1961 to 1990
x <- readLines("C:/Jonathanmodel/SST/HadISST1_SST_1961-1990.txt")
y <- gsub( "-3", " -3", x ) # needed because land data are stored as -32768, but not space delimited
for (year in 1961:1990)
{
 print(year)
 for (month in 1:12)
 {
  # read in file
  eff.month<-(year-1961)*12+month
  start.line<-(eff.month-1)*181+2
  end.line<-start.line+179
  mth<-y[start.line:end.line]
  fileout2<-paste("C:/Jonathanmodel/SST/monthly/HadISST1_SST_",year,"_",month,".txt",sep="")
  writeLines(mth,fileout2)
 }
}
# 1991 to 2003
x <- readLines("C:/Jonathanmodel/SST/HadISST1_SST_1991-2003.txt")
y <- gsub( "-3", " -3", x )
for (year in 1991:2003)
{
 print(year)
 for (month in 1:12)
 {
  # read in file
  eff.month<-(year-1991)*12+month
  start.line<-(eff.month-1)*181+2
  end.line<-start.line+179
  mth<-y[start.line:end.line]
  fileout2<-paste("C:/Jonathanmodel/SST/monthly/HadISST1_SST_",year,"_",month,".txt",sep="")
  writeLines(mth,fileout2)
 }
}
# 2004 to 2014
for (year in 2004:2014)
{
 filein<-paste("C:/Jonathanmodel/SST/HadISST1_SST_",year,".txt",sep="")
 x <- readLines(filein)
 # perform substitution to accoutn for fact that land cells and space delimited
 y <- gsub( "-3", " -3", x )
 for (month in 1:12){
     # read in file
     start.line<-(month-1)*181+2
     end.line<-start.line+179
     mth<-y[start.line:end.line]
     fileout2<-paste("C:/Jonathanmodel/SST/monthly/HadISST1_SST_",year,"_",month,".txt",sep="")
     writeLines(mth,fileout2)
 }
}
# 2015 (data to July)
year<-2015
for (month in 1:7){
     # read in file
     start.line<-(month-1)*181+2
     end.line<-start.line+179
     mth<-y[start.line:end.line]
     fileout2<-paste("C:/Jonathanmodel/SST/monthly/HadISST1_SST_",year,"_",month,".txt",sep="")
     writeLines(mth,fileout2)
 }
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This part of the programme crops to the relevent area
# and resamples data to resolution of 5km to match temperature data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(raster)
sst.crop<-function(month,year)
{
  # read in file and convert to raster
  filein<-paste("C:/Jonathanmodel/SST/monthly/HadISST1_SST_",year,"_",month,".txt",sep="")
  sst.m<-as.matrix(read.table(filein),nrow=180,ncol=360)
  sel<-which(sst.m==-32768)
  sst.m[sel]<-NA
  sst.m<-sst.m/100 # convert values to degrees
  sst.r<-raster(sst.m,xmn=-180,xmx=180,ymn=-90,ymx=90)
  # SW area
  e<-extent(c(-7,-2,49,52))
  sst.r<-crop(sst.r,e)
  # Set land areas to have same temperature as nearby sea
  sst.m<-getValues(sst.r,format="matrix")
  sst.m[1,5]<-(sst.m[1,4]+sst.m[2,5])/2
  sst.r<-raster(sst.m,template=sst.r)
  # Convert projection to OSGB
  projection(sst.r)<-"+init=epsg:4326"
  sst.rp<-projectRaster(sst.r,crs="+init=epsg:27700")
  # resample to 5km resolution
  template.m<-matrix(0,ncol=54,nrow=32)
  template.r<-raster(template.m,xmn=75000,xmx=345000,ymn=0,ymx=160000)
  sst.rs<-resample(sst.rp,template.r)
  tl<-paste("Year: ",year," Month: ",month,sep="")
  plot(sst.rs,main=tl)
  fileout<-paste("C:/Jonathanmodel/SST/monthlytiffs/sst_",year,"_",month,".tif",sep="")
  writeRaster(sst.rs,file=fileout,overwrite=T)
}
for (year in 1982:2014)
 {
 for (month in 1:12)
 {
  sst.crop(month,year)
 }
}
for (month in 1:7) sst.crop(month,2015)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This part of the programme interpolates hourly values
# NB at the moment, the files are output as tiff files
# Therefore for now, avoid setting time period as too long
# or too many files will be output!!!
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
nday.in.month <- function(date)
{
    m<-format(date, format="%m")
    while(format(date, format="%m") == m) date <- date + 1
    return(as.integer(format(date - 1, format="%d")))
}

sst.time.int<-function(year,month,day,hr)
{
 # How many days in month
 date<-as.Date(paste(year,"-",month,"-",day,sep=""),"%Y-%m-%d")
 dimth<-nday.in.month(date)
 # Read in SSSTs for the two months that lie either side of the date in question
 mth1<-month
 if (day+hr/24<dimth/2) mth1<-mth1-1  # reads in month before if date is in first half of month
 mth2<-mth1+1
 filein1<-paste("C:/Jonathanmodel/SST/monthlytiffs/sst_",year,"_",mth1,".tif",sep="")
 filein2<-paste("C:/Jonathanmodel/SST/monthlytiffs/sst_",year,"_",mth2,".tif",sep="")
 r1<-raster(filein1)
 r2<-raster(filein2)
 v1<-getValues(r1,format="matrix")
 v2<-getValues(r2,format="matrix")
 # work out weighting to attach
  # How many days in both months
 date1<-as.Date(paste(year,"-",mth1,"-",day,sep=""),"%Y-%m-%d")
 date2<-as.Date(paste(year,"-",mth2,"-",day,sep=""),"%Y-%m-%d")
 dimth1<-nday.in.month(date1)
 dimth2<-nday.in.month(date2)
  # time after dimth1
 tt<-day+hr/24+dimth1/2
 if (day+hr/24>dimth/2) tt<-day+hr/24-dimth1/2
 wgt<-tt/((dimth1+dimth2)/2)
 v<-v2*wgt+v1*(1-wgt)
 r.hr<-raster(v,template=r1)
 tl<-paste("Year: ",year," Month: ",month," Day: ",day," Hour: ",hr,sep="")
 plot(r.hr,main=tl)
 fileout<-paste("C:/Jonathanmodel/SST/hourlytiffs/sst_",year,"_",month,"_",day,"_",hr,".tif",sep="")
 writeRaster(r.hr,file=fileout,overwrite=T)
}
for (day in 1:30){
for (hr in 0:23){
  sst.time.int(2014,6,day,hr)
}}




