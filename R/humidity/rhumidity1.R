library(ncdf4)
library(raster)
library(rgdal)

# Import ncdf relative humidity files (organised by year)
dir_rh<-"C:/Data2015/RelHumidity/"

# Import from single file of 1980-2015 data
#infile.rh<-paste(dir_rh,"X144.173.78.19.237.7.10.4.nc",sep="")
infile.rh<-paste(dir_rh,"rhum.sig995.1992.nc",sep="")
print(infile.rh)
ncdf_rh<-nc_open(infile.rh) # summary of file variables,  dimensions, attributes
rh<-ncvar_get(ncdf_rh) # read data
d<-dim(rh) # dimension lengths of time
save(rh,file=paste(dir_rh,"rh.r",sep=""))



# ALTERNATIVE IF MULTIPLE FILES 
# INPUT ncdf files and create single array rh of all years of data
yearsin<-1980:2015
filenum<-length(yearsin)

for (i in 1:filenum){
  infile.rh<-paste(dir_rh,"rhum.sig995.",yearsin[i],".nc",sep="")
  ncdf_rh<-nc_open(infile.rh) # summary of file variables,  dimensions, attributes
  rh<-ncvar_get(ncdf_rh) # read data
  d<-dim(rh) # dimension lengths of time
  store[i]<-d[3] 
  
  # is it worth loading all files/years into single array???
  } # filenum

# Calculates and defines size od array for all data 
rh<-array(0,dim=c(d[1],d[2],sum(store))) # defines array size using sum of time and d1 and d2 of last file

for (i in 1:filenum){
  infile.rh<-paste(dir_rh,"rhum.sig995.",yearsin[i],".nc",sep="")
  ncdf_rh<-nc_open(infile.rh) 
  rh<-ncvar_get(ncdf_rh) 
  mx<-mn+store[i]-1
  print(paste("min=",mn," max=",mx,sep=""))
  rh[,,mn:mx]<-rh
  mn<-mn+store[i]
}
save(rh,file=paste(dir_rh,"rh.r",sep=""))

#sis.r<-raster(sis,xmn=-20,xmx=20,ymn=30,ymx=65,res )






