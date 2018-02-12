# See: http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm 

library(lib.loc = .Library)
library()
library(ncdf4)
library(chron)

windu<-nc_open("C:/Data/Wind/uwind_uk_2014.nc")
print(windu) # file info
windu$ndim # number of dimensions
windu$nvars # number of vars

windu$dim[1:3] # details of dimensions 1 to 3

windu.data<-ncvar_get(windu) # get data


windv<-nc_open("C:/Data/Wind/vwind_uk_2014.nc")
print(windv) # file info
windv$nvars # number of vars
windv.data<-ncvar_get(windv)

# get selected data by dimensions lat,lon,time
# calc time vals for 01-1965 to 01-1970
# calculate as hours since 01/01/1800 00:00


windv.1965<-ncvar_get(windv)

# Read lat, lon and time into vectors
lat<-ncvar_get(windu,"lat")
lon<-ncvar_get(windu,"lon")
time<-ncvar_get(windu,"time")
tunits<-ncatt_get(windu,"time","units") # gets units definition
u<-ncvar_get(windu,"uwnd")




# Convert time value - CORRECT
tustr <- strsplit(tunits$value, " ")
thstr <- strsplit(unlist(tustr)[4], ":")
tdstr <- strsplit(unlist(tustr)[3], "-")
thour = as.integer(unlist(thstr[1]))
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
tnew<-chron(time, origin = c(thour,tmonth, tday, tyear))


# *************************************

# Read in dimension values

windu.all<-nc_open("C:/Data/Wind/X144-1960-70-uwindx4.nc",readunlim=TRUE)

test<-ncvar_get(windu)

nc_close(windu)
nc_close(windv)