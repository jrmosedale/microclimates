# Open nscdf files of u/v wind data and store as 3D arrays
# OUtput: wind.v[lon,lat,time] wind.u[lon,lat,time]

library(ncdf4)

infile.u<-"C:/Data/Wind/uwind_uk_2014.nc"
infile.v<-"C:/Data/Wind/vwind_uk_2014.nc"

# Open U vector wind data ncdf4 file
ncdf_windu<-nc_open(infile.u) 
print(ncdf_windu) # file info

# Store dimensional categories as arrays
lat.u<-ncvar_get(ncdf_windu,"lat")
lon.u<-ncvar_get(ncdf_windu,"lon")
time.u<-ncvar_get(ncdf_windu,"time") #  as hours since 01/01/1800 00:00

# Store U vector data as 3d array wind.u[lon,lat,time]
wind.u<-ncvar_get(ncdf_windu)

# *****

# Open V vector wind data ncdf4 file
ncdf_windv<-nc_open(infile.v) 
print(ncdf_windv) # file info

# Store dimensional categories arrays
lat.v<-ncvar_get(ncdf_windv,"lat")
lon.v<-ncvar_get(ncdf_windv,"lon")
time.v<-ncvar_get(ncdf_windv,"time") #  as hours since 01/01/1800 00:00

# Store V vector data as 3d array wind.u[lon,lat,time]
wind.v<-ncvar_get(ncdf_windv)

# *****

# CHECK identical dimensions and ranges 
dim(wind.v); dim(wind.u)
lat.v; lat.u
lon.v;lon.u
min(time.v);min(time.v)
max(time.v); max(time.u)

