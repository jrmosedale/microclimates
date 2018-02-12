library(ncdf4)

#### Set directory with wind vector  data files (u and v) ###
### also used to store R wind arrays ###
dir_wind<-"C:/Data2015/Wind/"

######################################################################################
# # This bit of the programme converts netcdf files to a single array and stores array
######################################################################################
# define file time groupings
vars<-c(1960,1970,1980,1990,2000,2010)
filenum<-6

# u wind
# This bit works out what size the array should be
store<-0
for (i in 1:filenum){
    infile.u<-paste(dir_wind,"u_",vars[i],".nc",sep="")
    ncdf_windu<-nc_open(infile.u) # summary of file variables,  dimensions, attributes
    wind.u<-ncvar_get(ncdf_windu) # read data
    d<-dim(wind.u) # dimension lengths of wind.u
    store[i]<-d[3] 
}
wind_u<-array(0,dim=c(d[1],d[2],sum(store))) # defines array size using sum of time and d1 and d2 of last file
# This bit stores the data in an array
mn<-1
for (i in 1:filenum)
{
 infile.u<-paste(dir_wind,"u_",vars[i],".nc",sep="")
 ncdf_windu<-nc_open(infile.u)
 wind.u<-ncvar_get(ncdf_windu)
 mx<-mn+store[i]-1
 tp<-paste("min=",mn," max=",mx,sep="")
 print(tp)
 wind_u[,,mn:mx]<-wind.u
 mn<-mn+store[i]
}
save(wind_u,file=paste(dir_wind,"wind_u.r",sep=""))

# v wind
# This bit works out what size the array should be
store<-0
for (i in 1:filenum){
    infile.v<-paste(dir_wind,"v_",vars[i],".nc",sep="")
    ncdf_windv<-nc_open(infile.v)
    wind.v<-ncvar_get(ncdf_windv)
    d<-dim(wind.v)
    store[i]<-d[3]
}
wind_v<-array(0,dim=c(d[1],d[2],sum(store)))
# This bit stores the data in an array
mn<-1
for (i in 1:filenum)
{
 infile.v<-paste(dir_wind,"v_",vars[i],".nc",sep="")
 ncdf_windv<-nc_open(infile.v)
 wind.v<-ncvar_get(ncdf_windv)
 mx<-mn+store[i]-1
 tp<-paste("min=",mn," max=",mx,sep="")
 print(tp)
 wind_v[,,mn:mx]<-wind.v
 mn<-mn+store[i]
}
save(wind_v,file=paste(dir_wind,"wind_v.r",sep=""))
     