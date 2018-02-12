# Downscale SST
# Input: single ncdf file of Met Office Hadley Centre HadISST1 data 
# Output: Downscaled SST to hourly 100m data
# Issues/Problems: interpolation/resampling - in effect across N/S coasts

library(R.utils) # IMPORTANT - hides raster functions requiring raster::
library(ncdf4)
library(raster)
library(rgdal)
####################################################

# Input directories
#dir_sst<-"~/Documents/Exeter/Data2015/sst/"
dir_sst<-"C:/Data2015/SST/"  
#gunzip(paste(dir_sst,"HadISST_sst.nc.gz",sep=""))
dir_sstm<-paste(dir_sst,"monthly/",sep="")
dir_ssth<-paste(dir_sst,"hourly/",sep="")
dir_grids<-"~/Documents/Exeter/Data2015/Templates/"
#dir_grids<-"C:/Data2015/Templates/"

####################################################
# Read in Digital Eelevation data - DECISION: what area to cover
# **** ENSURE here that dem is rounded to 10km ****

#dem<-raster("C:/Data2015/DEM100/dem_sw_x60-420k_y-10-180k.tif")
dem<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif", crs=("+init=epsg:27700"))
#plot(dem,main="DEM-full")
extent(dem)
#e.dem<-extent(c(70000,420000,0,180000)) # includes scilly isles
e.dem<-extent(c(120000,420000,0,180000)) # excludes scilly isles
dem<-crop(dem,e.dem)
plot(dem,main="DEM-sw")

# set buffer eg to 10km
buffer<-30000
e.buf<-extent(xmin(e.dem)-buffer,xmax(dem)+buffer,ymin(dem)-buffer,ymax(dem)+buffer)

# Set ukcp gridcell frame and extent of dem
grid.file<-paste(dir_grids,"ukhistmask.grd",sep="") # = historic grid (less land)
#grid.file<-paste(dir_grids,"ukcpmask.grd",sep="") # =ukcp09 grid

print(grid.file)
gridmask.r<-raster(grid.file)
#plot(gridmask.r)

dsgrid.r<-crop(gridmask.r,e.buf)
# To convert to 100m resolution
#dsgrid.r<-disaggregate(dsgrid.r,fact=c(50,50))
plot(dsgrid.r)

####################################################
# Functions to return year/month from sst t value
ttoyr<-function(t){
  yr<-ceiling(t/12)-1+1960
  return(yr)
}

ttomonth<-function(t){
  month<-t-((ttoyr(t)-1960)*12)
  return(month)
}

# Function to spatially downscale sst data for month t to dsgrid  
sst.spdownsc<-function(t,dsgrid){
      # Read ncdf direct to raster
      sst.world.r<-raster(in.file,band=t) # read all of level t of ncdf file
      #plot(sst.world.r)
      # Extract relevant 1 deg cells 353-360E, 48:53N, 1960-2015 
      long.mn<--7;long.mx<--1 # -7 to -1
      lat.mn<-49; lat.mx<-52 # 49 to 52
      e.sst<-extent(c(long.mn,long.mx,lat.mn,lat.mx))
      sst.r<-crop(sst.world.r,e.sst)
      #plot(sst.r,main="t")
      #plot(sst.r,main=paste("SST long: ",long.mn," to ",long.mx,", lat: ",lat.mn," to ",lat.mx,sep=""))

      # Set land cells to have same temperature as nearby sea
      #sst.r<-focal(sst.r,w=matrix(1,3,3),fun=mean,na.rm=TRUE,NAonly=TRUE)
      sst.m<-getValues(sst.r,format="matrix")
      sst.m[1,5]<-(sst.m[1,4]) # Bristol channel - count as sea = to adjacent cell (NOT mean including CHannel!)
      #sst.m[1,6]<-(sst.m[2,6]+sst.m[1,5])/2 # keep as NA all land cells
      sst.r<-raster(sst.m,template=sst.r)
      #plot(sst.r)
      
      # Convert projection to OSGB
      projection(sst.r)<-"+init=epsg:4326"
      sst.os<-projectRaster(sst.r,crs="+init=epsg:27700")

      # resample to 5km resolution - cf with ngb method
      sst5km.r<-raster::resample(sst.os,dsgrid.r,method="bilinear",na.rm=TRUE)
      #plot(sst5km.r)
      sst5km.r<-mask(sst5km.r,dsgrid.r,inverse=TRUE)
      #plot(sst5km.r)
      
      return(sst5km.r)
     
}# end of function

##########################################################################################
# Spatial downscaling of sst data to 5km grid cells matching ukcp

in.file<-paste(dir_sst,"HadISST_sst.nc",sep="")
ncfile<-nc_open(paste(dir_sst,"HadISST_sst.nc",sep="")) # summary of file variables,  dimensions, attributes

# Set time parameters
start.yr <- 2010 
start.month<-1
end.yr<-2010
end.month<-12

#time.mn<-((1960-1870)*12)
#time.mx<-d[3]
start.t<-(start.yr-1960)*12+start.month
end.t<-((end.yr-1960)*12)+end.month
t<-start.t

for (t in start.t:end.t){
  sst5km.r<-sst.spdownsc(t,dsgrid.r)
  tl<-paste("Year: ",ttoyr(t)," Month: ",ttomonth(t),sep="")
  plot(sst5km.r,main=tl)
  fileout<-paste(dir_sstm,"sst_",ttoyr(t),"_",ttomonth(t),".grd",sep="")
  print(fileout)
  writeRaster(sst5km.r,file=fileout,overwrite=T)
} 


##########################################################################################
# 2. TEMPORAL DOWNSCALING to DAILY/HOURLY time period
# Input next month data as above
# For loop for writing hourly sst data

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
  filein1<-paste(dir_sstm,"sst_",year,"_",mth1,".grd",sep="")
  filein2<-paste(dir_sstm,"sst_",year,"_",mth2,".grd",sep="")
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
 
return(r.hr)
} # end function


year<-2010
month<-6


for (month in 1:12)
date<-as.Date(paste(year,"-",month,"-",1,sep=""),"%Y-%m-%d")
days<-nday.in.month(date)


for (day in 1:days){
    hr<-0
    sst.hr.r<-sst.time.int(year,month,day,hr)
    #sst.hr.r<-mask(sst.hr.r,gridmask.r)
    tl<-paste("Year: ",year," Month: ",month," Day: ",day," Hour: ",hr,sep="")
    plot(sst.hr.r,main=tl)
    fileout<-paste(dir_ssth,"sst_",year,"_",month,"_",day,"_00h.grd",sep="")
    print(fileout)
    writeRaster(sst.hr.r,file=fileout,overwrite=T)
  } # end day

} # end month






# CUT OUTS
# Crop to dem and check land cells (!NA in dem) = NA 
sst.sw.r<-crop(sst100.r,e.dem)
sst.sw.r <- overlay(x = dem,y = sst.sw.r,fun = function(x,y) {
  y[!is.na(x)] <- NA
  return(y)
} ) 
plot(sst.sw.r)




# ALTERNATIVE SPATIAL DOWNSCALING BUT not rounded!!!
sst100.r<-projectRaster(sst.r,crs="+init=epsg:27700", res=c(100,100), method='bilinear')
plot(sst100.r) 
sst.sw.r<-crop(sst100.r,e.dem)
sst.sw.r <- overlay(x = dem,y = sst.sw.r,fun = function(x,y) {
  y[!is.na(x)] <- NA
  return(y)
} ) 


xdisag<-res(sst.os)[1]/5000
ydisag<-res(sst.os)[2]/5000
sst5km.r<-disaggregate(sst.os,fact=c(xdisag,ydisag), method="bilinear")

e.grid<-extent(xmin(dsgrid),xmax(dsgrid),ymin(dsgrid),ymax(dsgrid))
sst5km.r<-resample(sst.os,e.grid,ncol=ncol(dsgrid),nrow=nrow(dsgrid))


