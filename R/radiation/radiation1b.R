library(ncdf4)
library(raster)
library(rgdal)
library(RAtmosphere)
library(insol)
library(sp)

#dir_sis<-"C:/Data2015/CMSAF-SIS/extract/"
#dir_dni<-"C:/Data2015/CMSAF-DNI/extract/"
#dir_sis_out<-"C:/Data2015/CMSAF-SIS/swfiles/"
#dir_dni_out<-"C:/Data2015/CMSAF-DNI/swfiles/"

dir_sis<-"~/Documents/Exeter/Data2015/CMSAF-SIS/extract/"
dir_dni<-"~/Documents/Exeter/Data2015/CMSAF-DNI/extract/"

# FUNCTION - returns number of days in a month for any year between 1980-2015 (accounts for leap years)
monthdays<-function(year,month){
  y=year-1979 # so 1=1980
  feb.d<-c(29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28,29,28,28,28,
           29,28,28,28,29,28,28,28,29,28,
           28,28,29,28,28,28)
  monthdays<-c(31,feb.d[y],31,30,31,30,31,31,30,31,30,31)
  return(monthdays[month])
}

#### FUNCTION for calculating solar zenith angle using sunvector and sunpos in insol package
sza.angle<- function(jd,lt,ln,tz=0){
  sv<-sunvector(jd,lt,ln,tz)
  sza<-sunpos(sv)[,2]
  return(sza)
}

# Function for printing multiple maps of 1+ layers in same stack. 
# Argumennts: stack, vector of layer names, hour (for label) 
plot.stack<-function(stk,lyrs,hour) {
  # Define common colour scheme and scale for plots
  par(mfrow=c(2,2))
  brk<-c(-50,0,100,200,300,400,500,600,700,800,900,1000)
  col<-rev(rainbow(11,start=1/6,end=4/6))
  for (map in 1:length(lyrs)){
    text<-paste(lyrs[map]," at ", hour, ":00",sep="")
    plot(x=stk,lyrs[map],col=col,breaks=brk, main="")
    title(main=text)
  } # for loop
} # function

#######################################################################################
#### Read in Digital Eelevation data - DECISION: what area we want to cover ####
#######################################################################################

#dem<-raster("C:/Data2015/DEM100/dem_sw_x60-420k_y-10-180k.tif")
dem<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif")
plot(dem,main="DEM-full")
extent(dem)
# Crop areas so that it neatly has whole 100m grid cells around the edge and save extent as e.dem for use later
#xmn=floor(e@xmin/100)*100 # CHANGED from Round to Floor - to ensure cropping
#xmx=floor(e@xmax/100)*100
#ymn=floor(e@ymin/100)*100
#ymx=floor(e@ymax/100)*100
#e.dem<-extent(c(xmn,xmx,ymn,ymx))
e.dem<-extent(c(113000,420000,0,174000))
dem<-crop(dem,e.dem)
plot(dem,main="DEM-sw")
e.dem <-extent(dem)
crs(dem)<-"+init=epsg:27700"
plot(dem,main="DEM-OScrs")

#######################################################################################
### Define time period - day/hour loops
#######################################################################################

# Set year and month - could be For loops
year<-2013
month<-6

for (day in 1:monthdays[year,month]) {
  # define arrays for missing value checks
  sismissing<-array(0, dim = c(24))
  sisneg1<-array(0, dim = c(24))
  siszero<-array(0, dim = c(24))
  dnimissing<-array(0, dim = c(24))
  dnineg1<-array(0, dim = c(24))
  dnizero<-array(0, dim = c(24))
  
  for (hr in 0:23) {
    mn<-0
    datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":",sprintf("%02d",mn,sep=""),sep="")
    print(datetime)
    
    
    #############################################################
    ###### 1. INPUT extracted ncdf files to matrix          #####
    ###### And correct any NA to equal '0'                  #####
    #############################################################
    
    # input sis data from extracted ncdf file
    infile.sis<-paste(dir_sis,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),sprintf("%02d",mn,sep=""),"002UD1000101UD.nc",sep="")
    ncdf_sis<-nc_open(infile.sis)
    sis<-ncvar_get(ncdf_sis) # input as matrix
    d.sis<-dim(sis)
    
    # EITHER 1. count missing and zero values and express as % of total elements
   # sismissing[hr+1]<-length(which(is.na(sis))) # calculate number of missing values 
  #  siszero[hr+1]<-length(which(sis==0))# calculate number of 0 values 
   # n<-d.sis[1]*d.sis[2]
  #  q<-array(n, dim = c(24))
   # sismissing[hr+1]<- (sismissing[hr+1]/q[hr+1])*100
    #siszero[hr+1]<- (siszero[hr+1]/q[hr+1])*100
    
    # OR 2.Set  NA values to zero - assumming this to be the case from studied NA/0 values 
    # i.e. NA only during nightime hours etc
    sel<-which(is.na(sis))
    sis[sel]<-0
    
    # input dni data
    infile.dni<-paste(dir_dni,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),sprintf("%02d",mn,sep=""),"002UD1000101UD.nc",sep="")
    ncdf_dni<-nc_open(infile.dni)
    dni<-ncvar_get(ncdf_dni) # input as matrice
    d.dni<-dim(dni)
    
    # EITHER 1. count missing and zero values and express as % of total elements
    #dnimissing[hr+1]<-length(which(is.na(dni))) # calculate number of missing values 
    #dnizero[hr+1]<-length(which(dni==0))# calculate number of 0 values  
    #n<-d.sis[1]*d.sis[2]
    #q<-array(n, dim = c(24))
    #dnimissing[hr+1]<- (dnimissing[hr+1]/q[hr+1])*100
    #dnizero[hr+1]<- (dnizero[hr+1]/q[hr+1])*100
    
    # OR 2.Set  NA values to zero - assumming this to be the case from studied NA/0 values 
    # i.e. NA only during nightime hours etc
    sel<-which(is.na(dni))
    dni[sel]<-0
    
    #############################################################
    ###### 2. Convert to long/lat res raster                  
    ######  Calculate SID from DNI and solar zenith angle
    #############################################################
    
    #### Convert to long/lat raster  ###
    res=c(0.05,0.05) # Resolution of sis/dni data = 0.05 degrees
    sis.r<-raster(sis,xmn=-20,xmx=20,ymn=30,ymx=65,res ) # set max/min values to match extent of sis/dni data
    crs(sis.r)<-"+init=epsg:4326"
    
    # repeat for dni
    res=c(0.05,0.05)
    dni.r<-raster(dni,xmn=-20,xmx=20,ymn=30,ymx=65,res ) # change max/min values ,"+proj=longlat+datum=WGS84"
    crs(dni.r)<-"+init=epsg:4326"
    
    
    #############################################################
  
    # Create vectors of long/lat and Julian times
    cells<-ncell(dni.r)
    d<-dim(dni.r)
    # longitudes - convert -ve
    ln<-xFromCell(dni.r,c(1:cells))
        sel<-which(ln<0)
    ln[sel]<-ln[sel]+360
    # lattitude
    lt<-yFromCell(dni.r,c(1:cells))
    # julian day/time
    jday<-JDymd(year,month,day,hr)
    jdays<-rep(jday,length(ln))  # Convert datetime to Julian time and hr to decimal fraction
    
    # Calculate SZA and create raster of SZA values
    sza<-sza.angle(jdays,lt,ln)/360  # in decial degrees!!!
    #csza<-cos(sza)
    
    sza.m<-matrix(sza,nrow=nrow(dni.r),byrow=TRUE) # convert array to matrix
    #csza.m<-matrix(csza,nrow=nrow(dni.r),byrow=TRUE) # convert array to matrix
    
    sza.r<-raster(sza.m, template=dni.r) # convert matrix to raster
    #csza.r<-raster(csza.m, template=dni.r) # convert matrix to raster

    # Calculate SID as dni * cosine(sza)
    sid.r<-overlay(x = dni.r,y = sza.r,fun = function(x,y) {
      z<-x*cos(y) 
      return(z)
    } ) 
  
    # Calculate indirect from SIS and SID
    ind.r<-overlay(x = sis.r,y = sid.r,fun = function(x,y) {
        z<-x-y 
        return(z)
    } ) 
  
    dif.r<-overlay(x = dni.r,y = sid.r,fun = function(x,y) {
      z<-x-y 
      return(z)
    } ) 
    dif2.r<-overlay(x = sis.r,y = dni.r,fun = function(x,y) {
      z<-x-y 
      return(z)
    } )  
    
    #############################################################
    ###### 3. Resample to 100m OS raster using stack
    ######   Crop to DEM and convert all sea cells to NA                                               #####
    #############################################################
    # Derive template from dem
    xmn=round(e.dem@xmin/100)*100
    xmx=round(e.dem@xmax/100)*100
    ymn=round(e.dem@ymin/100)*100
    ymx=round(e.dem@ymax/100)*100
    template<-raster(xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,resolution=100)
    
    # Convert to raster layers of single stack
    stack.r<-stack(sis.r,dni.r,sid.r,ind.r)
    #names(stack.r)<-c("sis","dni","sid","ind")
    # change co-ord to OS 
    stack.r<-projectRaster(stack.r,crs="+init=epsg:27700")
    #Resample and crop to dem 100metre template
    stack.r<-resample(stack.r,template)# use same template calulated for sis
    # Convert all sea cells to NA for all layers (based on dem)    
    stack.r <- overlay(x = dem,y = stack.r,fun = function(x,y) {
     y[is.na(x)] <- NA
     return(y)
    } )   
  
    stack.r<-addLayer(stack.r,dem)
    names(stack.r)<-c("sis","dni","sid","ind","dem")
  
  # Plot and save Layers
  plot.stack(stack.r,c("sis","dni","sid","ind"),hr)
  minvals<-paste("Minimum values - sis: ", minValue(sis.r), " dni: ", minValue(dni.r)," sid: ",minValue(sid.r)," ind: ",minValue(ind.r))
  maxvals<-paste("Maximum values - sis: ", maxValue(sis.r), " dni: ", maxValue(dni.r)," sid: ",maxValue(sid.r)," ind: ",maxValue(ind.r))
  print(minvals)
  print(maxvals)
  
  } # end of hour FOR loop
  
  #print("SIS missing: "); sismissing
  #print("SIS zero: "); siszero
  #print("DNI missing: "); dnimissing
  #print("DNI zero: "); dnizero
  
  
} # end of day FOR loop
    
    









######### PLOTS #############
par(mfrow=c(2,2))
plot(dem, main="dem")
plot(sis_100_land, main="sis")
plot(dni_100_land, main="dni")
plot(indirect_100_land, main="sis-dni")

#############################

    

######### PLOT TIME SERIES #############
par(mfrow=c(1,1))
year<-2013; month<-6; day<-30

# calculate and save plot per hr

# Define common colour scheme and scale for plots
brk<-c(0,100,200,300,400,500,600,700,800,900,1000)
col<-rev(rainbow(11,start=1/6,end=4/6))
#col<-rev(heat.colors(11))

# save plots as list for one day
sis.plot<-list(0:23)
for (hr in 0:23) {
  file.name<-paste(dir_sis_out,"sis_",year,"_",month,"_",day,"_",hr,".r",sep="")
  load(file=file.name)
  sis.r<-raster(sis.out, template=dem)
  t<-paste("SIS ",day,"-",month,"-",year," Hour: ", hr, ":00",sep="")
  sis.plot[hr+1]<-plot(sis.r, col=col,breaks=brk,main=t)
  }


# print plot per hr if required
for (hr in 0:23) {
  sis.plot[hr+1]
}


# plot( d, col=rev( rainbow( 99, start=0,end=1 ) ), breaks=seq(min(minValue( d )),max(maxValue(d)),length.out=100) )
  
#####################################################################
# Write outputs file for every hour - re-write for using STACK
#####################################################################

sis.out<-getValues(sis_100_land,format="matrix")
fileout.sis<-paste(dir_sis_out,"sis_",year,"_",month,"_",day,"_",hr,".r",sep="")
save(sis.out,file=fileout.sis)

dni.out<-getValues(dni_100_land,format="matrix")
fileout.dni<-paste(dir_dni_out,"dni_",year,"_",month,"_",day,"_",hr,".r",sep="")
save(dni.out,file=fileout.dni)

sid.out<-getValues(dni_100_land,format="matrix")
fileout.sid<-paste(dir_sid_out,"sid_",year,"_",month,"_",day,"_",hr,".r",sep="")
save(sid.out,file=fileout.sid)

#####################################################################




### IDEAS  ####

# Calculate gdd for 5km cells over a year - run as series of rasters ---> movie
# Could add predicted flowering etc !
# Calculate gdd and risks for different 'types' of growing season - good/average/poor
# to compare spatial variation under different types of seasonal conditions


















### CUT OUTS
# Reproject in OSGB projection
sis_osgb<-projectRaster(sis.r,crs="+init=epsg:27700")

# Interpolate to 100m cells and crop to fit dem.Default: simple bilinear
sis_100<-resample(sis_osgb,template)

sis_100<-crop(sis_100,e.dem) 

# Set sea cells (NA in dem) to NA 
sis_100_land <- overlay(x = dem,y = sis_100,fun = function(x,y) {
  y[is.na(x)] <- NA
  return(y)
} )  

# repeat for dni - DISCARD later
dni_osgb<-projectRaster(dni.r,crs="+init=epsg:27700")

dni_100<-resample(dni_osgb,template)# use same template calulated for sis
dni_100<-crop(dni_100,e.dem) 
dni_100_land <- overlay(x = dem,y = dni_100,fun = function(x,y) {
  y[is.na(x)] <- NA
  return(y)
} )   

# repeat for SID 
sid_osgb<-projectRaster(sid.r,crs="+init=epsg:27700")

sid_100<-resample(sid_osgb,template)# use same template calulated for sis
sid_100<-crop(sid_100,e.dem) 
sid_100_land <- overlay(x = dem,y = sid_100,fun = function(x,y) {
  y[is.na(x)] <- NA
  return(y)
} ) 


#######################################





##########################################################
## Unpack ncgz files for a day and save ncdf files to different directory   
##########################################################


# connect and read gz file
dir_dni_ncgz<-"C:/Data2015/CMSAF-DNI/"

dir_dni_ncgz<-"C:/Data2015/CMSAF-DNI/ncgz/"
dir_dni_ncdf<-"C:/Data2015/CMSAF-DNI/swfiles"
dni.extractfile<-paste(dir_dni_ncdf,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),min,"002UD1000101UD.nc",sep="")

# set date/time variables for file to be unpacked
year<-2000

# Loop here ?
# open file
dni.gzfile<-paste(dir_dni_ncgz,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),min,"002UD1000101UD.nc.gz",sep="")
print(dni.gzfile)

unzip(dni.gzfile,exdir=dir_dni_ncdf)
unzip("C:/Data2015/CMSAF-DNI/DNIhm200006301300002UD1000101UD.nc.gz", exdir="C:/Data2015/CMSAF-DNI/swfiles",unzip="C:/Program Files/7-Zip/7z.exe")

gzfile("C:/Data2015/CMSAF-DNI/DNIhm200006301300002UD1000101UD.nc.gz")
gzfile(description, open = "", encoding = getOption("encoding"),
       compression = 6)


system("7z C:/Data2015/CMSAF-DNI/DNIhm200006301300002UD1000101UD.nc.gz")




#gunzip(dni.gzfile)
#gunzip((dni.gzfile, destname=gsub("[.]gz$", "", filename, ignore.case = TRUE),
temporary=FALSE))

#untar(dni.gzfile,list=TRUE, compressed=TRUE, exdir="../swfiles")

test<-gzfile(dni.gzfile)
open(test,"rb")
test.sis<-readLines(test)
close(test)
unlink(infile.sis)

# Trim to desired area and resample the data at a 100m resolution. Interpolation is set at default simple bilinear
e<-extent(sis_osgb) 
xmn=round(e@xmin/100)*100
xmx=round(e@xmax/100)*100
ymn=round(e@ymin/100)*100
ymx=round(e@ymax/100)*100
template<-raster(xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,resolution=100)
sis_100<-resample(sis_osgb,template)

sis_100_b<-crop(sis_100,e.dem) # crop to match DEM map

### PLOTS
image(dem, zlim=c(1,300)) # low lying areas
image(dem, zlim=c(-10,1337))
image(dem, zlim=is.na)


########################################



JD<-function(day,month,year){
  a<-(14-month)/12
  y<-year+4800-a
  m<-month+12*a-3
  JDN<-floor((153*m+2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045 + day
  JDN
}

# calculates sunrise and sunset
suntimes<-function(day,month,year,Lat,Long,Timezone=0,DST=0){
  J<-JD(day,month,year)
  lw<-Long*-1 # converts all longitude to -ve assuming West of GMT
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
