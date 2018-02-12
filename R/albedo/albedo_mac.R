args <-commandArgs(trailingOnly = TRUE)
image <- args[1] 
print(image)

options(rasterTmpDir='/home/ISAD/jm622/rscripts/temporary/')

source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

dir_landsat<-"/home/ISAD/jm622/Data2015/Albedo/landsat/extract/"
dir_albedo<-"/home/ISAD/jm622/Data2015/Albedo/landsat/albedo/"
#library(raster)

###############################################
# FUNCTIONS USED
# converts any pixels with value not equal to 0 in mask to NA in raster
cloud2NA <- function(x, y){
  #x[y>1] <- NA # leaves water pixels
  x[y != 0] <- NA # sets all flags to NA incl water pixels
  return(x)
}

Planck<-function(delta)
{
  h=6.62606957*10^(-34)
  cc=299792458
  T=5250+273.16
  k=1.3806488*10^(-23)
  B=2*h*cc^2/delta^5*(1/(exp((h*cc)/(k*T*delta))-1))
  B
}

# ranges of bands from: http://landsat.usgs.gov//band_designations_landsat_satellites.php for Landsat 5,7  
albedo<-function(image,landsat)
{ 
  # Define bands and load band data
  image<-paste(image,"_sr_band",sep="")
  
  if (landsat==5|landsat==7){
    blue.range=c(450,520)
    green.range=c(520,600)
    red.range=c(630,690)
    NIR.range=c(770,900)
    
    print(paste("Load rasters ",dir_landsat,image,1,".tif"," ",
                dir_landsat,image,2,".tif"," ", dir_landsat,image,3,".tif"," ",
                dir_landsat,image,4,".tif",sep=""))
    
    # Landsat 5, 7 band numbers: blue (1), green(2)), red(3), NIR(4)
    Blue<-raster(paste(dir_landsat,image,1,".tif",sep=""))
    Green<-raster(paste(dir_landsat,image,2,".tif",sep=""))
    Red<-raster(paste(dir_landsat,image,3,".tif",sep=""))
    NIR<-raster(paste(dir_landsat,image,4,".tif",sep=""))
  }
  if (landsat==8){
    blue.range=c(450,510)
    green.range=c(520,590)
    red.range=c(640,670)
    NIR.range=c(850,880)
    
    print(paste("Load rasters ",dir_landsat,image,2,".tif"," ",
                dir_landsat,image,3,".tif"," ", dir_landsat,image,4,".tif"," ",
                dir_landsat,image,5,".tif",sep=""))
    
    # Landsat 8 band numbers: blue (2), green(3), red(4), NIR(5)
    Blue<-raster(paste(dir_landsat,image,2,".tif",sep=""))
    Green<-raster(paste(dir_landsat,image,3,".tif",sep=""))
    Red<-raster(paste(dir_landsat,image,4,".tif",sep=""))
    NIR<-raster(paste(dir_landsat,image,5,".tif",sep=""))
  }
  print("Calculate weightings..")
  d<-c(1:3000)*10^(-9)
  p<-Planck(d)
  # Calculates weighting for each band
  weight.blue<-sum(p[blue.range[1]:blue.range[2]])
  weight.green<-sum(p[green.range[1]:green.range[2]])
  weight.red<-sum(p[red.range[1]:red.range[2]])
  weight.NIR<-sum(p[NIR.range[1]:NIR.range[2]])
  weight.all<-weight.blue+weight.green+weight.red+weight.NIR
  rat.blue<-weight.blue/weight.all
  rat.green<-weight.green/weight.all
  rat.red<-weight.red/weight.all
  rat.NIR<-weight.NIR/weight.all
  print("CAlculate albedo...")
  #albedo.r<-(Blue*rat.blue)+(Green*rat.green)+(Red*rat.red)+(NIR*rat.NIR)
  albedo.r<-overlay(Blue,Green,Red,NIR,fun=function(w,x,y,z){w*rat.blue+x*rat.green+y*rat.red+z*rat.NIR})
  max.val<-10000
  albedo.r<-albedo.r/max.val  
  return(albedo.r)
}

#########################################################

# Select landsat image
dir_landsat<-"~/Documents/Exeter/Data2015/Albedo/landsat"
#dir_landsat<-"C:/Data2015/Albedo/landsat/extract/"
dir_albedo<-"~/Documents/Exeter/Data2015/Albedo/landsat_results/"
#image<-"LT52040242011118KIS00_sr_band"
# Test files
#dir_landsat<-"C:/Data2015/Albedo/landsat/2013-04-06/"
#image<-"LE72030242003081SGS00"
#image<-"LE72030252005198EDC00"
#image<-"LE72030252007156EDC00"
#image<-"LT52040242011118KIS00"
#image<-"LC82030242014103LGN00"
#image<-"LC82030252014183LGN00"
image<-"LC82030252014183LGN00"
#########################################################

# Define input output files
print(paste("In files from image: ", image,sep=""))
landsat<-substr(image,3,3)
print(landsat)

# Load cloud mask - for selection/elimation of pixels
cloudfile<-paste(dir_landsat,image,"_cfmask.tif",sep="")
cloudmask.r<-raster(cloudfile)

# Calculate albedo from bands
albedo.r<-albedo(image,landsat)
plot(albedo.r,main=paste("Albedo of image ",image,sep=""))



# Create stack
albedo.s<-stack(albedo.r,cloudmask.r)
names(albedo.s)<-c("albedo","cloudmask")



# Write raster and cloudmask
outfile<-paste(dir_albedo,substr(image,1,22),"_albedo.tif",sep="")
print(paste("Out file= ", outfile,sep=""))
writeRaster(albedo.osgb[1],file=outfile,overwrite=TRUE)

outfile<-paste(dir_albedo,substr(image,1,22),"_cloudmask.tif",sep="")
print(paste("Out file= ", outfile,sep=""))
writeRaster(albedo.osgb[2],file=outfile,overwrite=TRUE)


remove(albedo.r,albedo.osgb)



ndvi <- overlay(GewataB4, GewataB3, fun=function(x,y){(x-y)/(x+y)})