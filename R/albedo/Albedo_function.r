args <-commandArgs(trailingOnly = TRUE)
print(args)
image <- args[1] 


# FUNCTIONS USED
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
  if (landsat==5|landsat==7){
    blue.range=c(450,520)
    green.range=c(520,600)
    red.range=c(630,690)
    NIR.range=c(770,900)
    
    # Landsat 5, 7 band numbers: blue (1), green(2)), red(3), NIR(4)
    Blue<-raster(paste(dir_landsat,image.list[n],1,".tif",sep=""))
    Green<-raster(paste(dir_landsat,image.list[n],2,".tif",sep=""))
    Red<-raster(paste(dir_landsat,image.list[n],3,".tif",sep=""))
    NIR<-raster(paste(dir_landsat,image.list[n],4,".tif",sep=""))
  }
  if (landsat==8){
    blue.range=c(450,510)
    green.range=c(520,590)
    red.range=c(640,670)
    NIR.range=c(850,880)
    
    # Landsat 8 band numbers: blue (2), green(3), red(4), NIR(5)
    Blue<-raster(paste(dir_landsat,image.list[n],2,".tif",sep=""))
    Green<-raster(paste(dir_landsat,image.list[n],3,".tif",sep=""))
    Red<-raster(paste(dir_landsat,image.list[n],4,".tif",sep=""))
    NIR<-raster(paste(dir_landsat,image.list[n],5,".tif",sep=""))
  }
  
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
  albedo.r<-(Blue*rat.blue)+(Green*rat.green)+(Red*rat.red)+(NIR*rat.NIR)
  #albedo.r<-overlay(Blue,Green,Red,NIR,fun=function(w,x,y,z)
   # {w*rat.blue+x*rat.green+y*rat.red+z*rat.NIR})
  max.val<-10000
  albedo.r<-albedo.r/max.val  
  return(albedo.r)
}

#########################################################

# Select landsat image
dir_landsat<-"/Data2015/albedo/landsat/extract/"
#dir_landsat<-"~/Documents/Exeter/Data2015/Albedo/landsat/"
#image<-"LT52040242011118KIS00_sr_band"

print(paste("In files from image: ", image,sep=""))
landsat<-substr(image,3,3)
outfile<-paste(dir_landsat,substr(image,1,22),"albedo.tif",sep="")
print(paste("Out file= ", outfile,sep=""))

albedo.r<-albedo(image,landsat)
par(mfrow=c(1,1))
plot(albedo.r,main="Albedo")

writeRaster(albedo.r,file=outfile,overwrite=TRUE)

# load metadata
#con<-file(paste(dir_landsat,image.list[n],"_MTL.txt",sep=""),'r' )
#input<- readLines(con, n=-1)
#image.date<-substr(input[21],21,30)
# Get coords for image - x,y of each corner???
#image.e<-c(substr(input[31],29,37),substr(input[33],29,37),substr(input[34],29,37),substr(input[],29,37))

