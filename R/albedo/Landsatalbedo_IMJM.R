# Calculate albedo of Landsat image
library(rgdal)
library(raster)
# # # # # # # # # # # # # #  # #
# Calculate Albedo from Landsat
# # # # # # # # # # # # #  # # #
# # # # # Calculates Black body radiation as a function of wavelength
Planck<-function(delta) # delta is Wavelength of radiance (nm)
{
  h=6.62606957*10^(-34)
  cc=299792458
  T=5250+273.16
  k=1.3806488*10^(-23)
  B=2*h*cc^2/delta^5*(1/(exp((h*cc)/(k*T*delta))-1))
  B
}
# Band1..7 = pixel values for individual bands
# max.val  =  maximum value (255)
# band1..7.range = wavelengths (in nanometers) of individual bands
# Check vs: https://landsat.usgs.gov/what-are-band-designations-landsat-satellites
# Max value = 10000 from SR products
albedo<-function(Band1,Band2,Band3,Band4,Band5,Band7,max.val=10000,
                 band1.range=c(450,520),
                 band2.range=c(520,600),
                 band3.range=c(630,690),
                 band4.range=c(770,900),
                 band5.range=c(1055,1750),
                 band7.range=c(2090,2350))
{
  d<-c(1:3000)*10^(-9)
  p<-Planck(d)
  # Calculates weighting for each band
  weight.1<-mean(p[band1.range[1]:band1.range[2]])
  weight.2<-mean(p[band2.range[1]:band2.range[2]])
  weight.3<-mean(p[band3.range[1]:band3.range[2]])
  weight.4<-mean(p[band4.range[1]:band4.range[2]])
  weight.5<-mean(p[band5.range[1]:band5.range[2]])
  weight.7<-mean(p[band7.range[1]:band7.range[2]])
  weight.all<-weight.1+weight.2+weight.3+weight.4+weight.5+weight.7
  rat.1<-weight.1/weight.all
  rat.2<-weight.2/weight.all
  rat.3<-weight.3/weight.all
  rat.4<-weight.4/weight.all
  rat.5<-weight.5/weight.all
  rat.7<-weight.7/weight.all
  albedo<-Band1*rat.1+Band2*rat.2+Band3*rat.3+Band4*rat.4+Band5*rat.5+Band7*rat.7
  albedo<-albedo/max.val
  albedo
}
# # # # #  Calculate albedo
root<-"~/Documents/Exeter/Data2015/"
dir_landsat<-paste(root,"Albedo/landsat/extract/",sep="")
dir_albedo<-paste(root,"Albedo/landsat/albedo/",sep="")

dir.in<-dir_landsat

# nm<-"LE72040251999205AGS01" ; # for testing without a loop

dirs<-list.files(dir_landsat)
#dirs<-substring(filenames,1,21)
#nms<-unique(nms)
print(dirs)
bandin<-c(1,2,3,4,5,7)

for (n in 1:length(dirs)){
    dir.in<-paste(dir_landsat,dirs[n],"/",sep="")
    nms<-list.files(dir.in)
    nm<-substr(nms[1],1,40)
    print(paste(dir.in,nm,sep=""))
    b1<-getValues(raster(paste(dir.in,nm,"_sr_band",bandin[1],".tif",sep="")),format="matrix")
    b2<-getValues(raster(paste(dir.in,nm,"_sr_band",bandin[2],".tif",sep="")),format="matrix")
    b3<-getValues(raster(paste(dir.in,nm,"_sr_band",bandin[3],".tif",sep="")),format="matrix")
    b4<-getValues(raster(paste(dir.in,nm,"_sr_band",bandin[4],".tif",sep="")),format="matrix")
    b5<-getValues(raster(paste(dir.in,nm,"_sr_band",bandin[5],".tif",sep="")),format="matrix")
    b7<-getValues(raster(paste(dir.in,nm,"_sr_band",bandin[6],".tif",sep="")),format="matrix")
    alb<-albedo(b1,b2,b3,b4,b5,b7)
    
    # Convert to raster
    r<-raster(paste(dir.in,nm,"_sr_band",bandin[1],".tif",sep="")) # # raster template
    r.alb<-raster(alb,template=r) # albedo as raster
    plot(r.alb,main=nm)
    writeRaster(r.alb,file=paste(dir_albedo,nm,"-albedo-v2.tif",sep=""),overwrite=T)
}
################################
# Reset masked cells to NA values - NOT USED
################################
# converts any pixels with value not equal to 0 in mask to NA in raster
cloud2NA <- function(x, y){
  #x[y>1] <- NA # leaves water pixels
  x[y != 0] <- NA # sets all flags to NA incl water pixels
  return(x)
}

trim.image<-function(image.r)
{
  image.vals<-getValues(image.r,format="matrix")
  for (n in 1: ncol(image.r)){
    NonNAindex <- which(!is.na(image.vals[,n]))
    if (length(NonNAindex)!=0){
      firstNonNA <- min(NonNAindex)
      lastNonNA<-max(NonNAindex)
      x1<-firstNonNA+180
      x2<-lastNonNA-280
      image.vals[1:x1,n]<-NA
      image.vals[x2:lastNonNA,n]<-NA
    }
  } # end for
  result.r<-raster(image.vals,template=image.r)
  return(result.r)
} # end function

dir_qa<-paste(dir_albedo,"pixelqa/",sep="")

for (n in 1:length(dirs)){
  dir.in<-paste(dir_landsat,dirs[n],"/",sep="")
  nms<-list.files(dir.in)
  nm<-substr(nms[1],1,40)
  #print(paste(dir.in,nm,sep=""))
  # Apply cloud mask
  #nm<-nms[n]
  print(nm)
  albedo.r<-raster(paste(dir_albedo,nm,"-albedo-v2.tif",sep=""))
  
  # Load cloud mask
  #cloudmask<-paste(dir_qa,nm,"_pixel_qa.tif",sep="")
  #print(cloudmask)
  #cloudmask.r<-raster(cloudmask)
  #plot(cloudmask.r)
  
  new.r<-trim.image(albedo.r)
  
  # Apply cloud mask
  #new.r<-overlay(x=albedo.r, y=cloudmask.r, fun = cloud2NA)
  plot(new.r,main=nm)
  # Overwrite old with new raster
  writeRaster(new.r,file=paste(dir_albedo,nm,"-albedo-v2.tif",sep=""),overwrite=T)
  
}

################################

################################
# Reproject to OSGB and coarsen to 100m after cropping to someting sensible
# NB saves a bit of memory if you open and close R here
################################
library(rgdal)
library(raster)


for (n in 1:length(dirs)){
  dir.in<-paste(dir_landsat,dirs[n],"/",sep="")
  nms<-list.files(dir.in)
  nm<-substr(nms[1],1,40)
  print(nm)

  r.alb<-raster(paste(dir_albedo,nm,"-albedo-v2.tif",sep=""))
  rp<-projectRaster(r.alb, crs="+init=epsg:27700")
  rp<-shift(rp,x=-8.88,y=9.58)
  
  e.dem<-extent(c(60000,420000,0,180000 )) # includes scilly isles
  #e<-extent(c(73200,325200,0,160200))
  rp<-crop(rp,e.dem)
  plot(rp,main="Reprojected & cropped - 30m")
  
  # create template for resample
  #m<-array(0,dim=c(1602,2520))
  #m<-array(0,dim=c(1800,3600))
  #r<-raster(m,xmn=73200,xmx=325200,ymn=0,ymx=160200)
  rp<-resample(rp,dem)
  plot(rp,main=nm)
  
  writeRaster(rp,file=paste(dir_albedo,nm,"-albedo100.tif",sep=""), overwrite=T)

}


################################


###################
# CREATE SINGLE MAP
###################

albedo.map<-raster(ext=e.dem,res=c(100,100),crs="+init=epsg:27700")

for (n in 1:length(nms)){
  nm<-nms[n]
  print(nm)
  r<-raster(paste(dir_albedo,nm,"-albedo100.tif",sep=""))

  # Crop then extend to fit into area of interest e.dembuf
  new.r<-extend(r,e.dem) 
  albedo.map<-mosaic(albedo.map,new.r,fun=mean)
  plot(albedo.map,main=paste("Albedo map from ",n," images",sep=""))
  remove(new.r,cloudmask.r,albedo.r)
}
outfile<-paste(dir_albedo,"albedomap_mean.tif",sep="")
print(outfile)
writeRaster(albedo.map,file=outfile,format="GTiff",overwrite=TRUE)

## Set to default value of 0.2 any cells without values and not sea
albedo.map[is.na(albedo.map)]<-0.2
albedo.map<-mask(albedo.map,dem)

