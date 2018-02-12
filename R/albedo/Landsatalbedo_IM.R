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
albedo<-function(Band1,Band2,Band3,Band4,Band5,Band7,max.val=255,
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

n1<-"LE72040251999205AGS01"
n2<-"L72204025"

nms<-c(n1,n1,n1,n1,n1,n1)
bandin<-c(10,20,30,40,50,70)
b1<-getValues(raster(paste(dir.in,nms[1],"_sr" 02520100604_B",bandin[1],".TIF",sep="")),format="matrix")
b2<-getValues(raster(paste(dir.in,nms[2],"_02520100604_B",bandin[2],".TIF",sep="")),format="matrix")
b3<-getValues(raster(paste(dir.in,nms[3],"_02520100604_B",bandin[3],".TIF",sep="")),format="matrix")
b4<-getValues(raster(paste(dir.in,nms[4],"_02520100604_B",bandin[4],".TIF",sep="")),format="matrix")
b5<-getValues(raster(paste(dir.in,nms[5],"_02520100604_B",bandin[5],".TIF",sep="")),format="matrix")
b7<-getValues(raster(paste(dir.in,nms[6],"_02520100604_B",bandin[6],".TIF",sep="")),format="matrix")
alb<-albedo(b1,b2,b3,b4,b5,b7)

# Convert to raster
r<-raster(paste(dir.in,nms[1],"_02520100604_B",bandin[1],".TIF",sep=""))  # # raster template
r.alb<-raster(alb,template=r) # albedo as raster
plot(r.alb)
writeRaster(r.alb,file="C:/Jonathanmodel/Landsat/albedo.tif",overwrite=T)

# Reproject to OSGB and coarsen to 100m after cropping to someting sensible
# NB saves a bit of memory if you open and close R here
library(rgdal)
library(raster)
r.alb<-raster("C:/Jonathanmodel/Landsat/albedo.tif")
rp<-projectRaster(r.alb, crs="+init=epsg:27700")
rp<-shift(rp,x=-8.88,y=9.58)
e<-extent(c(73200,325200,0,160200))
rp<-crop(rp,e)
# create template for resample
m<-array(0,dim=c(1602,2520))
r<-raster(m,xmn=73200,xmx=325200,ymn=0,ymx=160200)
rp<-resample(rp,r)
writeRaster(rp,file="C:/Jonathanmodel/Landsat/albedo100.tif",overwrite=T)










