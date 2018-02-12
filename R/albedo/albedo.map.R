# Creates single map of albedo from mosaic of albedo images dervied from landsat data
# Uses cloudmask 
# Mosaic function when overlapping images = MEAN
# Projection, extent and masking to match dembuf 
# There remains some NA values for land aeas

#source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

# converts any pixels with value not equal to 0 in mask to NA in raster
cloud2NA <- function(x, y){
  #x[y>1] <- NA # leaves water pixels
  x[y != 0] <- NA # sets all flags to NA incl water pixels
  return(x)
}
# TO TRIM landsat 5 images of boundary 
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

# load stack of albedo rasters
#dir_input<-"~/Documents/Exeter/Data2015/Albedo/landsat/albedo/"
dir_input<-"/home/ISAD/jm622/Data2015/Albedo/landsat/albedo/"
filenames<-list.files(dir_input)
filenames<-paste(rep(dir_input,length(filenames)),filenames,sep="")
print(filenames)
# load cloud files
#dir_clouds<-"~/Documents/Exeter/Data2015/Albedo/landsat/cloudmask/"
dir_clouds<-"/home/ISAD/jm622/Data2015/Albedo/landsat/cloudmask/"
cloudnames<-list.files(dir_clouds)
cloudnames<-paste(rep(dir_clouds,length(cloudnames)),cloudnames,sep="")
print(cloudnames)


# Pre-selected files
filenames<-c("LE72030242003081SGS00-albedo.tif",
             "LE72030252002046SGS00-albedo.tif",
             "LE72030252005198EDC00-albedo.tif",
             "LE72030252007156EDC00-albedo.tif",
             "LE72040242001210EDC00-albedo.tif",
             "LE72040242006160EDC01-albedo.tif",
             "LE72040242009152ASN00-albedo.tif",
             "LE72040242011078ASN00-albedo.tif",
             "LE72040251999205AGS01-albedo.tif",
             "LE72040252005173ASN00-albedo.tif",
             "LE72040252007083ASN00-albedo.tif",
             "LE72040252011078ASN00-albedo.tif",
             "LE72050252015112NSG00-albedo.tif")
cloudnames<-c("LE72030242003081SGS00_cfmask.tif",
             "LE72030252002046SGS00_cfmask.tif",
              "LE72030252005198EDC00_cfmask.tif",
              "LE72030252007156EDC00_cfmask.tif",
              "LE72040242001210EDC00_cfmask.tif",
             "LE72040242006160EDC01_cfmask.tif",
              "LE72040242009152ASN00_cfmask.tif",
              "LE72040242011078ASN00_cfmask.tif",
             "LE72040251999205AGS01_cfmask.tif",
              "LE72040252005173ASN00_cfmask.tif",
             "LE72040252007083ASN00_cfmask.tif",
             "LE72040252011078ASN00_cfmask.tif",
              "LE72050252015112NSG00_cfmask.tif")

landsat<-substring(filenames,3,3)
filenames<-paste(rep(dir_input,length(filenames)),filenames,sep="")
print(filenames)
cloudnames<-paste(rep(dir_clouds,length(cloudnames)),cloudnames,sep="")
print(cloudnames)
print(landsat)

# Create blank raster of dem.buffer area with same projection and resolution as landsat images 
landsat.proj<-"+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
e.dembuf<-extent(269085,611835,5514015,5728575)  # calculated from extent(dembuf) to align with 30m resolution 
#blank.r<-projectExtent(e.dembuf,landsat.proj)
#res(blank.r)<-c(30,30)

albedo.map<-raster(ext=e.dembuf,res=c(30,30),crs=landsat.proj)
for (n in 1:length(filenames)){
  # Apply cloud mask
  print(filenames[n])
  albedo.r<-raster(filenames[n])
  #plot(albedo.r,main=filenames[n])
  print(cloudnames[n])
  cloudmask.r<-raster(cloudnames[n])
  #plot(cloudmask.r,main=cloudnames[n])
  # If landsat 5 trim image boundaries
  #if (landsat[n]==5) albedo.r<-trim.image(albedo.r)
  albedo.r<-trim.image(albedo.r)
  # Apply cloud mask
  new.r<-overlay(x=albedo.r, y=cloudmask.r, fun = cloud2NA)
  # Crop then extend to fit into area of interest e.dembuf
  new.r<-extend(crop(new.r,e.dembuf),e.dembuf) 
  albedo.map<-mosaic(albedo.map,new.r,fun=mean)
  plot(albedo.map,main=paste("Albedo map from ",n," images",sep=""))
  remove(new.r,cloudmask.r,albedo.r)
}
outfile<-paste(dir_albedo,"albedomap_mean.tif",sep="")
print(outfile)
writeRaster(albedo.map,file=outfile,format="GTiff",overwrite=TRUE)


# PLot and process map
brk<-c(0,0.5,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,1)
col<-rev(rainbow(14,start=1/6,end=4/6))
plot(albedo.map,col=col,main=paste("Albedo map from ",n," images",sep=""))

# Resample etc map to dembuf
albedo.map.os<-projectRaster(albedo.map,dembuf)
projection(albedo.map.os)<-"+init=epsg:27700"
plot(albedo.map.os)
albedo.dembuf<-resample(albedo.map.os,dembuf)
albedo.dembuf[is.na(albedo.dembuf)]<-0.2  # sets any NA values to 0.2 (will include sea - nbut these maksed later)
plot(albedo.dembuf)
albedo.dembuf<-mask(albedo.dembuf,dembuf)
projection(albedo.dembuf) = CRS("+init=epsg:27700")
plot(albedo.dembuf,col=col)  

outfile<-paste(dir_albedo,"albedo_mean_dembuf.tif",sep="")
print(outfile)
writeRaster(albedo.dembuf,file=outfile,format="GTiff",overwrite=TRUE)


writeRaster(adls48_31_21,filename="adls48_31_21",format="GTiff",dataType="INT2U",NAflag=0,overwrite=T)
