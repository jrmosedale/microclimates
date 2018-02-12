library(R.utils)

extract.tar.to.ncdf<-function(dir_tar,dir_ncgz,dir_nc)
{
  #dir_tar<-"~/Documents/Exeter/Data2015/CMSAF-CAL/tar/"
  #dir_ncgz<-"/Users/jonathanmosedale/Documents/Exeter/Data2015/CMSAF-CAL/ncgz/"
  #dir_nc<-"~/Documents/Exeter/Data2015/CMSAF-CAL/extract/"  # location of fully extracted files to be used
  
  # Extract .tar files to .ncgz files
  tarfiles<-list.files(path=dir_tar,include.dirs=FALSE)
  print(tarfiles)
  
  for (n in 1:length(tarfiles)){
    tarfile<-paste(dir_tar,tarfiles[n],sep="")
    print(paste("Extracting... ",tarfile,sep=""))
    untar(tarfile=tarfile,exdir=dir_ncgz)
  }
  
  # Extract from ncgz files and delete original to save file space 
  
  gzfiles<-list.files(path=dir_ncgz,include.dirs=FALSE)
  print("Extracting ncgz files to nc files...")
  #print(gzfiles)
  
  for (n in 1:length(gzfiles)){
    gzfile<-paste(dir_ncgz,gzfiles[n],sep="")
    outfile<-paste(dir_nc,substr(gzfiles[n],1,34),sep="") 
    print(paste("Infile: ",gzfile,"  Outfile: ",outfile, sep=""))
    gunzip(filename=gzfile, destname=outfile, overwrite=TRUE, remove=TRUE)
  }
} # end function


# Extract CAL tar files to ncdf
extract.tar.to.ncdf(dir_tar,dir_ncgz,dir_nc)



#Extract, reproject and resample to 100m for Lizard
# Geog extent
e.lizard<-extent(160000,185000,10000,35000)

# Time period 
year<-  ;   month<-  ; day<-  ; 

for (hr in 0:23)
{
    mn<-0
    datetime<-paste(year,"/",sprintf("%02d",month,sep=""),"/",sprintf("%02d",day,sep="")," ",sprintf("%02d",hr,sep=""),":",sprintf("%02d",mn,sep=""),sep="")
    print(datetime)
    
    # Read in data from ncdf file
    infile.cal<-paste(dir_dni,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sprintf("%02d",hr,sep=""),sprintf("%02d",mn,sep=""),"002UD1000101UD.nc",sep="")
    ncdf_cal<-nc_open(infile.cal)
    cal.r<-raster(infile.cal)
    e<-extent(c(-7,-2,49,52))
    cal.r<-crop(cal.r,e)
    projection(cal.r)<-"+init=epsg:4326"
    
    # reproject to OSGB and set extent to e.lizard
    cal.r<-projectRaster(cal.r,crs="+init=epsg:27700")
    cal.r<-crop(cal.r,e.lizard)
    plot(cal.r,main=paste("CAL for Lizard area on ",day,"/",month,"/",year," at ",hr,":00",sep=""))
}
