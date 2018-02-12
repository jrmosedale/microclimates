# EXTRACT tar files to ncdf files via gz files 
# Leaves tar files but deletes gz files
extract.tar.to.ncdf<-function(dir_tar,dir_ncgz,dir_nc)
{
    #dir_tar<-"~/Documents/Exeter/Data2015/CMSAF-DNI/tar/"
    #dir_ncgz<-"/Users/jonathanmosedale/Documents/Exeter/Data2015/CMSAF-DNI/ncgz/"
    #dir_nc<-"~/Documents/Exeter/Data2015/CMSAF-DNI/extract/"  # location of fully extracted files to be used
    library(R.utils) # IMPORTANT - hides raster functions requiring raster
  
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




# Test output readable as ncdf file
#ncdf_test<-nc_open(outfile)
#test.r<-raster(outfile)
#e<-extent(c(-7,-2,49,52))
#r<-crop(test.r,e)
#plot(r)
