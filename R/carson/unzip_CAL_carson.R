##########################################################################################
# UNZIP INPUT DATA
# Unzips SSPressure and Temperature files 
##########################################################################################
# to be used via job array submission

args <-commandArgs(trailingOnly = TRUE)
root <- args[1] 
tarfile<-args[2]
print(paste("Root directory= ",root,sep=""))
print(paste("Infile= ",tarfile,sep=""))

# Dir
dir_caltar<-paste(root,"CMSAF-CAL/tar/",sep="")
dir_calncgz<-paste(root,"CMSAF-CAL/ncgz/",sep="")
dir_cal<-paste(root,"CMSAF-CAL/extract/",sep="")

library(R.utils) # IMPORTANT - hides raster functions requiring raster


# 1 Extract CMSAF tarfile
print("Extracting CAL tar files")
untar(tarfile=tarfile,exdir=dir_calncgz)

# 2 Extract ncgz files created from tar file (deletes ncgz file once extracted

gzfiles<-list.files(path=dir_calncgz,include.dirs=FALSE)
print(paste("Extracting ",length(gzfiles)," ncgz files to nc files...",sep=""))
  
for (n in 1:length(gzfiles)){
    gzfile<-paste(dir_calncgz,gzfiles[n],sep="")
    outfile<-paste(dir_cal,substr(gzfiles[n],1,34),sep="") 
    print(paste("Infile: ",gzfile,"  Outfile: ",outfile, sep=""))
    gunzip(filename=gzfile, destname=outfile, overwrite=TRUE, remove=TRUE)
}


