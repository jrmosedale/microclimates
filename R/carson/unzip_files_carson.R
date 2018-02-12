##########################################################################################
# UNZIP INPUT DATA
# Unzips SSPressure and Temperature files 
##########################################################################################
# FUNCTIONS USED
##########################################################################################

# EXTRACT CMSAF (SIS, DNI, CAL) tar files to ncdf files via gz files 
# Leaves tar files but deletes gz files
extract.tar.to.ncdf<-function(dir_tar,dir_ncgz,dir_nc)
{
  #dir_tar<-"~/Documents/Exeter/Data2015/CMSAF-DNI/tar/"
  #dir_ncgz<-"/Users/jonathanmosedale/Documents/Exeter/Data2015/CMSAF-DNI/ncgz/"
  #dir_nc<-"~/Documents/Exeter/Data2015/CMSAF-DNI/extract/"  # location of fully extracted files to be used
  #library(R.utils) # IMPORTANT - hides raster functions requiring raster
  
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
  print(paste("Extracting ",length(gzfiles),"ncgz files to nc files...",sep=""))
  
  for (n in 1:length(gzfiles)){
    gzfile<-paste(dir_ncgz,gzfiles[n],sep="")
    outfile<-paste(dir_nc,substr(gzfiles[n],1,34),sep="") 
    print(paste("Infile: ",gzfile,"  Outfile: ",outfile, sep=""))
    gunzip(filename=gzfile, destname=outfile, overwrite=TRUE, remove=TRUE)
  }
} # end function



# Unzip yearly temperature files in dir_zip required for time period and save daily files in dir_temp
unzip_allfiles<-function(dir_zip,dir_temp)
{
  # Unzip year files of temperature data here
  zipfiles<-list.files(path=dir_zip,include.dirs=FALSE)
  print("Extracting .zip files...")

  for (n in 1:length(zipfiles)){
    zip.file<-paste(dir_zip,zipfiles[n], sep="")
    print (paste("Unzipping ",zip.file,sep=""))
    unzip(zip.file, exdir=dir_temp)
   }
} # end function

##########################################################################################
# CALL FUNCTIONS
##########################################################################################
args <-commandArgs(trailingOnly = TRUE)
root <- args[1] 
print(paste("Root directory= ",root,sep=""))

# Radiation budget
dir_dnitar<-paste(root,"CMSAF-DNI/tar/",sep="")
dir_dnigz<-paste(root,"CMSAF-DNI/ncgz/",sep="")
dir_sistar<-paste(root,"CMSAF-SIS/tar/",sep="")
dir_sisgz<-paste(root,"CMSAF-SIS/ncgz/",sep="")
dir_caltar<-paste(root,"CMSAF-CAL/tar/",sep="")
dir_calgz<-paste(root,"CMSAF-CAL/ncgz/",sep="")

dir_sis<-paste(root,"CMSAF-SIS/extract/",sep="")
dir_dni<-paste(root,"CMSAF-DNI/extract/",sep="")
dir_cal<-paste(root,"CMSAF-CAL/extract/",sep="")

# Sea Pressure
dir_pressure<-paste(root,"SeaPressure/",sep="")

library(R.utils) # IMPORTANT - hides raster functions requiring raster

# 1. Unzip data files 
# Various programs to be run only once (eg unzip, download, basic provessing of data sources)

# SEA SURFACE PRESSURE data - daily 0.25 deg - load and crop raster bands to area of interest
#gzfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc.gz",sep="")
#ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")
#print(paste("Extract ",gzfile," to ",ncfile,sep=""))
#gunzip(filename=gzfile, destname=ncfile, overwrite=TRUE)

# Extract CMSAF Radiation datafiles. Prog: extract_tar_ncdf_function
# WARNING - currently extracts ALL tar files - produces 8750 files per year per factor
print("Extracting radiation files")
extract.tar.to.ncdf(dir_dnitar,dir_dnigz,dir_dni)
extract.tar.to.ncdf(dir_sistar,dir_sisgz,dir_sis)
extract.tar.to.ncdf(dir_caltar,dir_calncgz,dir_cal)

# Unzip TEMPERATURE files required - start.jd -1?? 
#unzip_allfiles(dir_zip,dir_temp)

# WInd already extracted to r files


