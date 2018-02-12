##########################################################################################
# UNZIP INPUT DATA
# Unzips SSPressure and Temperature files 
##########################################################################################
# FUNCTIONS USED
##########################################################################################

# EXTRACT landsat tar files to ncdf files via gz files 
# Leaves tar files but deletes gz files
extract.tar<-function(dir_tar,dir_nc)
{
  #dir_tar<-"~/Documents/Exeter/Data2015/CMSAF-DNI/tar/"
  #dir_nc<-"~/Documents/Exeter/Data2015/CMSAF-DNI/extract/"  # location of fully extracted files to be used
  #library(R.utils) # IMPORTANT - hides raster functions requiring raster
  
  # Extract .tar files to .ncgz files
  tarfiles<-list.files(path=dir_tar,include.dirs=FALSE)
  print(tarfiles)  
  for (n in 1:length(tarfiles)){
    tarfile<-paste(dir_tar,tarfiles[n],sep="")
    print(paste("Extracting... ",tarfile,sep=""))
    untar(tarfile=tarfile,exdir=dir_nc)
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
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

args <-commandArgs(trailingOnly = TRUE)
root <- args[1] 
print(paste("Root directory= ",root,sep=""))

library(R.utils) # IMPORTANT - hides raster functions requiring raster

# Extract CMSAF Radiation datafiles. Prog: extract_tar_ncdf_function
# WARNING - currently extracts ALL tar files - produces 8750 files per year per factor
print("Extracting landsat files")
#extract.tar(dir_landsattar,dir_landsatnc)


