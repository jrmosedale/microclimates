
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


args <-commandArgs(trailingOnly = TRUE)
dir_zip <- args[1] 
dir_temp<-args[2]
print(dir_zip)
print(dir_temp)

# Unzip TEMPERATURE files required - start.jd -1?? 
unzip_allfiles(dir_zip,dir_temp)

# qsub -v dirzip dirout