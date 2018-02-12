
## TILE RASTER
## get filesnames (assuming the datasets were downloaded already. 
## please see http://thebiobucket.blogspot.co.at/2013/06/use-r-to-bulk-download-digital.html on how to download high-resolution DEMs)
setwd("D:/GIS_DataBase/DEM")
files <- dir(pattern = ".hgt")

files<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_100m.tif", sep="") 

## function for single file processing mind to replace the PATH to gdalinfo.exe!
## s = division applied to each side of raster, i.e. s = 2 gives 4 tiles, 3 gives 9, etc.
split_raster <- function(file, s = 2) {
  
  filename <- gsub(".hgt", "", file)
  gdalinfo_str <- paste0("\"C:/OSGeo4W64/bin/gdalinfo.exe\" ", file)
  
  # pick size of each side
  x <- as.numeric(gsub("[^0-9]", "", unlist(strsplit(system(gdalinfo_str, intern = T)[3], ", "))))[1]
  y <- as.numeric(gsub("[^0-9]", "", unlist(strsplit(system(gdalinfo_str, intern = T)[3], ", "))))[2]
  
  # t is nr. of iterations per side
  t <- s - 1
  for (i in 0:t) {
    for (j in 0:t) {
      # [-srcwin xoff yoff xsize ysize] src_dataset dst_dataset
      srcwin_str <- paste("-srcwin ", i * x/s, j * y/s, x/s, y/s)
      gdal_str <- paste0("\"C:/OSGeo4W64/bin/gdal_translate.exe\" ", srcwin_str, " ", "\"", file, "\" ", "\"", filename, "_", i, "_", j, ".tif\"")
      system(gdal_str)
    }
  }
}

## process all files and save to same directory
mapply(split_raster, files, 2) 