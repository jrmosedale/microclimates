# Prepare Sea level Pressure data - extract data at 100m resolution for time jd
# Requires JD functions
# Works for whole area if block=dembuf or similar or for small block
# Write raster for area grid

downscale.pressure<-function(block,ncfile,jd,write.file=FALSE)
{
  # unzip .gz file if necessary
  #gzfile<-paste(dir_pressure025,"pp_0.25deg_reg_v11.0.nc.gz",sep="")
  #ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")
  #gunzip(filename=gzfile, destname=ncfile, overwrite=TRUE)
  #ncdf_test<-nc_open(ncfile)
  
  # Select time bands 
  # Base jd value = 1/1/1950
  jd.base<-JDdmy(1,1,1950)
  jd.band<-jd-jd.base
  
  # Load raster from ncdf file 
  pressure.r<-raster(ncfile,band=jd.band)
  #plot(pressure.r)
  
  # Reproject cropped version to OSGB projection 
  pressure.r<-projectRaster(crop(pressure.r,c(-7,0,49,52)),crs="+init=epsg:27700")
  #plot(pressure.r)
  
  # Resample to 5km
  p5km.r<-resample(pressure.r,grid5km.r)
  p5km.r<-mask(p5km.r,grid5km.r)
  
  # The convert to 100m grid using nearest neighbour method 
  p.block<-ref5km.to.block100m(dem.block,p5km.r)
  
  # Resample & crop to 100m resolution of dembuf
  #p100m.r<-resample(pressure.r,block)
  #p100m.r<-mask(p100m.r,block)
  
  # Crop and write pressure raster (all times) for area of interest
  if(write.file){
  plot(p.block,main=paste("Sea level Pressure 100m cell for ",DMYjd(jd)$day,"/",DMYjd(jd)$month, "/",DMYjd(jd)$year,sep=""))
  outfile<-paste(dir_pressure,"pressure.tif",sep="")  
  writeRaster(file=outfile,p.block,overwrite=TRUE )
  }
  return(p.block)
} # end function

# Correct pressure for elevation according to hypsometric formula
# Inputs: sea level pressure, tmp at location, elevation above sea level
# References:  http://keisan.casio.com/exec/system/1224579725, http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric%20pressure
correct.pressure<-function(sl.pressure,tmp,elevation)
{
  kelvins<-273+tmp
  pressure<-sl.pressure*((kelvins-(0.0065*elevation)) / kelvins )^5.25588
  return(pressure)
}