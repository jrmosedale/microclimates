# Assign temperature values from 5km cells to 100m cells in block

# Function returns 100m raster from 5km raster for block
# Compare percentland_map_function & elevdif_map_function

tref.block<-function(dem.block,tref5km.r)
  {
  tref.block<-resample(tref5km.r,dem.block,method="ngb")
  tref.block<-mask(tref.block,dem.block)
  plot(tref.block,main="Tref block") # will generally be 
  return(tref.block)
} # end function 
  

