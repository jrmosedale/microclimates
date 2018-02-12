# Assign temperature values from 5km cells to 100m cells in block

# Function returns 100m raster from 5km raster for block
# Compare percentland_map_function & elevdif_map_function

ref5km.to.block100m<-function(dem.block,ref5km.r)
  {
  ref.block<-resample(ref5km.r,dem.block,method="ngb")
  ref.block<-mask(ref.block,dem.block)
  plot(ref.block,main="ref block")
  return(ref.block)
} # end function 
  

