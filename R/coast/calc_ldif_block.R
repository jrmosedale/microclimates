# Calculates single raster  ldif.block
# Input: ldif.stack of ldif for each wind direction
#        wind direction for block
#        direction interval for which ldif calculated

calc.ldif.block<-function(ldif.stack,wdir.block,interval=10){
  ldif.block<-raster
  wdir.layer<-round(wdir.block/interval)+1 # creates raster holding layer in ldif.stack for wind direction
  ldif.block<-stackSelect(ldif.stack,wdir.layer)
  return(ldif.block)
} # end function