# Calculate Ldif maps for different wind directions
# Input:  rasters of % land within radius at 100m 
#         rasters of inverse land-sea ratio (L) for dif wind directions at 100m
# Output: raster of Ldif (%land index - inv_ls index) for each wind direction
# Prev Programs: perc_land_maps, inlsratio_jm2


library(raster)
library(rgdal)

dir_percland<-"~/Documents/Exeter/Data2015/CoastEffect/percland/"
dir_grids<-"~/Documents/Exeter/Data2015/Templates/"
dir_lsratio<-"~/Documents/Exeter/Data2015/CoastEffect/lsratio/"
dir_ldif<-"~/Documents/Exeter/Data2015/CoastEffect/ldif/"
#dir_grids<-"C:/Data2015/Templates/"
#dir_percland<-"C:/Data2015/CoastEffect/percland/"

# Load % land maps
radius<-10000
in.file<-paste(dir_percland,"landin_",radius/1000,"km_complete_100m.r",sep="")
print(in.file)
pland100m.r<-raster(in.file)

#angle<-270

# For each wind direction and coastal effect map load inv_ls ratio maps
for (angle in seq(0,350,10))
{
# Load coast index maps and extract values for centre of 5km cells
in.file<-paste(dir_lsratio,"invratio_",direction,"deg.tif",sep="")
print(in.file)
inv.lsratio<-raster(in.file) # load raster of lsratio
# Calculate and save Ldif
Ldif.r<- pland100m.r-inv.lsratio
plot(Ldif.r)
out.file<-paste(dir_ldif,"ldif_from_pland_in_",radius/1000,"km.tif",sep="")
writeRaster(Ldif.r,file=out.file)

}
