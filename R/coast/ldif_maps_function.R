# Calculate Ldif maps for different wind directions
# Input:  Lref rasters of % land of 5km referencecells within radius at resolution of 100m
#         rasters of inverse land-sea ratio (L) for dif wind directions at 100m
# Output: raster of Ldif (%land index - inv_ls index) for each wind direction
# Prev Programs: perc_land_maps_function, inv_lsratio_maps_function - ASSUMES maps of same dimensions and resolution

ldif.maps<-function(interval,radius,dir_percland,dir_lsratio,dir_ldif){
    # Load % land maps
    #radius<-10000
    in.file<-paste(dir_percland,"percent_land_",radius/1000,"km_lref_100mgrid.tif",sep="")
    print(in.file)
    pland100m.r<-raster(in.file)
    
    # For each wind direction and coastal effect map load inv_ls ratio maps
    for (direction in seq(0,350,interval))
    {
    # Load coast index maps and extract values for centre of 5km cells
    in.file<-paste(dir_lsratio,"invratio_",direction,"deg.tif",sep="")
    print(in.file)
    inv.lsratio<-raster(in.file) # load raster of lsratio
    # Calculate and save Ldif
    Ldif.r<- pland100m.r-inv.lsratio
    plot(Ldif.r,main=paste("Ldif for wind direction ",direction,sep=""))
    out.file<-paste(dir_ldif,"ldif_",direction,"deg_from_percland_in_",radius/1000,"km.tif",sep="")
    print(out.file)
    writeRaster(Ldif.r,file=out.file,overwrite=TRUE)
    }# end for
} # end function