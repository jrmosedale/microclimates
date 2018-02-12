dir_data<-paste(root,"proxyt100/",sep="")
dir_results<-paste(root,"proxyt100/riskmaps/",sep="")

r.s<-brick(paste(dir_data,"devon-tmax-1989.tif",sep=""))
plot(raster(r.s,layer=183))

r<-raster(paste(dir_results,"gdd10_gs_2009_cornwall.tif",sep=""))
r<-raster(paste(dir_results,"lastspfr_doy_2009_cornwall.tif",sep=""))
r<-raster(paste(dir_results,"fl_tmean_2009_cornwall.tif",sep=""))
r<-raster(paste(dir_results,"tmin_year_2009_cornwall.tif",sep=""))
r<-raster(paste(dir_results,"t25_gsdays_2009_cornwall.tif",sep=""))


plot(r)
plot(crop(r,extent(188000,207000,65000,80000)))

plot(crop(r,extent(-560000,-520000,6520000,6560000)))


rad<-deg*(pi/180)



filein<-paste(dir_data,"cornwall-tmax-2011.tif",sep="")
r.b<-brick(filein)
r.b
plot(raster(r.b,layer=50))
