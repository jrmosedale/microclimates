

rad_downscale_100m<-function(dem.buffer,dnr.24h,sis24h)
{
xdim<-(xmax(dem.buffer)-xmin(dem.buffer))
ydim<-(ymax(dem.buffer)-ymin(dem.buffer))
e.tps<-extent(xmin(dem.buffer)-xdim,xmax(dem.buffer)+xdim,ymin(dem.buffer)-ydim,ymax(dem.buffer)+ydim)# Run setup programs for creating constant raster maps etc

sis.stack<-stack()
dnr.stack<-stack()
for (lyr in 1:24){
  dnr.r<-raster(dnr.24h,layer=lyr); #plot(dnr.r)
  sis.r<-raster(sis.24h,layer=lyr); #plot(sis.r)
  
  sis.r<-crop(sis.r,e.tps)
  sis.hr<-tps.resample(sis.r,dem.buffer,FALSE)
  
  dnr.r<-crop(dnr.r,e.tps)
  dnr.hr<-tps.resample(dnr.r,dem.buffer,FALSE)

  sis.stack<-stack(sis.stack,sis.hr)
  dnr.stack<-stack(dnr.stack,dnr.hr)
}
rad.24h<-c(dnr.stack,sis.stack)
return(rad.24h)
} # end function