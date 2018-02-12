# Extract dembuf and dem for isles of scilly
e.sc<-extent(c(80000,100000,5000,20000 ))
dem.sc<-crop(demuk,e.sc)
plot(dem.sc)
demscilly.file<-paste(dir_dem,"demscilly.tif",sep="")
print(demscilly.file)
writeRaster(dem.sc,file=demscilly.file)

buffer<-20000 # 20km
e.scbuf<-extent(xmin(dem.sc)-buffer,xmax(dem.sc)+buffer,ymin(dem.sc)-buffer,ymax(dem.sc)+buffer)# Run setup programs for creating constant raster maps etc
dem.scbuf<-crop(demuk,e.scbuf)
plot(dem.scbuf)
dembufscilly.file<-paste(dir_dem,"dembufscilly.tif",sep="")
print(dembufscilly.file)
writeRaster(dem.scbuf,file=dembufscilly.file)
