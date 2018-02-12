# Code to compare rasters 
# Raster groupings:
#1 dem elevdif albedo ldif prevtref tic
#2 altdif 
#3 twi coldair.block - missing all coastal cells

r1<-dem.block
r2<-crop(raster(paste(dir.basinmap,"flowacc-all.tif",sep="")),dem.block)
plot(r1)
plot(r2)
test1<-calc(r1,fun=function(x){ifelse(!is.na(x),1000,0)})
plot(test1)
test2<-calc(r2,fun=function(x){ifelse(!is.na(x),500,0)})
plot(test2)
test3<-test1+test2
plot(test3)
