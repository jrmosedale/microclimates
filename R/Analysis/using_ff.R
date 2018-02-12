

library(ff)

# create
fileout<-paste(dir_finalt,"block20km-",sprintf("%03d",ukcpcell,sep=""),"-",year,sep="")
print(fileout)
mat <- ff(vmode="single",dim=c(200,200,max.hr),filename=fileout)


print(max.hr)

# load one hr of data
test<-raster(mat[,,17],template=dem.block)


ffsave(mat,file=fileout)

close(mat)


ffload(file=fileout)

# extract
open.ff()
t<-x[hr]




ID_Raster <- raster(STACK[[1]])
ID_Raster[] <- 1:ncell(STACK[[1]])

Now I can use the extract function on this raster to identify the correct cell and the extract the corresponding values from the ff matrix, with the following lines:
  
  ext_ID <- extract(ID_Raster,MapUTM)
ext2 <- mat[as.numeric(ext_ID),]

