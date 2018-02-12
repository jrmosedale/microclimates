library(raster)
dir.in<-"C:/Jonathanmodel/coldairdrainage/datain/"
dir.out<-"C:/Jonathanmodel/coldairdrainage/dataout/"
# Read in basin data
library(raster)
dem<-raster(paste(dir.in,"demsw.asc",sep=""))
basins.r<-raster(paste(dir.out,"basinmap-all.tif",sep=""))


m<-getValues(dem,format="matrix")
flow.dir<-ifelse(is.na(m),NA,0)
# put buffer around m
m2<-array(-9999,dim=c(dim(dem)[1]+2,dim(dem)[2]+2)) # set boundary to very low value - to ensure flow direction of edge cells is out
m2[2:(dim(dem)[1]+1),2:(dim(dem)[2]+1)]<-m
m<-m2
done.m<-m # will be assigned as NA when analysed
bm<-getValues(basins.r,format="matrix")
bm2<-array(NA,dim=c(dim(bm)[1]+2,dim(bm)[2]+2))
bm2[2:(dim(bm)[1]+1),2:(dim(bm)[2]+1)]<-bm
u<-unique(bm)
# Creates flow direction matrix
  # 1 = NE
  # 2 = E
  # 3 = SE
  # 4 = N
  # 5 = 0
  # 6 = S
  # 7 = NW
  # 8 = W
  # 9 = SW
for (b in 1:max(u,na.rm=T))
{
   m3<-array(NA,dim=c(dim(m)))
   sel<-which(bm2==b)
   m3[sel]<-m[sel]
   sel2<-sel-1
   y<-sel2%/%dim(m)[1]+1
   x<-sel2%%dim(m)[1]+1
   for (ii in 1:length(y))
   {
      focalcell<-m3[x[ii],y[ii]]
       if (is.na(focalcell)==FALSE)
       {
         done.m[x[ii],y[ii]]<-NA
         m9<-matrix(c(m3[x[ii]-1,y[ii]-1],m3[x[ii],y[ii]-1],m3[x[ii]+1,y[ii]-1],m3[x[ii]-1,y[ii]],focalcell,
                      m3[x[ii]+1,y[ii]],m3[x[ii]-1,y[ii]+1],m3[x[ii],y[ii]+1],m3[x[ii]+1,y[ii]+1]),nrow=3)
         flow.dir[(x[ii]-1),(y[ii]-1)]<-which(m9==min(m9,na.rm=T))[1]
       }  # end if
   }  # end ii
   if (b%%1000==0 | b==max(u,na.rm=T)) # Plot and save every 1000 basins or if last basin
   {
      plot(raster(flow.dir,template=dem),main="Flow direction")
      print(paste("basin: ",b," of ",max(u,na.rm=T)," analysed",sep=""))
      fileout<-paste(dir.out,"flowdir.tif",sep="")
      flowdir.r<-raster(flow.dir,template=dem)
      writeRaster(flowdir.r,file=fileout,overwrite=TRUE)
   } # end if
}

#####################################
#Calculate flow accumulation
#####################################
library(raster)
dir.in<-"C:/Jonathanmodel/coldairdrainage/datain/"
dir.out<-"C:/Jonathanmodel/coldairdrainage/dataout/"
dem<-raster(paste(dir.in,"demsw.asc",sep=""))
basins.r<-raster(paste(dir.out,"basinmap-all.tif",sep=""))
flowdir.r<-raster(paste(dir.out,"flowdir.tif",sep=""))

dm<-getValues(dem,format="matrix")
bm<-getValues(basins.r,format="matrix")
fd<-getValues(flowdir.r,format="matrix")
# but buffer around flow direction
fa<-ifelse(is.na(dm),NA,1)
fa2<-array(NA,dim=c(dim(fa)[1]+2,dim(fa)[2]+2)) # set boundary to very high value - to ensure basins do not extend beyond area of interest
fa2[2:(dim(fa)[1]+1),2:(dim(fa)[2]+1)]<-fa
fa<-fa2
u<-unique(bm)
for (b in 1:max(u,na.rm=T))
{
   # select dem pixels for basin b
   sel<-which(bm==b)
   d<-dm[sel]
   # order  dem for selected basin from highest to lowest
   o<-order(d,decreasing = T)
   for (i in 1:length(o))
   {
      # get row and column of nth pixel of basin
      sel2<-sel[o[i]]-1
      cl<-sel2%/%dim(dm)[1]+1
      rw<-sel2%%dim(dm)[1]+1
      # get flow direction of focal cell
      fdr<-fd[rw,cl]
      if(is.na(fdr)==F)
      {
         if (fdr==1) {rw2<-rw-1; cl2<-cl-1}
         if (fdr==2) {rw2<-rw;   cl2<-cl-1}
         if (fdr==3) {rw2<-rw+1; cl2<-cl-1}
         if (fdr==4) {rw2<-rw-1; cl2<-cl}
         if (fdr==4) {rw2<-rw;   cl2<-cl}
         if (fdr==6) {rw2<-rw+1; cl2<-cl}
         if (fdr==7) {rw2<-rw-1; cl2<-cl+1}
         if (fdr==8) {rw2<-rw;   cl2<-cl+1}
         if (fdr==9) {rw2<-rw+1; cl2<-cl+1}
         fa[rw2+1,cl2+1]<-fa[rw2+1,cl2+1]+1
      }
   }
   if (b%%1000==0)  # Plot and save every 1000 basins
   {
      fac<-fa[2:(dim(fa)[1]-1),2:(dim(fa)[2]-1)]
      plot(raster(fac,template=dem),main=paste("basin: ",b,sep=""))
      fileout<-paste(dir.out,"flowacc.tif",sep="")
      flowacc<-raster(fac,template=dem)
      writeRaster(flowacc,file=fileout,overwrite=TRUE)

   }
}
fac<-fa[2:(dim(fa)[1]-1),2:(dim(fa)[2]-1)]
plot(raster(fac,template=dem),main=paste("basin: ",b,sep=""))
fileout<-paste(dir.out,"flowacc.tif",sep="")
flowacc<-raster(fac,template=dem)
writeRaster(flowacc,file=fileout,overwrite=TRUE)

# Log topographic wetness index
library(raster)
dir.in<-"C:/Jonathanmodel/coldairdrainage/datain/"
dir.out<-"C:/Jonathanmodel/coldairdrainage/dataout/"
dem<-raster(paste(dir.in,"demsw.asc",sep=""))
crs(dem)<-"+init=epsg:27700"
flowacc<-raster(paste(dir.out,"flowacc.tif",sep=""))

slope<-terrain(dem,opt='slope', unit='radians')
# set slope of zero to 0.001
sm<-getValues(slope,format="matrix")
sel<-which(sm<0.00001); sm[sel]<-0.00001
fm<-getValues(flowacc,format="matrix")
tm<-fm/tan(sm)
topidx<-raster(log(tm),template=dem)
plot(topidx)
fileout<-paste(dir.out,"log_topidx.tif",sep="")
writeRaster(topidx,file=fileout,overwrite=T)




