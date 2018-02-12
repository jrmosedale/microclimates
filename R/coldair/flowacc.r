flowacc<-function(oneb)
{
   test=0
   counter=1
   while(test==0)
   {
      sel.test<-which(oneb$dun==0)
      oneb2<-oneb[sel.test,]
      # select highest pixel
      selm<-which(oneb2$hgt==max(oneb2$hgt))
      maxb<-oneb2[selm[1],]
      # set to done
      selm<-which(oneb$x==maxb$x & oneb$y==maxb$y)
      oneb$dun[selm[1]]<-1
      # get pixels in each of 8 compass directions
      xs=maxb$x; ys=maxb$y+100 #N
      xs[2]=maxb$x+100; ys[2]=maxb$y+100 #NE
      xs[3]=maxb$x+100; ys[3]=maxb$y #E
      xs[4]=maxb$x+100; ys[4]=maxb$y-100 #SE
      xs[5]=maxb$x; ys[5]=maxb$y-100 #S
      xs[6]=maxb$x-100; ys[6]=maxb$y-100 #SW
      xs[7]=maxb$x-100; ys[7]=maxb$y #W
      xs[8]=maxb$x-100; ys[8]=maxb$y+100 #NW
      surround<-data.frame(x=xs,y=ys,hgt=NA,flowacc=1)
      for (i in 1:8)
      {
         sel3<-which(oneb$x==xs[i] & oneb$y==ys[i])
         if(length(sel3)>0) surround$hgt[i]<-oneb$hgt[sel3]
      }
      ht<-which(is.na(surround$hgt)==F)
      if (length(ht)>0)
      {
         sel4<-which(surround$hgt==min(surround$hgt,na.rm=T))
         surround$flowacc[sel4]<-surround$flowacc[sel4]+maxb$flowacc/length(sel4)
         for (i in 1:length(sel4))
         {
            sel<-which(oneb$x==surround$x[sel4[i]] & oneb$y==surround$y[sel4[i]] & oneb$dun==0)
            oneb$flowacc[sel]<-surround$flowacc[sel4[i]]
         }
      }
      counter<-counter+1 ; print(counter)
      if (counter>1000000)
      {
        tp<-paste("too many iterations in basin ",u[xx],sep="")
        stop(tp)
      }
      # select only those not done
      sel.test<-which(oneb$dun==0)
      test<-ifelse(length(sel.test)>0,0,1)
    }
oneb
}

library(raster)
dir.flowstrips<-paste(root,"basins/flowacc/strips/",sep="")
dir.flowmap<-paste(root,"basins/flowacc/",sep="")


u<-u[o]
basins<-raster(paste(dir.basinmap,"all.tif",sep=""))
basins<-shift(basins,x=-0.5,y=-0.5)

dem<-raster("C:/Jonathanmodel/DEM/demsw.asc")
e<-extent(c(79400,350000,0, 160000))
basins<-extend(basins,e)
dem<-extend(dem,e)

# convert to data.frame
xyb<-data.frame(rasterToPoints(basins))
xyb$x<-floor(xyb$x/100)*100
xyb$y<-floor(xyb$y/100)*100
names(xyb)[3]<-"basin"
# convert to data.frame
xyd<-data.frame(rasterToPoints(dem))
xyd$x<-floor(xyd$x/100)*100
xyd$y<-floor(xyd$y/100)*100
names(xyd)[3]<-"hgt"
# merge
xydb<-merge(xyb,xyd,by=c("x","y"))
# unique basins
u<-unique(xydb$basin)
sel<-which(is.na(u)==F)
u<-u[sel]
o<-order(u)
for (block in 17:17)
{
  print(block) 
  xx<-1+block*1000
   # select each basin in turn
   sel<-which(xydb$basin==u[xx])
   oneb<-xydb[sel,]
   oneb$flowacc<-1
   oneb$dun<-0
   master<-flowacc(oneb)
   for (ii in 2:1000)
   {
     print(ii)
     xx<-ii+block*1000
     tp1<-ii/50
     tp2<-floor(ii/50)
     sel<-which(xydb$basin==u[xx])
     oneb<-xydb[sel,]
     oneb$flowacc<-1
     oneb$dun<-0
     fa<-flowacc(oneb)
     master<-rbind(master,fa)
     if (tp1==tp2)
     {
        xyz<-data.frame(x=master$x,y=master$y,z=master$flowacc)
        r<-rasterFromXYZ(xyz)
        plot(r,main=ii)
     }
   }
   xyz<-data.frame(x=master$x,y=master$y,z=master$flowacc)
   r<-rasterFromXYZ(xyz)
   plot(r,main=block)
   rout<-paste(dir.flowstrips,"strip_",block,".tif",sep="")
   writeRaster(r,file=rout,overwrite=T)
}

# tile together
library(raster)
block=0
r.master<-raster(paste(dir.flowstrips,"strip_",block,".tif",sep=""))
for (block in 1:17)
{
  r.one<-raster(paste(dir.flowstrips,"strip_",block,".tif",sep=""))
  r.master<-mosaic(r.master,r.one,fun=mean)
  plot(r.master,main=block)
}
writeRaster(r.master,file=paste(dir.flowmap,"all.tif",sep=""),overwrite=T)

# # # topgraphic wetness index
projection(dem)<-"+init=epsg:27700"
slope<-terrain(dem, opt='slope', unit='radians', neighbors=8)
tanb<-tan(slope)
e<-extent(dem)
r.master<-extend(r.master,e)
r.master<-shift(r.master,-50,-50)
topidx<-log(r.master/(tanb+0.05))
writeRaster(topidx,file=paste(dir.basinmap,"topidx.tif",overwrite=T))


## Basin altitudinal range
# i.e. Difference between highest and lowest point in basin
# #
basins<-raster(paste(dir.basinmap,"all.tif",sep=""))
basins<-shift(basins,x=-0.5,y=-0.5)

dem<-raster("C:/Jonathanmodel/DEM/demsw.asc")
e<-extent(c(79400,350000,0, 160000))
basins<-extend(basins,e)
dem<-extend(dem,e)

# convert to data.frame
xyb<-data.frame(rasterToPoints(basins))
xyb$x<-floor(xyb$x/100)*100
xyb$y<-floor(xyb$y/100)*100
names(xyb)[3]<-"basin"
# convert to data.frame
xyd<-data.frame(rasterToPoints(dem))
xyd$x<-floor(xyd$x/100)*100
xyd$y<-floor(xyd$y/100)*100
names(xyd)[3]<-"hgt"
# merge
xydb<-merge(xyb,xyd,by=c("x","y"))
xydb$altrange<-0
# unique basins
u<-unique(xydb$basin)
sel<-which(is.na(u)==F)
u<-u[sel]
o<-order(u)
u<-u[o]
for (ii in 1:length(u))
{
  tp1<-ii/50
  tp2<-floor(ii/50)
  # select each basin in turn
  sel<-which(xydb$basin==u[ii])
  oneb<-xydb[sel,]
  mn<-min(oneb$hgt,na.rm=T)
  mx<-max(oneb$hgt,na.rm=T)
  xydb$altrange[sel]<-(mx-mn)
  if (tp1==tp2)
  {
        xyz<-data.frame(x=xydb$x,y=xydb$y,z=xydb$altrange)
        r<-rasterFromXYZ(xyz)
        plot(r,main=ii)
  }
}
xyz<-data.frame(x=xydb$x,y=xydb$y,z=xydb$altrange)
r<-rasterFromXYZ(xyz)
writeRaster(r,file=paste(dir.basinmap,"altrange.tif",sep=""),overwrite=T)

## Basin altitudinal difference
# i.e. Difference between each grid cell and highest in basin
# #

library(raster)
basins<-raster(paste(dir.basinmap,"all.tif",sep=""))
basins<-shift(basins,x=-0.5,y=-0.5)
dem<-raster("C:/Jonathanmodel/DEM/demsw.asc")
e<-extent(c(79400,350000,0, 160000))
basins<-extend(basins,e)
dem<-extend(dem,e)
# convert to data.frame
xyb<-data.frame(rasterToPoints(basins))
xyb$x<-floor(xyb$x/100)*100
xyb$y<-floor(xyb$y/100)*100
names(xyb)[3]<-"basin"
# convert to data.frame
xyd<-data.frame(rasterToPoints(dem))
xyd$x<-floor(xyd$x/100)*100
xyd$y<-floor(xyd$y/100)*100
names(xyd)[3]<-"hgt"
# merge
xydb<-merge(xyb,xyd,by=c("x","y"))
# aggregate by basin to get max height
mx<-function(x)
{
   y<-max(x,na.rm=T)
   y
}
agg<-aggregate(xydb$hgt,by=list(xydb$basin),mx)
agg$basin<-agg$Group.1; agg$Group.1<-NULL
agg$maxhgt<-agg$x; agg$x<-NULL
xydb<-merge(xydb,agg,by=c("basin"),all=T)
xydb$dif.hgt<-xydb$maxhgt-xydb$hgt
xyz<-data.frame(x=xydb$x,y=xydb$y,z=xydb$dif.hgt)
r<-rasterFromXYZ(xyz)
plot(r)
writeRaster(r,file=paste(dir.basinmap,"altdif.tif",sep=""),overwrite=T)
## Basin altitudinal difference