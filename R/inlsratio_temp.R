library(raster)
inv.dist<-function(x)
{
 d<-c(1:length(x))
 d<-1/d
 id<-d*x
 ido<-sum(id)
 ido<-ido/5.187378 # Source of value???
 ido
}

inv.ls<-function(landsea,direction)
{
 store<-array(0,dim=c(100,100,101))
 store[,,1]<-v[101:200,101:200]
 for (i in 1:100)
 {
  xshift<-round(i*sin(direction*(pi/180)),0)
  yshift<-round(i*cos(direction*(pi/180)),0)
  yshift<-yshift*(-1)
  store[,,(i+1)]<-v[(101+yshift):(200+yshift),(101+xshift):(200+xshift)]
 }
 storev<-array(store,dim=c(100*100,101))
 distance<-apply(storev,1,inv.dist)
 id<-array(distance,dim=c(100,100))
 id
}

# read in dem for SW as template
dem.template<-raster("C:/Jonathanmodel/wind/demsw.asc")
# read in dem for the whole of the UK
dem<-raster("C:/Jonathanmodel/wind/demoriginal.tif")
# select roughly area (matched to 10 km square)
e<-extent(c(65000,355000,0,170000))
dem<-crop(dem,e)
e<-extent(c(65000,355000,-10000,170000))
dem<-extend(dem,e)
for (dd in 0:35)
{
 direction<-dd*10
 # create a template raster for storing inverse land sea files
 et<-extent(c(75000,345000,10000,160000))
 inv.ls.ratio<-crop(dem,et)
 inv.ls.ratio<-inv.ls.ratio*0
 # NB dataset too big to work out inverse land -sea ratio in entirity, so does for each 10km square in turn (selected a
 for (b.col in 0:26){
 for (b.row in 0:15){
   
     xmn<-65000+(b.col*10000)
     xmx<-xmn+30000
     ymn<-(b.row*10000)-10000
     ymx<-ymn+30000
     eb<-extent(c(xmn,xmx,ymn,ymx))
   
     # test for single block
     xmn<-250000;xmx<-280000
     ymn<-30000; ymx<-60000
     eb<-extent(c(xmn,xmx,ymn,ymx))
     # end test code
   
     r.block<-crop(dem,eb)
     
     # convert to land = 1, sea = 0
     v<-getValues(r.block,format="matrix")
     sel1<-which(is.na(v)==T)
     sel2<-which(is.na(v)==F)
     v[sel1]<-0
     v[sel2]<-1
     
     id<-inv.ls(v,direction) # map of index - sea cells=0
     e.centre<-extent(xmn+10000,xmx-10000,ymn+10000,ymx-10000)
     r2<-raster(id,xmn=xmn+10000,xmx=xmx-10000,ymn=ymn+10000,ymx=ymx-10000)
     landsea.block<-crop(landsea.r,e.centre)
     
     #test plots
     final.r<-mask(r2, landsea.block, maskvalue=0)
     plot(crop(r.block,e.centre))
     plot(final.r)
     #end test
     
     inv.ls.ratio<-mosaic(inv.ls.ratio,r2,fun=max)
     plot(inv.ls.ratio)
 }}
 ecr<-extent(dem.template)
 inv.ls.ratio<-crop(inv.ls.ratio,ecr)
 # convert sea to NA
 crop.nas<-dem.template*0
 inv.ls.ratio<-inv.ls.ratio+crop.nas
 plot(inv.ls.ratio)
 fileout<-paste("C:/Jonathanmodel/landseartiodata/invratio_",direction,"deg.tif",sep="")
 writeRaster(inv.ls.ratio,file=fileout,overwrite=T)
}



