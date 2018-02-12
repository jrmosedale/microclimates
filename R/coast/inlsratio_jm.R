library(raster)
dir_lsratio<-"~/Documents/Exeter/Data2015/CoastEffect/lsratio/"
  
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


# read in dem for the whole of the UK
demuk<-rasterdemuk<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif")
# Define sw dem of interest (no buffer)
#e.dem<-extent(c(70000,420000,0,180000 )) # includes scilly isles
e.dem<-extent(c( 130000,400000,10000,180000 )) # excludes scilly isles
dem<-crop(demuk,e.dem)
plot(dem)

# Create buffered south west dem raster required - USED to create coastal effect index maps
buffer<-10000 # include 10km buffer area around zone of interest and 10km of NA
e.buf<-extent(c( (xmin(dem)-buffer),(xmax(dem)+buffer),(ymin(dem)-buffer),(ymax(dem)+buffer) )) 
dembuf<-crop(demuk,e.buf)
# ensure divisible by 30km
xmax30<-(ceiling((xmax(e.buf)-xmin(e.buf))/30000)*30000)+xmin(dembuf)
ymax30<-(ceiling((ymax(e.buf)-ymin(e.buf))/30000)*30000)+ymin(dembuf)
e.buf<-extent(c(xmin(dembuf),xmax30,ymin(dembuf),ymax30)) 
dembuf<-extend(dembuf,e.buf)
plot(dembuf,main="DEM-buffered")

for (dd in 0:35)
{
 begin<-proc.time()
 direction<-dd*10
 # create a template raster for storing inverse land sea files of same dimensions as dem
 inv.lsratio<-raster(extent(dembuf),res=res(dembuf),nrows=nrow(dembuf),ncols=ncol(dembuf))

 # NB dataset too big to work out inverse land -sea ratio in entirity, so does for each 10km square in turn (selected a
 for (b.col in 0:( ((xmax(inv.lsratio)-xmin(inv.lsratio)) /10000 )-2)  ){
 for (b.row in 0:( ((ymax(inv.lsratio)-ymin(inv.lsratio)) /10000 )-2)  ){
   
   # define 30km sq block
     xmn<-xmin(inv.lsratio)+(b.col*10000)
     xmx<-xmn+30000
     ymn<-(b.row*10000)+ymin(inv.lsratio)
     ymx<-ymn+30000
     e.block<-extent(c(xmn,xmx,ymn,ymx))
     
     # test for single block
     #xmn<-250000;xmx<-280000
     #ymn<-30000; ymx<-60000
     #e.block<-extent(c(xmn,xmx,ymn,ymx))
     # end test code
   
     r.block<-crop(dembuf,e.block)
     
     # convert to land = 1, sea = 0
     v<-getValues(r.block,format="matrix")
     v<-ifelse(is.na(v),0,1)
     
     id<-inv.ls(v,direction) # map of index - sea cells=0
     r2<-raster(id,xmn=xmn+10000,xmx=xmx-10000,ymn=ymn+10000,ymx=ymx-10000)
     #plot(r2)
     
     inv.lsratio<-mosaic(inv.lsratio,r2,fun=max)
     #plot(inv.lsratio)
 }}

# crop to dem and convert sea to NA
 e.dem<-extent(dem)
 inv.lsratio<-crop(inv.lsratio,e.dem)
 inv.lsratio<-mask(inv.lsratio,dem)
 plot(inv.lsratio)

# Write lsratio file
 fileout<-paste(dir_lsratio,"invratio_",direction,"deg.tif",sep="")
 writeRaster(inv.lsratio,file=fileout,overwrite=T)

 end<-proc.time()-begin
 print(end)
}






# CUT OUTS and TESTS
# COMPARISON 
plot(inv.lsratio-coast.r)

xmn<-250000;xmx<-280000
ymn<-30000; ymx<-60000
e.block<-extent(c(xmn,xmx,ymn,ymx))

plot(crop(coast.r,e.block))
plot(crop(inv.lsratio,e.block))
plot(crop(inv.lsratio,e.block)-crop(coast.r,e.block))