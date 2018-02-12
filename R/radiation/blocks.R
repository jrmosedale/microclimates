library(raster)
# read in DEM
r<-raster("C:/Jonathanmodel/wind/demsw.asc")
e<-extent(c(78000,346000,-1000,161000))
r<-extend(r,e)
store<-data.frame(block.x=0,block.y=0,val=0,
                  c.xmn=0,c.xmx=0,c.ymn=0,c.ymx=0,
                  b.xmn=0,b.xmx=0,b.ymn=0,b.ymx=0)
# extend dem so divisable by 1 km

i<-1
for (block.x in 79:345)
{
  for (block.y in 0:160)
  {
  tp<-paste("block.x = ",block.x," block.y= ",block.y,sep="")
  print(tp)
  # select one kilometre block
  c.xmn=block.x*1000
  c.xmx=c.xmn+1000
  c.ymn=block.y*1000
  c.ymx=c.ymn+1000
  # select for buffer
  b.xmn=c.xmn-1000
  b.xmx=c.xmx+1000
  b.ymn=c.ymn-1000
  b.ymx=c.ymx+1000
  e<-extent(c(c.xmn,c.xmx,c.ymn,c.ymx))
  crop.r<-crop(r,e)
  v<-getValues(crop.r)
  m<-mean(v,na.rm=T)
  store[i,1]<-block.x
  store[i,2]<-block.y
  store[i,3]<-m
  store[i,4]<-c.xmn
  store[i,5]<-c.xmx
  store[i,6]<-c.ymn
  store[i,7]<-c.ymx
  store[i,8]<-b.xmn
  store[i,9]<-b.xmx
  store[i,10]<-b.ymn
  store[i,11]<-b.ymx
  i<-i+1
  }
}
write.csv(store,file="C:/Jonathanmodel/blocks.csv",row.names=F)