# Read in dem and change to extent to multiple of 10000
library(raster)
root<-"~/Documents/Exeter/Data2015/"

dir.dems<-paste(root,"basins/blocks/dems/",sep="")
dir.basins<-paste(root,"basins/blocks/basins/",sep="")
dir.strips<-paste(root,"basins/strips/",sep="")
dir.basinmap<-paste(root,"basins/",sep="")

# Load whole UK file and buffer
demuk<-raster(paste(root,"DEM/demoriginal.tif",sep=""))
projection(demuk)<-"+init=epsg:27700"
e.ukexp<-c(0,7e+05,-10000,1200000) # expand to allow 20km buffer to south of area of interest - set to sea (NA)
demuk<-extend(demuk,e.ukexp,values=NA)

e.dem<-extent(c( 120000,420000,10000,180000 )) # excludes scilly isles
#e.dem<-c(70000,350000,0,160000)
dem<-crop(demuk,e.dem)

# Chop DEM into 100 x 100 cell blocks
for (x in (xmin(e.dem)/10000):((xmax(e.dem)/10000)-1) )
{
  print(x)
  for (y in (ymin(e.dem)/10000):((ymax(e.dem)/10000)-1) )
  {
    xmn=x*10000
    ymn=y*10000
    e<-c(xmn,xmn+10000,ymn,ymn+10000)
    demc<-crop(dem,e)
    plot(demc) 
    v<-getValues(demc)
    sel<-which(is.na(v)==F)
    if (length(sel)>0)
    {
      rout<-paste(dir.dems,"dem_",x,"_",y,".tif",sep="")
      writeRaster(demc,file=rout,overwrite=T)
    }
  }
}
# Delineate basins
basin=1
counter<-0
fls<-list.files(dir.dems)    
#for (ii in 1:length(fls))
for (ii in 101:262)
{
   dem.b<-raster(paste(dir.dems,fls[ii],sep=""))    
   m1<-getValues(dem.b,format="matrix")
   # put buffer around dem
   sel<-which(is.na(m1)==T)
   m1[sel]<-(-999)
   m2<-array(-999,dim=c(102,102))
   m2[2:101,2:101]<-m1
   # create raster template
   e<-extent(dem.b)
   dem.b<-raster(m2,xmn=e@xmin-100,xmx=e@xmax+100,ymn=e@ymin-100,ymx=e@ymax+100)
   na.num<-which(m2==-999)
   # convert to data frame
   xyd<-data.frame(rasterToPoints(dem.b))
   sel<-which(xyd$layer==-999)
   xyd$layer[sel]<-NA
   xyd$x<-floor(xyd$x)
   xyd$y<-floor(xyd$y)
   xyd$hgt<-xyd$layer; xyd$layer<-NULL
   xyd$basin<-NA
   xyd$dn<-0
   test0<-0
   while (test0==0)
   {
      selt<-which(is.na(xyd$basin)==T)
      if (length(selt)==length(na.num)) test0<-1
      if (length(selt)>0)
      {
         sel<-which(xyd$hgt==min(xyd$hgt,na.rm=T) & is.na(xyd$basin)==T)
         xyd$basin[sel[1]]<-basin
         # select 8 surrounding cells
         test=0
         while(test==0)
         {
            tp1<-counter/500
            tp2<-floor(counter/500)
            if (tp1==tp2)
            {
               xyz<-data.frame(x=xyd$x,y=xyd$y,z=xyd$basin)
               r<-rasterFromXYZ(xyz)
               plot(r,main=basin)
            }
            sel2<-which(xyd$basin==basin & xyd$dn==0)
            if (length(sel2)==0) test<-1
            if (length(sel2)>0)
            {
               hgt<-xyd$hgt[sel2[1]]
               xs=xyd$x[sel2[1]]; ys=xyd$y[sel[1]]+100 #N
               xs[2]=xyd$x[sel2[1]]+100; ys[2]=xyd$y[sel2[1]]+100 #NE
               xs[3]=xyd$x[sel2[1]]+100; ys[3]=xyd$y[sel2[1]] #E
               xs[4]=xyd$x[sel2[1]]+100; ys[4]=xyd$y[sel2[1]]-100 #SE
               xs[5]=xyd$x[sel2[1]]; ys[5]=xyd$y[sel2[1]]-100 #S
               xs[6]=xyd$x[sel2[1]]-100; ys[6]=xyd$y[sel2[1]]-100 #SW
               xs[7]=xyd$x[sel2[1]]-100; ys[7]=xyd$y[sel2[1]] #W
               xs[8]=xyd$x[sel2[1]]-100; ys[8]=xyd$y[sel2[1]]+100 #NW
               for (i in 1:8)
               {
                  sel3<-which(xyd$x==xs[i] & xyd$y==ys[i])
                  if (is.na(xyd$hgt[sel3])==F)
                  {
                     if (xyd$hgt[sel3]>=hgt) xyd$basin[sel3]<-basin
                  }
               }
               xyd$dn[sel[1]]<-1
               xyd$dn[sel2[1]]<-1
               counter<-counter+1
            }
         }
      }
      sel5<-which(xyd$basin==basin)
      xyd$hgt[sel5]<-9999
      basin<-basin+1
   }
   xyz<-data.frame(x=xyd$x+0.5,y=xyd$y+0.5,z=xyd$basin)
   u<-unique(xyd$basin)
   tp<-paste("basins:",length(u),sep="")
   r<-rasterFromXYZ(xyz)
   r<-crop(r,e)
   plot(r,main=tp)
   fileout<-paste(dir.basins,fls[ii],sep="")
   writeRaster(r,file=fileout,overwrite=T)
}
# # # # # # # # # # # # # # # # # 
# merge blocks into strips  # # # # # # # # #
# # # # # # # # # # # # # # # # #  
basinjoin<-function(in1,in2)
{
m<-array(0,dim=c(length(in1),length(in1)))
for (iii in 1:length(in1))
{
  for (jjj in 1:length(in1))
  {
    x1<-c(in1[iii],in2[iii])
    x2<-c(in1[jjj],in2[jjj])
    x<-unique(c(x1,x2))
    m[iii,jjj]<-ifelse(length(x)<(length(x1)+length(x2)),1,0)
  }
}
for (iii in 1:length(in1))
{
  x<-m[iii,]
  sel<-which(x==1)
  for (jjj in 1:length(sel))
  {
     y<-m[sel[jjj],]
     sel2<-which(y==1)
     sel<-c(sel,sel2)
  }
  sel<-unique(sel)
  m[iii,sel]<-1
}
out1<-in1
out2<-in2
for (iii in 1:length(out1))
{
  sel<-which(m[,iii]==1)
  out1[sel]<-out1[iii]
  out2[sel]<-out1[iii]
}
dout<-data.frame(ins=c(in1,in2),outs=c(out1,out2))
dout<-unique(dout)
dout
}


library(raster)
fls<-list.files(dir.basins)
lns<-substr(fls,5,6)
lns<-as.numeric(gsub("_","",lns))
u<-unique(lns)
o<-order(u)
u<-u[o]
for (ii in 1:length(u))
{
sel<-which(lns==u[ii])
strp.fls<-fls[sel]
# order strips from lowest to highest
l<-nchar(strp.fls)
num<-substr(strp.fls,l-5,l-4)
num<-as.numeric(gsub("_","",num))
o<-order(num)
strp.fls<-strp.fls[o]
if (length(sel)==1)
{
  r.master<-raster(paste(dir.basins,strp.fls[1],sep=""))   
  rout<-paste(dir.strips,"strip_",u[ii],".tif",sep="")
  writeRaster(r.master,file=rout,overwrite=T)
}
if (length(sel)>1)
{
   r.master<-raster(paste(dir.basins,strp.fls[1],sep=""))   
   for (jj in 2:length(sel))
   {
      r.one<-raster(paste(dir.basins,strp.fls[jj],sep=""))   
      m1<-getValues(r.master,format="matrix")
      m2<-getValues(r.one,format="matrix")
      dfr<-data.frame(r.m=m1[1,],r.o=m2[100,])
      # remova NAs
      sel<-which(is.na(dfr$r.m)==F)
      dfr<-dfr[sel,]
      sel<-which(is.na(dfr$r.o)==F)
      dfr<-dfr[sel,]
      dfr<-unique(dfr)
      if (length(dfr$r.o)>0)
      {
        dout<-basinjoin(dfr$r.m,dfr$r.o)
        # combine rasters
        r.master<-mosaic(r.master,r.one,fun=mean)
        m<-getValues(r.master,format="matrix")
        for (bb in 1:length(dout$ins))
        {
          sel<-which(m==dout$ins[bb])
          m[sel]<-dout$outs[bb]
        }
        r.master<-raster(m,template=r.master)
      }
      if (length(dfr$r.o)==0) r.master<-mosaic(r.master,r.one,fun=mean)
      plot(r.master,main=u[ii])
   }
   rout<-paste(dir.strips,"strip_",u[ii],".tif",sep="")
   writeRaster(r.master,file=rout,overwrite=T)
}
}
   
# # # # # # # # # # # # # # # # # 
# merge strips  # # # # # # # # #
# # # # # # # # # # # # # # # # #  
library(raster)      
fls<-list.files(dir.strips)
# order strips from lowest to highest
l<-nchar(fls)
num<-substr(fls,l-5,l-4)
num<-as.numeric(gsub("_","",num))
o<-order(num)
fls<-fls[o]
r.master<-raster(paste(dir.strips,fls[1],sep=""))   
for (ii in 2:length(fls))
{
   r.one<-raster(paste(dir.strips,fls[ii],sep=""))      
   # set correct extent
   e1<-extent(r.master)
   e2<-extent(r.one)
   ymax<-max(e1@ymax,e2@ymax)
   ymin<-min(e1@ymin,e2@ymin)
   e3<-extent(c(e1@xmin,e1@xmax,ymin,ymax))
   e4<-extent(c(e2@xmin,e2@xmax,ymin,ymax))
   r.master<-extend(r.master,e3)
   r.one<-extend(r.one,e4)
   m1<-getValues(r.master,format="matrix")
   d1<-dim(m1)
   m2<-getValues(r.one,format="matrix")    
   dfr<-data.frame(r.m=m1[,d1[2]],r.o=m2[,1])
   # remova NAs
   sel<-which(is.na(dfr$r.m)==F)
   dfr<-dfr[sel,]
   sel<-which(is.na(dfr$r.o)==F)
   dfr<-dfr[sel,]
   dfr<-unique(dfr) 
   if (length(dfr$r.o)>0)
   {
        dout<-basinjoin(dfr$r.m,dfr$r.o)
        # combine rasters
        r.master<-mosaic(r.master,r.one,fun=mean)
        m<-getValues(r.master,format="matrix")
        for (bb in 1:length(dout$ins))
        {
          sel<-which(m==dout$ins[bb])
          m[sel]<-dout$outs[bb]
        }
        r.master<-raster(m,template=r.master)       
       
     }   
     if (length(dfr$r.o)==0) r.master<-mosaic(r.master,r.one,fun=mean)
     plot(r.master,main=ii)
}
# give basins sensible numbers
m<-getValues(r.master,format="matrix")
v<-getValues(r.master)
u<-unique(v)
sel<-which(is.na(u)==F)
u<-u[sel]
max(u)
for (ii in 1:length(u))
{
  tp1<-ii/500
  tp2<-floor(ii/500)
  sel<-which(m==u[ii])
  m[sel]<-ii+20000
  if (tp1==tp2)
  {
   r<-raster(m,template=r.master)
   plot(r,main=ii)
  }
}
r.master<-raster(m,template=r.master)

writeRaster(r.master,file=paste(dir.basinmap,"all.tif",sep=""),overwrite=T)


