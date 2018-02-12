library(raster)
library(rgdal)
ukgrid <- "+init=epsg:27700"

dir_tellus<-"C:/Data2015/Tellus_dtm/"
dir_sw<-paste(dir_tellus,"SW/",sep="")
dir_ss<-paste(dir_tellus,"SS/",sep="")
dir_sx<-paste(dir_tellus,"SX/",sep="")

# set here
dir_used<-dir_sx

# Using sw files/directories

# for each subdirectory in sw
subdir<-list.files(path=dir_used,include.dirs=TRUE)
print(subdir)

for (sd in 2:length(subdir.sw))
{
    sdfiles<-list.files(path=paste(dir_used,subdir[sd],sep=""))
    print(sdfiles)
    #alltiles.r<-raster(nrow=5000,ncol=5000,res=1,, crs=(ukgrid))
    #r.stack<-stack()
    
    for (i in 1:length(sdfiles))
    {
      infile<-paste(dir_used,subdir[sd],"/",sdfiles[i],sep="")
      print(paste("In file: ",infile,sep=""))
      if (i==1){alltiles.r<-raster(infile, crs=(ukgrid) )
      } else {tile.r<-raster(infile, crs=(ukgrid) )}
    
      if (i>1){alltiles.r<-merge(alltiles.r,tile.r)} 
    } # end for i in sdfiles[sd]
    
    #r<-raster(infile, crs=(ukgrid))
    plot(alltiles.r)
    
    outfile<-paste(dir_used,subdir[sd],"all.tif",sep="")
    print(outfile)
    writeRaster(alltiles.r,file=outfile)
    remove(alltiles.r)
}# end for sd in subdir
    

#all10m.r<-aggregate(alltiles.r,10)
#plot(all10m.r)

# Create arrays of xy extent for each tile and filename
num.tiles<-
ymn<-rep(NA,num.tiles)
ymx<-ymn; xmn<-ymn; xmx<-ymn
infile<-rep("",num.tiles)

for (dir_used in c()) {
  for (tile in list (files in dir)){
    infile[tile]<-paste()
    tile.r<-raster(file=infile[tile])
    xmn[tile]<-xmin(tile.r)
    xmx[tile]<-xmax(tile.r)
    ymn[tile]<-ymin(tile.r)
    ymx[tile]<-ymax(tile.r)
  }
}
tileinfo<-cbind(infile,xmn,xmx,ymn,ymx)
outfile<-paste()
save(tileinfo,file=outfile)

