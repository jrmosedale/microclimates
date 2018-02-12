# Calculates inv_land:sea ratio for different wind directions
# Adapted from Ilya's program
# Input:wind direction, sst for time t - sst must include 10kmbuffer around area of interest
# Output: rasters of ls ratio at 100m for each wind direction  

#####################################################################
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
 store[,,1]<-landsea[101:200,101:200]
 for (i in 1:100)
 {
  xshift<-round(i*sin(direction*(pi/180)),0)
  yshift<-round(i*cos(direction*(pi/180)),0)
  yshift<-yshift*(-1)
  store[,,(i+1)]<-landsea[(101+yshift):(200+yshift),(101+xshift):(200+xshift)]
 }
 storev<-array(store,dim=c(100*100,101))
 distance<-apply(storev,1,inv.dist)
 id<-array(distance,dim=c(100,100))
 id
}

#####################################################################
# Main Function called
# Uses min 10km buffer around area of interest
#####################################################################
inv.lsmaps<-function (dembuf,dem,interval,dir_lsratio){
  
    for (direction in seq(0,360,interval))
    {
     print(direction)  # wind direction (origin)
     # create a template raster for storing inverse land sea files of same dimensions as dem
     inv.lsratio<-raster(extent(dem),res=res(dem),nrows=nrow(dem),ncols=ncol(dem))
     inv.lsratio<-setValues(inv.lsratio,rep(0,ncell(inv.lsratio)))
     
     # NB dataset too big to work out inverse land -sea ratio in entirity, so does for each 10km square in turn 
     numcols<-( ((xmax(inv.lsratio)-xmin(inv.lsratio)) /10000 )-1) 
     numrows<-( ((ymax(inv.lsratio)-ymin(inv.lsratio)) /10000 )-1)
     for (b.col in 0:numcols ){
     for (b.row in 0:numrows ){
       # define 30km sq block
         xmn<-xmin(inv.lsratio)+(b.col*10000)-10000
         xmx<-xmn+30000
         ymn<-(b.row*10000)+ymin(inv.lsratio)-10000
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
         r2<-raster(id,xmn+10000,xmx-10000,ymn+10000,ymx-10000)
         #plot(r2)
         if ( !all(v==0)){ 
         inv.lsratio<-raster::mosaic(inv.lsratio,r2,fun=max) 
          } # end if       
      } # end b.row 
      #plot(inv.lsratio)
      print(paste("Col: ",b.col," Row ",b.row,sep=""))
      }# end b.col
    
     #crop to dem and convert sea to NA
     #plot(inv.lsratio)
     inv.lsratio<-mask(crop(inv.lsratio,dem),dem)
     plot(inv.lsratio)
    
    # Write lsratio raster file
     fileout<-paste(dir_lsratio,"invratio_",direction,"deg.tif",sep="")
     print(fileout)
     writeRaster(inv.lsratio,file=fileout,overwrite=T)
    } # end dd loop
    
} # end inv_lsmaps function



