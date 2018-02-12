library(raster)
# Basin metrics:
 # Maximum elevation of basin
 # Maximum elevation of basin - pixel elevation (as a measure of how cold could air could be)
 # Elevation of lowest boundary basin - pixel elevation (as a measure of sink potential)
dir_in<-"C:/Jonathanmodel/coldairdrainage/datain/"
dir_out<-"C:/Jonathanmodel/coldairdrainage/dataout/"
dem<-raster(paste(dir_in,"demsw.asc",sep=""))
basins.r<-raster(paste(dir_out,"basinmap-all.tif",sep=""))

dm<-getValues(dem,format="matrix")
bm<-getValues(basins.r,format="matrix")
m1<-ifelse(is.na(dm),NA,0)
m2<-ifelse(is.na(dm),NA,0)
m3<-ifelse(is.na(dm),NA,0)
# set elevation of sea to zero
dsea<-dm
dm<-ifelse(is.na(dm),0,dm)

# put one metre buffer around bm
bm2<-array(-1,dim=c(dim(bm)[1]+2,dim(bm)[2]+2)) # set boundary to very high value - to ensure basins do not extend beyond area of interest
bm2[2:(dim(bm)[1]+1),2:(dim(bm)[2]+1)]<-bm
for (b in 1:max(bm,na.rm=TRUE))
{
   print(paste("b=",b))
   sel<-which(bm==b)
   # max elevation of basin
   m1[sel]<-max(dm[sel],na.rm=T)
   # max elevation of basin - pixel elevation
   m2[sel]<-max(dm[sel],na.rm=T)-dm[sel]
   bdry<-0
   for (i in 1:length(sel))
   {
      print(paste("i=",i))
      # get surrounding pixels form focal pixel
      sel2<-sel[[i]]-1
      cl<-sel2%/%dim(bm)[1]+2
      rw<-sel2%%dim(bm)[1]+2
      m9<-matrix(c(bm2[rw-1,cl-1],bm2[rw,cl-1],bm2[rw+1,cl-1],bm2[rw-1,cl],
                 bm2[rw,cl],bm2[rw+1,cl],bm2[rw-1,cl+1],bm2[rw,cl+1],bm2[rw+1,cl+1]),nrow=3)

      # find out whether boundary or not
      m92<-ifelse(m9[2,2]==m9,1,0) # detects if boundary to another basin
      #print(m92)
      bdry[i]<-min(m92,na.rm=T); print(paste("bdry=",bdry)) 
      
      # ADDITION get heights with sea marked NA and corrects if any sea present
      #dem9<-matrix(c(dsea[rw-2,cl-2],dsea[rw-1,cl-2],dsea[rw,cl-2],dsea[rw-2,cl-1],dsea[rw-1,cl-1],dsea[rw,cl-1],dsea[rw-2,cl],dsea[rw-1,cl],dsea[rw,cl]),nrow=3)
     # if (TRUE %in% is.na(dem9)) bdry[i]<-0  
      
   }
   # select heights of boundaries
   sel3<-which(bdry==0)
   ds<-dm[sel[sel3]]; print(ds)
   for (i in 1:length(sel)) m3[sel[i]]<-min(ds,na.rm=T)-dm[sel[i]]
   
   # PLot and write rasters
   if (b%%100==0 | b==max(bm,na.rm=TRUE))
   {
      r11<-raster(m1,template=dem)
      r12<-raster(m2,template=dem)
      r13<-raster(m3,template=dem)
      par(mfrow=c(2,2))
      plot(r11,main=paste("basin: ",b," max elev",sep=""))
      plot(r12,main=paste("basin: ",b," elev dif",sep=""))
      plot(r13,main=paste("basin: ",b," bound dif",sep=""))
      #writeRaster(r1,file=paste(dir_out,"maxbasinelevation.tif",sep=""),overwrite=T)
      #writeRaster(r2,file=paste(dir_out,"elevationdif.tif",sep=""),overwrite=T)
      #writeRaster(r3,file=paste(dir_out,"boundarydif.tif",sep=""),overwrite=T)
   }
}











