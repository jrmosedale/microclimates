# Calculate SST for each 5km cell given wind direction and sst map
# Input:  5km grid land/sea - core and buffered regions
#         raster of sst for t
#         wind direction for t
# Output: raster of SST-Tref at t

#####################################################################
# Buffered and normal 5km landsea grid = grid5kmbuf.r, grid5km.r
# Buffered and normal 100m  dem = dembuf, dem
#####################################################################

findsst<-function(x)
{
  is.sea<-ifelse(x==0,NA,1)
  nearest.sea<-match(1,is.sea)
  if (is.na(nearest.sea)){
    seacell<-0 
  } else {seacell<-x[nearest.sea]}  
  return(seacell)
}

# Finds nearest upwind sea cell within 'buffer' (20km) for cells 
# If no sea cell found returns value of...0 sea cells returned = NA
# Otherwise records upwind cell number in gridmask.r
# Input:  sst.r - sea surface temperatures 
#         grid5kmbuf.r (including buffer)
#         grid5km.r - defines ara of interest within buffer
#         direction assumed to be constant across area
#         distance = max distance a search for nearest sea cell (= buffer of 20km)

upwind.sst<-function(sst.r,grid5kmbuf.r,grid5km.r,direction,distance=20000)
{
  x<-dim(sst.r)[1]
  y<-dim(sst.r)[2]
  step<-distance/res(sst.r)[1] # max number of cells from focal cell to be searched
  
  # create matrix holding cell numbers for sea cells and 0 for landcells
  vals<-getValues(grid5kmbuf.r,format="matrix") # land/sea 1/NA values
  cells<-getValues(sst.r,format="matrix")
  seacells<-ifelse(is.na(vals),cells,0) # matrix of -999 if land or sea temperature if sea - could be changed to sst values as long as land =0
  #plot(raster(seacells,template=grid5kmbuf.r))
  
  store<-array(0,dim=c(dim(seacells)[1]-(2*step),dim(seacells)[2]-(2*step),step+1)) # 3d array to hold values for non-buffered region for every 'step'
  store[,,1]<-seacells[(step+1):(dim(seacells)[1]-step),(step+1):(dim(seacells)[2]-step)]
  #plot(raster(store[,,1],template=grid5km.r))
  
  for (i in 1:step)
  {
    xshift<-round(i*sin(direction*(pi/180)),0) 
    yshift<-round(i*cos(direction*(pi/180)),0)
    print(paste("i: ",i, " xshift: ",xshift," yshift: ",yshift,sep=""))
    yshift<-yshift*(-1)
    store[,,(i+1)]<-seacells[(step+1+yshift):(dim(seacells)[1]-step+yshift),(step+1+xshift):(dim(seacells)[2]-step+xshift)]
    
  } # end
  
  # use first/last to find nearest sea cell??
  storev<-array(store,dim=c((dim(store)[1]*dim(store)[2]),step+1))
  
  upwind.sst<-apply(storev,1,findsst) # if no sea within buffer then = NA?
  upwind.sst<-matrix(upwind.sea,nrow=nrow(store),ncol=ncol(store)) 
  
  upwind.sst.r<-raster(upwind.sst,template=grid5km.r)
  upwind.sst.r<-mask(upwind.sst,grid5km.r)
  plot(upwind.sst.r,main=paste("nearest seacell where direction= ",direction,sep=""))
  
  return(upwind.sst.r)
}# end function

#####################################################################
# Function to create maps recording nearest sea cell for different wind directions
# could create 5km raster block 1-36 layers according to wind dir
# record cell number of nearest sea cell!!
#PROBLEM - check cell values point to correct cell in matrix/raster
infile.sst<- paste(dir_sst,"hourly/sst_2010_2_3_7.tif",sep="")
infile.sst<-"~/Documents/Exeter/Data2015/sst/hourly/sst_2010_2_2_12h.tif"
sst.r<-raster(infile.sst)
#compareRaster(sst.r,grid5kmbuf.r)

write.upwind.maps<-function(dem.buffer,dem.block,dir_upwindsea)
{
  for (direction in seq(0,350,10)){
    #print(direction)
    #plot(dem.buffer)
    upwindsea.r<-upwind.sst(dem.buffer,dem.block,direction)
    # upwindsea100.r<-upwind.sea(gbuf100.r,direction,gmask100.r,buffer)
    out.file<-paste(dir_upwindsea,"upwind_5kmcell_in_dir_",direction,"_dist_",distance,".tif",sep="")
    print(out.file)
    writeRaster(upwindsea.r,file=out.file,overwrite=TRUE)
  }
}# end function

#####################################################################
# Function to load upwind maps into 3d matrix for block cells
block.upwind.sst<-function(dem.block,dir_upwindsea,interval=10)
{
  distance<-20000
  dem.m<-getValues(dem.block,format="matrix")
  upwind.block<-array(NA, dim=c((360/interval),nrow(dem.m),ncol(dem.m)))
  for (i in 1:(360/interval)){
    in.file<-paste(dir_upwindsea,"upwind_5kmcell_in_dir_",direction,"_dist_",distance,".tif",sep="")
    print(in.file)
    upwind.r<-raster(in.file) 
    upwind.r<-crop(upwind.r,dem.block)
    upwind.block[i,,]<-getValues(upwind.r,format="matrix") # fills shelter[i,1:end,1] to shelter[i,1:10,end]
    print(paste("i= ",i," dir= ",dir))
  }
  return(upwind.block)
}# end function

#####################################################################
