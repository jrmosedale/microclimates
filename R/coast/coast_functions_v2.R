
##########################################################################################
# Calculates single raster  ldif.block
# Input: ldif.stack of ldif for each wind direction
#        wind direction for block
#        direction interval for which ldif calculated

calc.ldif.block<-function(ldif.stack,wdir.block,interval=10){
  ldif.block<-raster
  wdir.layer<-round(wdir.block/interval) # creates raster holding layer in ldif.stack for wind direction
  wdir.layer<-calc(wdir.layer,fun=function(x){ifelse(x==0,36,x)}) # corrects layer 0 to 36 (angle 0 to 360)
  ldif.block<-stackSelect(ldif.stack,wdir.layer)
  return(ldif.block)
} # end function
##########################################################################################

# Functions to calculate Upwind SST and SST-Tref for each 100m cell in block given wind direction and sst map
# Input:  100m dem grid -  block and buffered regions
#         raster of sst for buffer attime  t
#         raster of wind direction for block at time t
# Output: raster of SST and SST-Tref at t

#####################################################################
# FUNCTIONS 
# Used by upwind.sst
# Records temperature of nearest sea cell or records NA if no seacell withiin buffer
findsst<-function(x)
{
  is.sea<-ifelse(x==0,NA,1)
  nearest.sea<-match(1,is.sea) 
  if (is.na(nearest.sea)){
    seatemp<-NA 
  } else {seatemp<-x[nearest.sea]}  
  return(seatemp)
}

# Used by upwind.sst
# Records SST  temperature of nearest sea cell x distance/maxdistance
finddist<-function(x)
{
  is.sea<-ifelse(x==0,NA,1)
  nearest.sea<-match(1,is.sea)
  if (is.na(nearest.sea)) seadist<-0 
  if (is.na(nearest.sea)==FALSE) seadist<-(length(x)-nearest.sea)/(length(x)-1)   
  return(seadist)
}
#####################################################################
# UPWIND>SST BLOCK 
# Finds nearest upwind sea cell within 'buffer' (20km) for a block of cells for a SINGLE wind direction 
# If no sea cell found returns value of...0 sea cells returned = NA
# Otherwise records upwind sst held in sst.r
# Input:  sstbuffer.r - sea surface temperatures for buffer region
#         gridbuffer.r (including buffer) - could be dem
#         gridblock.r - defines area block of interest within buffer - could be dem
#         direction assumed to be constant across area
#         distance = max distance a search for nearest sea cell (= buffer of 20km)
# called by: upwind.sst.block

upwind.sst<-function(sst.buffer,dem.buffer,dem.block,direction,distance=20000)
{
  x<-dim(sst.buffer)[1]
  y<-dim(sst.buffer)[2]
  step<-distance/res(sst.buffer)[1] # max number of cells from focal cell to be searched
  
  # create matrix holding cell numbers for sea cells and 0 for landcells
  vals<-getValues(dem.buffer,format="matrix") # land/sea 1/NA values
  cells<-getValues(sst.buffer,format="matrix")
  seacells<-ifelse(is.na(vals),cells,0) # matrix of -999 if land or sea temperature if sea - could be changed to sst values as long as land =0
  #plot(raster(seacells,template=dem.buffer))
  
  store<-array(0,dim=c(dim(seacells)[1]-(2*step),dim(seacells)[2]-(2*step),step+1)) # 3d array to hold values for non-buffered region for every 'step'
  store[,,1]<-seacells[(step+1):(dim(seacells)[1]-step),(step+1):(dim(seacells)[2]-step)]
  #plot(raster(store[,,1],template=gridblock.r))
  
  for (i in 1:step)
  {
    xshift<-round(i*sin(direction*(pi/180)),0) 
    yshift<-round(i*cos(direction*(pi/180)),0)
    #print(paste("i: ",i, " xshift: ",xshift," yshift: ",yshift,sep=""))
    yshift<-yshift*(-1)
    store[,,(i+1)]<-seacells[(step+1+yshift):(dim(seacells)[1]-step+yshift),(step+1+xshift):(dim(seacells)[2]-step+xshift)]
    
  } # end
  
  # use first/last to find nearest sea cell??
  storev<-array(store,dim=c((dim(store)[1]*dim(store)[2]),step+1))
  
  upwind.sst<-apply(storev,1,findsst) # if no sea within buffer then = NA?
  upwind.sst<-matrix(upwind.sst,nrow=nrow(store),ncol=ncol(store)) 
  upwind.sst.r<-raster(upwind.sst,template=dem.block)
  upwind.sst.r<-mask(upwind.sst.r,dem.block)
  
  upwind.dist<-apply(storev,1,finddist) # if no sea within buffer then = NA?
  upwind.dist<-matrix(upwind.dist,nrow=nrow(store),ncol=ncol(store)) 
  upwind.dist.r<-raster(upwind.dist,template=dem.block)
  upwind.dist.r<-mask(upwind.dist.r,dem.block)
  #plot(upwind.sst.r,main=paste("nearest seacell where direction= ",direction,sep=""))
  #plot(upwind.dist.r,main=paste("Dist to seacell where direction= ",direction,sep=""))
  upwind.results<-c(upwind.sst.r,upwind.dist.r)
  return(upwind.results)
}# end function

#####

# Function to return upwind.sst for block given w.dir for block
upwind.sst.block<-function(wdir.block,sst.buffer,dem.buffer,dem.block, plotresult=FALSE)
{
  # define output rasters
  sst.vals<-array(NA,dim=c(nrow(wdir.block),ncol(wdir.block)))
  sst.dist<-array(NA,dim=c(nrow(wdir.block),ncol(wdir.block)))
  
  # Round wind direction to nearest degree
  wdir.block<-round(wdir.block,0)
  # Calc number of unique wind dir - create raster stack of upwind.sst for each unique w.direction
  wdir.vals<-unique(round(getValues(wdir.block)[which(!is.na(getValues(wdir.block)))])) # every unique wind.val excluding NA and rounding to nearest degree
  
  # Create stack of upwind sst for block for every unique wind direction during time t
  upwind.temp<-array(0,dim=c(nrow(wdir.block),ncol(wdir.block),length(wdir.vals)))  
  upwind.dist<-array(0,dim=c(nrow(wdir.block),ncol(wdir.block),length(wdir.vals)))  
  for (i in 1:length(wdir.vals)){
    upwind.results<-upwind.sst(sst.buffer,dem.buffer,dem.block,wdir.vals[i],distance=20000)
    upwind.temp[,,i]<-getValues(upwind.results[[1]], format="matrix")
    upwind.dist[,,i]<-getValues(upwind.results[[2]], format="matrix")
  }
  # assign sst according to w.dir of cell
  wdir.class<-cbind(wdir.vals,1:length(wdir.vals)) # create 2col vector for reclassifying
  wdir.class.r<-reclassify(wdir.block,wdir.class)  #reclass so that values = layer of stack holding sst values 
  class.vals<-getValues(wdir.class.r, format="matrix")
  
  for (i in 1:length(wdir.vals)){
    sel<-which(class.vals==i, arr.ind=TRUE)
    upwind.i<-upwind.temp[,,i]
    dist.i<-upwind.dist[,,i]
    sst.vals[sel]<-upwind.i[sel]
    sst.dist[sel]<-dist.i[sel]
    #plot(raster(sst.vals,template=wdir.block))
  }
  sst.r<-raster(sst.vals,template=wdir.block)
  seadist.r<-raster(sst.dist,template=wdir.block)
  if (plotresult) {
    plot(sst.r)
    plot(seadist.r)}
  
  results<-c(sst.r,seadist.r)
  return(results)
} # end function

#######
# Function to calculate (sst-tref) * distance index (seadistance/maxdist(10/20km)) 
sst.tref<-function(sst.r,tref.r,seadist.r)
{
  if (compareRaster(sst.r,tref.r)!=TRUE){warning("!!sst-tref rasters not comparable!!")}
  sst.tref.r<-overlay(sst.r,tref.r,seadist.r,fun=function(x,y,z){ifelse(is.na(x)&!is.na(y),0,(x-y)*z)})
  return(sst.tref.r)
} # end function