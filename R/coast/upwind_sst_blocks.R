# Functions to calculate Upwind SST and SST-Tref for each 100m cell in block given wind direction and sst map
# Input:  100m dem grid -  block and buffered regions
#         raster of sst for buffer attime  t
#         raster of wind direction for block at time t
# Output: raster of SST and SST-Tref at t

#####################################################################
# FUNCTIONS 
# Used by upwind.sst
findsst<-function(x)
{
  is.sea<-ifelse(x==0,NA,1)
  nearest.sea<-match(1,is.sea)
  if (is.na(nearest.sea)){
    seacell<-0 
  } else {seacell<-x[nearest.sea]}  
  return(seacell)
}
#####################################################################
# Finds nearest upwind sea cell within 'buffer' (20km) for a block of cells for a SINGLE wind direction 
# If no sea cell found returns value of...0 sea cells returned = NA
# Otherwise records upwind sst held in sst.r
# Input:  sstbuffer.r - sea surface temperatures for buffer region
#         gridbuffer.r (including buffer) - could be dem
#         gridblock.r - defines area block of interest within buffer - could be dem
#         direction assumed to be constant across area
#         distance = max distance a search for nearest sea cell (= buffer of 20km)
# called by: upwind.sst.block
upwind.sst<-function(sstbuffer.r,gridbuffer.r,gridblock.r,direction,distance=20000)
{
  x<-dim(sstbuffer.r)[1]
  y<-dim(sstbuffer.r)[2]
  step<-distance/res(sstbuffer.r)[1] # max number of cells from focal cell to be searched
  
  # create matrix holding cell numbers for sea cells and 0 for landcells
  vals<-getValues(gridbuffer.r,format="matrix") # land/sea 1/NA values
  cells<-getValues(sstbuffer.r,format="matrix")
  seacells<-ifelse(is.na(vals),cells,0) # matrix of -999 if land or sea temperature if sea - could be changed to sst values as long as land =0
  #plot(raster(seacells,template=gridbuffer.r))
  
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
  
  upwind.sst.r<-raster(upwind.sst,template=gridblock.r)
  upwind.sst.r<-raster::mask(upwind.sst.r,gridblock.r)
  #plot(upwind.sst.r,main=paste("nearest seacell where direction= ",direction,sep=""))
  
  return(upwind.sst.r)
}# end function

#####################################################################

# Function to return upwind.sst for block given w.dir for block
upwind.sst.block<-function(wdir.block,sst.buffer,dem.buffer,dem.block)
{
  # define output raster
  sst.vals<-array(NA,dim=c(nrow(wdir.block),ncol(wdir.block)))
  # Round wind direction to nearest degree
  wdir.block<-round(wdir.block,0)
  # Calc number of unique wind dir - create raster stack of upwind.sst for each unique w.direction
  wdir.vals<-unique(round(getValues(wdir.block)[which(!is.na(getValues(wdir.block)))])) # every unique wind.val excluding NA and rounding to nearest degree
  
  # Create stack of upwind sst for block for every unique wind direction udirng time t
  upwind.vals<-array(0,dim=c(nrow(wdir.block),ncol(wdir.block),length(wdir.vals)))  
  for (i in 1:length(wdir.vals)){
    upwind.r<-upwind.sst(sst.buffer,dem.buffer,dem.block,wdir.vals[i],distance=20000)
    upwind.vals[,,i]<-getValues(upwind.r, format="matrix")
  }
  # assign sst according to w.dir of cell
  wdir.class<-cbind(wdir.vals,1:length(wdir.vals)) # create 2col vector for reclassifying
  wdir.class.r<-reclassify(wdir.block,wdir.class)  #reclass so that values = layer of stack holding sst values 
  class.vals<-getValues(wdir.class.r, format="matrix")

  for (i in 1:length(wdir.vals)){
    sel<-which(class.vals==i, arr.ind=TRUE)
    upwind.i<-upwind.vals[,,i]
    sst.vals[sel]<-upwind.i[sel]
    #plot(raster(sst.vals,template=wdir.block))
  }
  result<-raster(sst.vals,template=wdir.block)
  plot(result)
  return(result)
} # end function

#####################################################################
# Function to calculate sst-tref 
sst.tref<-function(sst.r,tref.r)
{
    if (compareRaster(sst.r,tref.r)!=TRUE){warning("!!sst-tref rasters not comparable!!")}
    sst.tref.r<-sst.r-tref.r
    return(sst.tref.r)
} # end function

#####################################################################

#####################################################################
# Input files for this hour at 100m resolution 

# If new day then input SST data - 5km hrly - use single file for each day NOT hour
if (paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")!=infile.sst){ 
  infile.sst<- fileout<-paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")
  sst5km.r<-raster(infile.sst)  
  sst.r<-resample(sst5km.r,dem.buffer)
}

# Input wind dir data - 100m hrly for block - convert to 5km
infile.wind<-paste(dir_winddirection,"direction_",year,"_",month,"_",day,"_",hr,".tif",sep="")
#infile.wind<-paste(dir_winddirection,"direction_2010_1_2010_12.tif",sep="")
wdir.r<-raster(infile.wind)
#wdir5km.r<-aggregate(wdir.r,50) # convert to 5km grid

# Input temperature data - whole area 5km hrly matrix
infile.tmp<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"-",sprintf("%02d",hr,sep=""),"00.tif", sep="") # define file name from year,month,day,hr
tref5km.r<-raster(infile.tmp)
# downscale to 100m

# Function calls
sst.block<-upwind.sst.block(wdir.block,sst.buffer,dem.buffer,dem.block)
sst.tref.block<-sst.tref(sst.block,tref.block)





######## DRAFT FUNCTION - NOT WORKING ########
# Function to return upwind.sst for block given w.dir for block
upwind.sst.block2<-function(wdir.block,sst.buffer,dem.buffer,dem.block)
{
  # Round wind direction to nearest degree
  wdir.block<-round(wdir.block,0)
  # Create stack of upwind sst for block for every unique wind direction udirng time t
  upwind.sst.stack<-stack()
  # Calc number of unique wind dir - create raster stack of upwind.sst for each unique w.direction
  wdir.vals<-unique(round(wdir.vals[which(!is.na(wdir.vals))])) # every unique wind.val excluding NA and rounding to nearest degree
  for (i in wdir.vals){
    upwind.sst.stack<-stack(upwind.sst.stack,upwind.sst(sst.buffer,dem.buffer,dem.block,i,distance=20000))
  }
  # plot(upwind.sst.stack)
  
  # assign sst according to w.dir of cell
  wdir.class<-cbind(wdir.vals,1:length(wdir.vals)) # create 2col vector for reclassifying
  wdir.class.r<-reclassify(wdir.block,wdir.class)  #reclass so that values = layer of stack holding sst values 
  result<-overlay(upwind.sst.stack,wdir.class.r,fun=function(x,y) { raster(x,layer=y) } ) # create upwindsst raster for block
  
  return(result)
} # end function

#####################################################################