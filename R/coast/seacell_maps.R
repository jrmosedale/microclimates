direction<-50
distance<-10000

# FUNCTIONS 
# Used to calculate nearest sea cell
findsst<-function(x)
{
  is.sea<-ifelse(x==-999,NA,1)
  nearest.sea<-match(1,is.sea)
  if (is.na(nearest.sea)){
    seacell<-NA 
  } else {seacell<-x[nearest.sea]}  
  return(seacell)
}


#####################################################################
# Function to create maps recording nearest sea cell for different wind directions
# could create 5km raster block 1-36 layers according to wind dir
# record cell number of nearest sea cell!!
#PROBLEM - check cell values point to correct cell in matrix/raster

# If no sea cell found returns value of...0 sea cells returned = NA
# Otherwise records upwind sst held in sst.r
# Input:  sstbuffer.r - sea surface temperatures for buffer region
#         gridbuffer.r (including buffer) - could be dem
#         gridblock.r - defines area block of interest within buffer - could be dem
#         direction assumed to be constant across area
#         distance = max distance a search for nearest sea cell (= buffer of 20km)
#compareRaster(sst.r,grid5kmbuf.r)
### NB: IMPORTANT cell values referes to upsea cell in dembuffer proportioned raster

write.upwind.map<-function(dembuf,dem,direction,distance=10000,dir_upwindsea) {

# Create land/sea buffered DEM
landsea.buffer<-calc(dembuf,function(x) ifelse(is.na(x),0,NA))
# Create raster of cell numbers
cellvals<-c(1:ncell(landsea.buffer))
cellvals.m<-matrix(c(1:ncell(landsea.buffer)),nrow=nrow(landsea.buffer),ncol=ncol(landsea.buffer),byrow=TRUE)
#cellvals.m<-ifelse(is.na(landsea.buffer),cellvals.m,-999)
cellref.r<-raster(cellvals.m,template=landsea.buffer)
cellref.r<-mask(cellref.r,landsea.buffer)
plot(cellref.r)

# create matrix holding cell numbers for sea cells and -999 for landcells
res<-100
step<-distance/res(cellref.r)[1] # max number of cells from focal cell to be searched
x<-dim(cellref.r)[1]
y<-dim(cellref.r)[2]
seacells<-getValues(cellref.r,format="matrix")
seacells<-ifelse(is.na(seacells),-999,seacells)

store<-array(0,dim=c(dim(seacells)[1]-(2*step),dim(seacells)[2]-(2*step),step+1)) # 3d array to hold values for non-buffered region for every 'step'
store[,,1]<-seacells[(step+1):(dim(seacells)[1]-step),(step+1):(dim(seacells)[2]-step)]

for (i in 1:step)
{
  xshift<-round(i*sin(direction*(pi/180)),0) 
  yshift<-round(i*cos(direction*(pi/180)),0)
  #print(paste("i: ",i, " xshift: ",xshift," yshift: ",yshift,sep=""))
  yshift<-yshift*(-1)
  store[,,(i+1)]<-seacells[(step+1+yshift):(dim(seacells)[1]-step+yshift),(step+1+xshift):(dim(seacells)[2]-step+xshift)]
  
} # end

# Find nearest cell
storev<-array(store,dim=c((dim(store)[1]*dim(store)[2]),step+1))

upwind.cell<-apply(storev,1,findsst) # if no sea within buffer then = NA?
upwind.cell.m<-matrix(upwind.cell,nrow=nrow(store),ncol=ncol(store))
e<-extent(xmin(dembuf)+10000,xmax(dembuf)-10000,ymin(dembuf)+10000,ymax(dembuf)-10000)
upwind.cell.r<-raster(upwind.cell.m,template=crop(dembuf,e))
#plot(upwind.cell.r)
upwind.cell.r<-raster::mask(upwind.cell.r,crop(dembuf,e))
#upwind.cell.x<-xFromCell(landsea.buffer,upwind.cell.m)
#upwind.cell.y<-yFromCell(landsea.buffer,upwind.cell.m)
#upwind.cell.5kmcell<-
#upwind.5kmcell.m<-ifelse(is.na(upwind.cell.m),NA,cellFromRowCol(grid5kmbuf.r, rowFromCell(grid5kmbuf.r,), colnr)(,)

plot(upwind.cell.r,main=paste("Ref to upwind sea cell, angle= ",direction))

# Convert to cell number of 5km grid sea temperature raster
cells.100m<-getValues(upwind.cell.r)
xy.100m<-xyFromCell(landsea.buffer,cells.100m)
cells.5km<-cellFromXY(grid5kmbuf.r,xy.100m)
upwind.5kmcell.r<-setValues(upwind.cell.r,cells.5km)

plot(upwind.5kmcell.r,main=paste("Ref to 5km upwind sea cell, angle= ",angle))
# cellStats(upwind.sst.r,max)
return(upwind.5kmcell.r)
} # end function 


##############################
# WRITE MAP FILES FOR EACH WIND DIRECTION
#  direction<-10
  for (direction in seq(10,360,10)){
    print(direction)
    distance<-10000
    upwind.seacell.r<-write.upwind.map(dembuf,dem,direction,distance)
    # upwindsea100.r<-upwind.sea(gbuf100.r,direction,gmask100.r,buffer)
    out.file<-paste(dir_upwindsea,"upwind_5kmcellref_in_dir_",direction,"_dist_",distance,".tif",sep="")
    print(out.file)
    writeRaster(upwind.seacell.r,file=out.file,overwrite=TRUE)
  }

#####################################################################

# Use map files to identify nearest upwind sea temperature
direction<-0
sst.buffer<-
dem.block<-
infile.sst<- paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")
sst5km.r<-raster(infile.sst)  
#direction<-cellStats(wdir.block,mean)
cellref.block<-crop(raster(paste(dir_upwindsea,"upwind_5kmcellref_in_dir_",direction,"_dist_",distance,".tif",sep="")),dem.block)
cellref.vals<-getValues(cellref.block)
sst.vals<-extract(sst5km.r,cellref.vals)
sst.block<-setValues(dem.block,sst.vals)
plot(sst.block)
