# Calculate SST-Tref for each 100m cell for time t
# Input:  raster of temperature data at t
#         raster of sst for t
#         wind direction for t
# Output: raster of SST-Tref at t


#####################################################################
# Buffered and normal 5km landsea grid = grid5kmbuf.r, grid5km.r
# Buffered and normal 100m  dem = dembuf, dem
#####################################################################

findsea<-function(x)
{
  is.sea<-ifelse(x==0,NA,1)
  nearest.sea<-match(1,is.sea)
  if (is.na(nearest.sea)){
    seacell<-0 
  } else {seacell<-x[nearest.sea]}  
  return(seacell)
}

# Set landsea.r to include buffer region 
# Finds nearest upwind sea cell within 'buffer' (30km) for cells in gridmask.r
# If no sea cell found returns value of...0
# Otherwise records upwind cell number in gridmask.r
# Input:  grid5kmbuf.r (including buffer)
#         grid5km.r - defines ara of interest within buffer
#         direction assumed to be constant across area
#         distance = max distance a search for nearest sea cell (= buffer of 20km)

upwind.sea<-function(grid5kmbuf.r,grid5km.r,direction,distance=20000)
{
  x<-dim(grid5kmbuf.r)[1]
  y<-dim(grid5kmbuf.r)[2]
  step<-distance/res(grid5kmbuf.r)[1] # max number of cells from focal cell to be searched
  
  # create matrix holding cell numbers for sea cells and 0 for landcells
  vals<-getValues(grid5kmbuf.r,format="matrix") # land/sea 1/NA values
  cells<-getValues( rasterFromCells(grid5kmbuf.r,1:ncell(grid5kmbuf.r),values=TRUE), format="matrix") 
  seacells<-ifelse(is.na(vals),cells,0) # matrix of 0 if land or cell number if sea
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
  
  upwind.sea<-apply(storev,1,findsea) # if no sea within buffer then = NA?
  upwind.sea<-matrix(upwind.sea,nrow=nrow(store),ncol=ncol(store)) 
  
  upwind.sea.r<-raster(upwind.sea,template=grid5km.r)
  upwind.sea.r<-mask(upwind.sea.r,grid5km.r)
  plot(upwind.sea.r,main=paste("nearest seacell where direction= ",direction,sep=""))
  
  return(upwind.sea.r)
}# end function

#####################################################################
# Function to create maps recording nearest sea cell for different wind directions
# could create 5km raster block 1-36 layers according to wind dir
# record cell number of nearest sea cell!!
#PROBLEM - check cell values point to correct cell in matrix/raster

write.upwind.maps<-function(grid5kmbuf.r,grid5km.r,dir_upwindsea)
{
    for (direction in seq(0,350,10)){
      #print(direction)
      upwindsea.r<-upwind.sea(grid5kmbuf.r,grid5km.r,direction)
      # upwindsea100.r<-upwind.sea(gbuf100.r,direction,gmask100.r,buffer)
      out.file<-paste(dir_upwindsea,"upwind_5kmcell_in_dir_",direction,"_dist_",distance,".grd",sep="")
      print(out.file)
      writeRaster(upwindsea.r,file=out.file,overwrite=TRUE)
    }
}# end function


#####################################################################
# Function to load upwind maps into 3d matrix for block cells
block.upwind.sst<-function(dem.block,dir_upwindsea,interval=10)
{
    distance<-30000
    dem.m<-getValues(dem.block,format="matrix")
    upwind.block<-array(NA, dim=c((360/interval),nrow(dem.m),ncol(dem.m)))
    for (i in 1:(360/interval)){
      in.file<-paste(dir_upwindsea,"upwind_5kmcell_dir_",direction,".tif",sep="")
      print(in.file)
      upwind.r<-raster(in.file) 
      upwind.r<-crop(upwind.r,dem.block)
      upwind.block[i,,]<-getValues(upwind.r,format="matrix") # fills shelter[i,1:end,1] to shelter[i,1:10,end]
      print(paste("i= ",i," dir= ",dir))
    }
      return(upwind.block)
}# end function

#####################################################################

# Define time period
day<-9; month<-6;year<-2010; hr<-1 

# Input files for this hour

# If new day then input SST data - 5km hrly - use single file for each day NOT hour
if (paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")!=infile.sst)
  { infile.sst<- fileout<-paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")
    sst5km.r<-raster(infile.sst)  }

# Input wind dir data - 100m hrly for block - convert to 5km
infile.wind<-paste(dir_winddirection,"direction_",year,"_",month,"_",day,"_",hr,".tif",sep="")
wdir.r<-raster(infile.wind)
wdir5km.r<-aggregate(wdir.r,50) # convert to 5km grid

# Input temperature data - whole area 5km hrly matrix
infile.tmp<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"-",sprintf("%02d",hr,sep=""),"00.tif", sep="") # define file name from year,month,day,hr
tref5km.r<-raster(infile.tmp)


compare(sst5km.r,tref5km.r,widr5km.r,gridmask.r)

#####################################################################





sstnear.m<-upwind.sst(sst.r,wdir5km.r)
sstnear.r<- raster(sstnear.m,template=gridmask.r)
sst-tref.r<-sstnear.r-tref5km.r


extract array of cells upwind of focal
find first cell from focal not.na(sst)
nearest.raster.point(x,y,)
#dissagregate to 100m
