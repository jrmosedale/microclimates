############################################
# FUNCTIONS
############################################
binary<-function(i) { # creates binary from integer
  a<-2^(0:9)
  b<-2*a
  sapply(i,function(x) sum(10^(0:9)[(x %% b)>=a]))
}

binary8<-function(i) { # creates text string of 8 bit binary from integer
  a<-2^(0:9)
  b<-2*a
  bin.c<-format(sapply(i,function(x) sum(10^(0:9)[(x %% b)>=a])),scientific=FALSE)
  if (nchar(bin.c)>8) warning("Integer > 8 bit binary")
  if (nchar(bin.c)<8) bin.c<-sprintf("%08s",bin.c,sep="") # add starting zeros to binary
  return(bin.c)
}

BinToDec <- function(x) 
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))

# Function to get arr.ind from matrix cell number
matrix.xy<-function(m,i){
  xy<-c(0,0)
  xy[2]<-ceiling(i/length(m)); print(xy[2])
  xy[1]<-i-(xy[2]-1)*length(m); print(xy[1])
  #print(paste(i,xy[1],xy[2]))  
  return(xy)
}

# Function return binary string from 3x3 matrix
lowerthan.int<-function(m){
  hgt<-rep(m[2,2],8)
  neighbours<-c(m[1,2],m[1,3],m[2,3],m[3,3],m[3,2],m[3,1],m[2,1],m[1,1])
  #print(neighbours)
  lower.than.neighbour<-ifelse(hgt>=neighbours,1,0)
  #print(lower.than.neighbour)
  int.code<-sum(2^(which(rev(unlist(strsplit(as.character(lower.than.neighbour), "")) == 1))-1)) # convert binary to integer
  #print(int.code)
} # end function

############################################
# Record cell height to 8 adjecent neighbours 
# Write raster decribing if higher/lower than neighbouring cells test for single block
# Records as integer of 8 bit binary number
e.dem<-extent(c( 120000,420000,10000,180000 )) # excludes scilly isles
#e.dem<-c(70000,350000,0,160000)
dem<-crop(demuk,e.dem)
m<-getValues(dem,format="matrix")

# put buffer around m
m2<-array(9999,dim=c(dim(dem)[1]+2,dim(dem)[2]+2)) # set boundary to very high value - to ensure basins do not extend beyond area of interest
m2[2:(dim(dem)[1]+1),2:(dim(dem)[2]+1)]<-m
m<-m2

# Create height comparison matrix
lowerthan.m<-matrix(rep(NA,length(dem)),nrow=dim(dem)[1])
for (y in 2:(dim(m2)[2]-1) ) {
  print(y)
  for (x in 2:(dim(m2)[1]-1) ) {
    focalcell<-m[x,y]
    if (!is.na(focalcell)==TRUE) {
      m9<-matrix(c(m[x-1,y-1],m[x,y-1],m[x+1,y-1],m[x-1,y],focalcell,m[x+1,y],m[x-1,y+1],m[x,y+1],m[x+1,y+1]),nrow=3)
      lowerthan.m[(x-1),(y-1)]<-lowerthan.int(m9)
    }# end if
  }# end for y
}# end for x

# Write raster
lowerthan.r<-raster(lowerthan.m,template=dem)
plot(lowerthan.r)
file.out<-paste(dir.basinmap,"lowerthanmap.tif",sep="")
writeRaster(lowerthan.r,file=fileout,overwrite=TRUE)

#############
# WRITE BASINS RASTER 
  
# Create vectors to modify row/col to that of each of 8 neighbouring cell (order=clockwise from north) 
ngb.row<-c(-1,-1,0,1,1,1,0,-1)
ngb.col<-c(0,1,1,1,0,-1,-1,-1)
  
# Create matrices from rasters
lowerthan.r<-raster(paste(dir.basinmap,"lowerthanmap.tif",sep=""))
e<-extent(120000,190000,10000,50000)
lowerthan.r<-crop(lowerthan.r,e)
lowerthan.m<-getValues(lowerthan.r,format="matrix")
dem.m<-getValues(crop(dem,e),format="matrix")
plot(lowerthan.r,main="lowerthan")
plot(crop(dem,e),main="dem") 

# Create results matrices 
basins.m<-ifelse(is.na(dem.m),NA,0) # NA=sea 0=unassigned to basin
done.m<-dem.m # will be assigned as NA when analysed
# plot(raster(basins.m,template=r))

# Set starting values etc
unassigned<-order(dem.m) # low to high
basin<-1
#cell<-1

while ( length(unassigned)>=1 & !is.na(done.m[unassigned[1]]) ) {
  basincell<-unassigned[1] # start with lowest unassigned cell

  while ( basincell!=-999)  { # while unchecked cells remain in basin
    i<-basincell # assign lowest remaining cell index to i
    
    # If cell is unassigned calculate basin
    if (!is.na(done.m[i]) ) { # if cell has not been focal cell
      # calculate matrix x,y for cell
      i.row<-arrayInd(i,dim(dem.m))[1]
      i.col<-arrayInd(i,dim(dem.m))[2]
      basins.m[i.row,i.col]<-basin
      print(paste("Cell=",cell," i=",i,"Basin= ",basin,"Height= ",dem.m[i]," Neighbour Index= ",lowerthan.m[i]))
      neighbours<-binary8(lowerthan.m[i])
      #print(neighbours)
      
      # assign neighbours of higher elevation NOT already assigned to same basin
      for (n in 1:8){
        if( (i.row+ngb.row[n])>0 & (i.row+ngb.row[n])<=dim(basins.m)[1] & (i.col+ngb.col[n])>0 & (i.col+ngb.col[n])<=dim(basins.m)[2] ) { # check within limits of raster
          if ( !is.na(basins.m[(i.row+ngb.row[n]),(i.col+ngb.col[n])])) { # check if neighbour = sea
          if (substr(neighbours,n,n)=="0" & basins.m[(i.row+ngb.row[n]),(i.col+ngb.col[n])]==0 ) {
            basins.m[(i.row+ngb.row[n]),(i.col+ngb.col[n])]<-basin
            } # end if
          }
        } # end if 
      } # end for
      
    } # if basin unassigned
    done.m[i]<-NA # record that cell has been treated as focal cell
    
    # Recalculate next basincell remaining to be checked and order by height
    if (length(which(basins.m==basin & !is.na(done.m)))>0)  basincell<-min(which(basins.m==basin & !is.na(done.m) )) else basincell<--999
    #print(paste("Basincell= ",basincell))
  } # end while unchecked cells remain

  unassigned<-order(done.m) # resorts unchecked cells
  if (basin%%50==0) plot(raster(basins.m,template=r)) ; # CHANGE THIS to DEM
  basin<-basin+1
} # end while unassigned cells exist new basin

basins.r<-raster(basins.m,template=crop(dem,e))
plot(basins.r,main="Basins")
file.out<-paste(dir.basinmap,"basinmap-test1.tif",sep="")
print(paste("Writing ",file.out,sep=""))
writeRaster(basins.r,file=fileout,overwrite=TRUE)


#####################################################
# Compare basin methods - e sets extent of plot
e<-extent(140000,150000,30000,40000)
e<-extent(120000,190000,10000,50000)

file.in<-paste(dir.basinmap,"all.tif",sep="")
print(file.in)
oldbasins.r<-raster(file.in)
plot(crop(oldbasins.r,e),main="old basins", col=rainbow(4000))
plot(crop(basins.r,e),main="Basins",col=rainbow(4000))

num.basins1<-length(unique(crop(oldbasins.r,e)))
num.basins2<-length(unique(crop(basins.r,e)))
print(paste("Number of old basins= ",num.basins1," Number of new basins= ",num.basins2))


  
