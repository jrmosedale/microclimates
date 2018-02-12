############################################
# FUNCTIONS
############################################
# Creates text string of an 8 bit binary from an integer
integer.to.binary8<-function(i) { 
  a<-2^(0:9)
  b<-2*a
  bin.c<-format(sapply(i,function(x) sum(10^(0:9)[(x %% b)>=a])),scientific=FALSE)
  if (nchar(bin.c)>8) warning("Integer > 8 bit binary")
  if (nchar(bin.c)<8) bin.c<-sprintf("%08s",bin.c,sep="") # add starting zeros to binary
  return(bin.c)
}

# Function returns integer describing lower/higher than neighbouring cells
# Input:  3x3 matrix 
# Describes difference between focal cell and each of 8 neighbours (clockwise from N) as 8 bit binary
# Converts binary to integer
hgt.to.ngbr<-function(m){
  hgt<-rep(m[2,2],8)
  neighbours<-c(m[1,2],m[1,3],m[2,3],m[3,3],m[3,2],m[3,1],m[2,1],m[1,1])
  #print(neighbours)
  higher.than.neighbour<-ifelse(hgt>=neighbours,1,0)
  #print(lower.than.neighbour)
  int.code<-sum(2^(which(rev(unlist(strsplit(as.character(higher.than.neighbour), "")) == 1))-1)) # convert binary to integer
  #print(int.code)
} # end function

############################################
# Record cell height to 8 adjacent neighbours 
# Writes raster decribing if higher/lower than neighbouring cells
# Records as integer of 8 bit binary number

root<-"~/Documents/Exeter/Data2015/"
dir.basinmap2<-paste(root,"basins/maps2/maps2",sep="")

e<-extent(c( 120000,420000,10000,180000 )) # excludes scilly isles
#e<-c(70000,350000,0,160000)
dem<-crop(demuk,e)
m<-getValues(dem,format="matrix")

# put buffer around m
m2<-array(9999,dim=c(dim(dem)[1]+2,dim(dem)[2]+2)) # set boundary to very high value - to ensure basins do not extend beyond area of interest
m2[2:(dim(dem)[1]+1),2:(dim(dem)[2]+1)]<-m
m<-m2

# Create height comparison matrix
updown.m<-matrix(rep(NA,length(dem)),nrow=dim(dem)[1])
for (y in 2:(dim(m2)[2]-1) ) {
  print(y)
  for (x in 2:(dim(m2)[1]-1) ) {
    focalcell<-m[x,y]
    if (is.na(focalcell)==FALSE) {
      m9<-matrix(c(m[x-1,y-1],m[x,y-1],m[x+1,y-1],m[x-1,y],focalcell,m[x+1,y],m[x-1,y+1],m[x,y+1],m[x+1,y+1]),nrow=3)
      updown.m[(x-1),(y-1)]<-hgt.to.ngbr(m9)
    }# end if
  }# end for y
}# end for x

# Write raster
updown.r<-raster(updown.m,template=dem)
plot(updown.r)
fileout<-paste(dir.basinmap,"updownmap.tif",sep="")
writeRaster(updown.r,file=fileout,overwrite=TRUE)

#############################################
# WRITE BASINS RASTER 
# Input: raster describing if up/down to 8 nearest neighbours; dem
#############################################
ptm<-proc.time()
# Create vectors to modify row/col to that of each of 8 neighbouring cell (order=clockwise from north) 
ngb.row<-c(-1,-1,0,1,1,1,0,-1)
ngb.col<-c(0,1,1,1,0,-1,-1,-1)
  
updown.r<-raster(paste(dir.basinmap,"updownmap.tif",sep=""))

# Set extent, crop and create matrices from rasters
#e<-extent(120000,190000,10000,50000) #(SW peninsular)
e<-extent(175000,215000,60000,90000)  # Camel valley
#e<-extent(c( 120000,420000,10000,180000 )) # excludes scilly isles

updown.r<-crop(updown.r,e)
updown.m<-getValues(updown.r,format="matrix")
dem.m<-getValues(crop(dem,e),format="matrix")

plot(updown.r,main="height to neighbours index")
plot(crop(dem,e),main="dem") 

# Create results matrices 
basins.m<-ifelse(is.na(dem.m),NA,0) # NA=sea 0=unassigned to basin
done.m<-dem.m # will be assigned as NA when analysed

# Set starting values etc
unassigned<-order(dem.m) # low to high
basin<-1

# Choose lowest unassigned cell  
while ( length(unassigned)>=1 & !is.na(done.m[unassigned[1]]) ) {
  basincell<-unassigned[1] 
  
  # While unanalysed cells in basin
  while ( basincell!=-999)  { 
    i<-basincell 
    
    # If unassigned then assign cell i to basin
    if (!is.na(done.m[i]) ) { 
      i.row<-arrayInd(i,dim(dem.m))[1]
      i.col<-arrayInd(i,dim(dem.m))[2]
      basins.m[i.row,i.col]<-basin
      print(paste("Cell=",cell," i=",i,"Basin= ",basin,"Height= ",dem.m[i]," Ngb Index= ",updown.m[i] ))

      # assign neighbours of higher elevation  to same basin
      # NOTE only does so if neighbour is NOT already assigned to previous basin
      neighbours<-integer.to.binary8(updown.m[i])
      for (n in 1:8){
        # Check if neighbour within limits of raster
        if( (i.row+ngb.row[n])>0 & (i.row+ngb.row[n])<=dim(basins.m)[1] & (i.col+ngb.col[n])>0 & (i.col+ngb.col[n])<=dim(basins.m)[2] ) { 
          # Check if neighbour = sea 
          if ( !is.na(basins.m[(i.row+ngb.row[n]),(i.col+ngb.col[n])])) { 
          # check if neighbour higher and unassigned to basin
          if (substr(neighbours,n,n)=="0" & basins.m[(i.row+ngb.row[n]),(i.col+ngb.col[n])]==0 ) {
            basins.m[(i.row+ngb.row[n]),(i.col+ngb.col[n])]<-basin
            } # end if
          } # end if
        } # end if 
      } # end for each neighbour
    } # if unassigned
    done.m[i]<-NA # record that cell has been used as focal cell
    
    # Identify next basincell (lowest) to be checked if any
    if (length(which(basins.m==basin & !is.na(done.m)))>0)  basincell<-min(which(basins.m==basin & !is.na(done.m) )) else basincell<--999
  } # end while unanalysed cells in basin

  # Prepare for analysis of new basin
  unassigned<-order(done.m) # resorts unchecked cells
  #if (basin%%100==0) plot(raster(basins.m,template=crop(dem,e))) ; # Plot every 100 basins
  basin<-basin+1
} # end while unassigned cells exist 

print(proc.time()-ptm)

basins.r<-raster(basins.m,template=crop(dem,e))
plot(basins.r,main="Basins")

fileout<-paste(dir.basinmap,"basinmap-all.tif",sep="")
print(paste("Writing ",fileout,sep=""))
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


  
