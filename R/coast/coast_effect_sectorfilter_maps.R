# Works out costal index based on wind direction and distance from coast
# Inputs: dem = a digital eleveation model  - and extent including  buffer region around area of interest (=coastal ef distance)
#         wind angle, sector width, coastal effect distance (sector radius), res of grid
# Outputs: maps of coastal effect index (0:1) for each wind direction
# IMPORTANT: calculates using euclidian geometry???? 

library(sp)
library(ncdf4)
library(raster)
library(rgdal)

#####################################################################
# Testing - do not run
#angles<-seq(1,360,1)
#testx<--50; testy<-0
#areClockwise2(testx,testy,test[,1],test[,2]); print(paste(testx,testy,sep=" , "))
#atan(testx/testy)/(pi/180)
#####################################################################

# FUNCTION -  calculate end x/y coordinates for vector 
# INPUT: angle (degrees) and radius of vector
# OUTPUT: xy coordinates 
vectorxy <- function(angle,radius) #
{ angle<-angle*(pi/180)
  x<-sin(angle)*radius
  y<-cos(angle)*radius
  xy<-cbind(x=x,y=y)
  return(xy)
}

# FUNCTION - to calculate if vector v2 clockwise from vector v1 defined using x/y coordinates
# INPUT: v1 x & y coord, v2 x & y 
# OUTPUT: TRUE=clockwise
# ASSUMPTION: assumes vectors within 180 degres - false results if obtuse angle
# Ref: http://stackoverflow.com/questions/13652518/efficiently-find-points-inside-a-circle-sector
areClockwise<-function(v1.x,v1.y,v2.x,v2.y){ # is v2  clockwise of  v1 
  clockwise<- (-v1.x*v2.y )+(v1.y*v2.x) >0 
  return(clockwise)
}

# FUNCTION - to calculate if vector v2 clockwise from vector v1 defined using x/y coordinates
# INPUT: v1 x & y coord, v2 x & y 
# OUTPUT: TRUE=clockwise  Identical vectors=FALSE
# ASSUMPTION: assumes vectors within 180 degres - false results if obtuse angle
# EXPANDED version of above
areClockwise2<-function(v1.x,v1.y,v2.x,v2.y){ # is v2 clockwise of  v1 
  nv.x<--v1.y
  nv.y<-v1.x
  projx<-v2.x*nv.x  ;#print(projx); plot(projx)
  projy<-v2.y*nv.y ;#print(projy); plot(projy)
  clockwise<- (projx+projy) <0 ; #print(projx+projy);plot(projx+projy)
  return(clockwise)
}

# FUNCTION - to calculate filter for coastal effect
# INPUT: blank matrix of results, radius (distance of coastal effect in m), start/end vectors
# Calculates inverse distance measure
# ALL cells for which any of the 4 corners fall within sector are included in filter
# OUTPUT: filled matrix with weighted index that sums to 1.0 for all m
# Ref: https://scrogster.wordpress.com/2012/10/05/applying-a-circular-moving-window-filter-to-raster-data-in-r/

invdist.sector.filter<-function(m, radius, res=100, v1, v2) {  
  # calculate end xy for vectors defining sector
  v1xy<-vectorxy(v1,radius) # start of sector vector
  v2xy<-vectorxy(v2,radius) #  end of sector (clockwise of start) vector
  #print(paste("sector = ",v1,v2,sep="  "))
  #print(paste("v1xy=",v1xy,sep=""))
  #print(paste("v2xy=",v2xy,sep=""))
  
  for (row in 1:nrow(m)){
    for (col in 1:ncol(m)){
      y<-as.numeric(dimnames(m)[[1]])[row] 
      x<-as.numeric(dimnames(m)[[2]])[col] 
      ymax<-y+(res/2) ; ymin<-y-(res/2) 
      xmax<-x+(res/2) ; xmin<-x-(res/2) 
        
      dist<-sqrt(y^2 + x^2)
      #cwv2<-areClockwise(v2xy[,"x"],v2xy[,"y"],x,y)
      #ccwv1<-areClockwise(x,y,v1xy[,"x"],v1xy[,"y"])
      # Check to see if angle with either v1 or v2 is >90degrees
      if (  c(v1xy[,"x"],v1xy[,"y"]) %*% c(x,y) <=0 |  c(v2xy[,"x"],v2xy[,"y"]) %*% c(x,y) <=0)  obtuse=TRUE else obtuse<-FALSE
      # if (floor(x/1000)==ceiling(x/1000) & floor(y/1000)==ceiling(y/1000)) print(paste(x," ", y,"  ", dist,"    ", obtuse))
      # Check if grid cell falls into sector
      if ( (dist<=radius) & (obtuse==FALSE) &
           ( areClockwise2(v2xy[,"x"],v2xy[,"y"],xmax,ymax ) |
             areClockwise2(v2xy[,"x"],v2xy[,"y"],xmax,ymin ) |
              areClockwise2(v2xy[,"x"],v2xy[,"y"],xmin,ymax ) |
               areClockwise2(v2xy[,"x"],v2xy[,"y"],xmin,ymin ) )  &
           ( !areClockwise2(v1xy[,"x"],v1xy[,"y"],xmax,ymax ) |
              !areClockwise2(v1xy[,"x"],v1xy[,"y"],xmax,ymin ) |
               !areClockwise2(v1xy[,"x"],v1xy[,"y"],xmin,ymax ) |
                !areClockwise2(v1xy[,"x"],v1xy[,"y"],xmin,ymin ) )
               )
      { m[row,col]<-1/(dist/radius)  }     
    }
  }
  #sect.cells<-length(which(m>0))
  #index.m<-m/sum(m) # sum of index.m=1.0 (if all sector cells qualify ie are sea)
  index.m<-m/sum(m) ; print(sum(m))
  return(index.m)
}


########################################################################################
# Define coastal effect in metres and resolution of rasters in metres
radius<-10000 # distance of coastal effect
res<-100 # resolution of raster
dir_coast<-"~/Documents/Exeter/Data2015/CoastEffect/"
########################################################################################
# Define DEM required for calculations
#dem<-raster("C:/Data2015/DEM100/demoriginal.tif")
demuk<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif", crs=("+init=epsg:27700"))

# Define sw dem of interest (no buffer)
#e.dem<-extent(c(70000,420000,0,180000 )) # includes scilly isles
e.dem<-extent(c( 130000,400000,10000,180000 )) # excludes scilly isles
dem<-crop(demuk,e.dem)

# Create buffered UK dem
# SIMPLIFICATION: REQUIRE land outline/raster for whole UK/N.France 
demukb<-extend(demuk,c((radius/res),(radius/res)),value=NA) # extend 20km of "sea" around dem

# 1.Create buffered south west dem raster required - USED to create coastal effect index maps
# Sea cells = NA, land cells = elevation
buffer<-20000 # include 20km buffer area around zone of interest
#e.dem<-extent(c((70000-buffer),(420000+buffer),(0-buffer),(180000+buffer) )) # includes scilly isles
e.buf<-extent(c( (xmin(dem)-buffer),(xmax(dem)+buffer),(ymin(dem)-buffer),(ymax(dem)+buffer) )) # excludes scilly isles
demsw<-crop(demukb,e.buf)
plot(demsw,main="DEM-sw")
#e.dem <-extent(demsw)

# 2. Create simplified land/sea raster

landsea.buffer<-calc(dembuf,function(x) ifelse(is.na(x),0,1))
plot(landsea.buffer)

# 3. For every wind direction (1-360  in steps of  10 degrees) calculate a coast index map

for (angle in seq(190,360,10))  {

#  angle<-50 # e.g direction from which wind blowing or... 
  
# Define sector width = 1+ degrees and resulting start/end vectors v1/v2
sect.width<-10
sect.rad<-sect.width/2
if (angle>5) {v2<-angle-sect.rad
} else {v2<-angle-sect.rad}
if (angle<356) {v1<-angle+sect.rad
} else {v1<-(angle-(360-sect.rad))}

#####################################################################
# Calculate filter matrix

# Define results matrix - ATTENTION - ORIENTATION!!!
m<-matrix(0,nrow=1+(2*radius/res),ncol=1+(2*radius/res))
  dimnames(m)[[1]]<-rev(seq(-radius,radius,by=res))
  dimnames(m)[[2]]<-seq(-radius,radius,by=res)

# Calculate filter/index
#filter<-sector.filter(m,radius,v1,v2) 
filter<-invdist.sector.filter(m,radius,res=100,v1,v2) 
filter.r<-raster(filter)
plot(filter.r,main=angle)
#####################################################################

# Apply to full raster and confirm sea cells as NA via mask
new.r<-focal(landsea.buffer,filter)
# new.r<-focal(land.buffer,filter,fun=function(x,y){ifelse(y!=0,grid)}
#plot(new.r); 
coast.r<-mask(new.r, dembuf)
coast.r<-crop(coast.r,dembuf)
plot(coast.r,main=paste("Coast Effect where wind dir= ",angle))

#####################################################################
# Write coastal index map as file
out.file<-paste(dir_coast,"Coast_Index_Sector_invdist_",sprintf("%03d",angle,sep=""),"_",(radius/res),"km_",sect.width,"deg.r",sep="")
print(out.file)
save(coast.r,file=out.file)

} # end for angle loop

#####################################################################
# END HERE




########################################################################################
# CUT OUTS and TESTING
########################################################################################
# Apply filter to test raster to check effects

# Create test square landraster
testm<-matrix(1,nrow=1+(4*radius/res),ncol=1+(4*radius/res))
n1<-((ncol(testm)-radius/res))
n0<-(radius/res)*1.2
testm[n0:n1,n0:n1]<-0
test.r<-raster(testm,xmn=0,xmx=1000000,ymn=0,ymx=1000000)
plot(test.r)
# Apply filter 
result<-focal(test.r,filter)
testmask.r<-test.r
testmask.r[testmask.r==1]<-NA
result2<-mask(result, testmask.r)
plot(result2)
#####################################################################
# Apply filter to reduced size map where sea=1, land=0
#e<-extent(110000,260000,-10000,110000)
#tmp.r<-crop(landsea.r,e)
#plot(tmp.r)
#new1.r<-focal(tmp.r,filter,fun=sum,na.rm=TRUE) # where sea=1, land=0 
#plot(new1.r)
#mini.r<-mask(new1.r, crop(demsw,e))
#plot(mini.r)  
#####################################################################



# new = old except where mask=maskvalue then become updateval
new.r<-mask(old.r,mask.r,maskvalue= , updatevalue=)
overlay(r1,r2,fun=function(r1,r2){ commands})
substitute(old.r,datafr,by=col , which= col, subsWithNA=FALSE (no match keeps old val), )
#replacement: r[cellnum]<-  

cellFromXY(raster,xy)
rasterFromXYZ(ci,res=100,)

# Trig / math etc
sector.area<-(pi/360)*pi*r^2
sector.area<-((arc.length/circumf)/2*pi*r)*pi*r^2

rasterize(spl[1],landsea.r,fun=sum)

# Function calculates coastal index from start/end xy and a background sea=1, land=0 raster map
# INPUT:  starting xy; ending xy, landsea raster
# OUTPUT: array of index values (0:1) in same order as start/end xy
coastindex <- function(start.xy,end.xy,landsea) {
  linelist<-vector("list",nrow(start.xy)) # create list to stores series of 2*2 matrices 
  for (l in 1:nrow(start.xy)){
    linelist[[l]] <- Lines(list(Line(rbind(start.xy[l,],end.xy[l,]))), as.character(l))
    #print(linelist[[l]])
  }
  spl<-SpatialLines(linelist)
  print("Line data complete - now calc index")
  # calculate coast index - number of sea cells on line/ total number of cells on line
  c.index<-extract(landsea,spl,fun=function(x,...)sum(x)/length(x))
  return(c.index)
}