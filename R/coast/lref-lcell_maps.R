
# TASK 2. Calc 5kmref-cell difference in coast effect using RASTERS

# FUNCTIONS

# NOT USED - Calculate Lref from 100m cell numbers corresponding to 5km centre xy coordinates
# Method: mean of four 100m adjoiing cells (removing NA)
# Problem: where all 4 central cells =NA so Output raster^=5km hsitoric temp raster
L4celltoLref<-function(inv.lsratio,gridmask.r)
{ # define output raster lref.r
  lref5km.r<-raster(extent(gridmask.r),res=res(gridmask.r),nrows=nrow(gridmask.r),ncols=ncol(gridmask.r))
  lref5km.r<-setValues(lref5km.r,rep(0,ncell(lref5km.r)))
  lref<-rep(NA,ncell(lref5km.r))
  # Get centre xy coordinates for 5km grid cells
  cell5kmxy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
  # remove xy for NA cells
  #sel<-which(!is.na(getValues(gridmask.r)))
  #cell5kmxy<-cell5kmxy[sel,1:2]
  cells<-fourCellsFromXY(inv.lsratio,cell5kmxy)
  # define vals matrix and get values of four 100m cells at centre of 5km cell or NA if all 100m cells = NA
  vals<-matrix(NA,nrow=dim(cells)[1],ncol=4)
  for(n in 1:dim(cells)[1]) {
    vals[n,]<-extract(inv.lsratio,c(cells[n,]))
    if (all(is.na(vals[n,]))){ lref5km[n]<-NA
    } else{ lref5km<-rowMeans(vals,na.rm=TRUE)}
  }
  lref5km.r<-setValues(lref5km.r,lref5km)
  plot(lref5km.r)
  return(lref5km.r) 
} # end function

         

#####################################################################
#get valid 5km gridcells
dir_grids<-"~/Documents/Exeter/Data2015/Templates/"
#dir_grids<-"C:/Data2015/Templates/"
dir_lsratio<-"~/Documents/Exeter/Data2015/CoastEffect/lsratio/"
e.dem<-extent(dem)

in.file<-paste(dir_grids,"ukhistmask.grd",sep="")
print(in.file)
gridmask.r<-raster(in.file) #  1=valid cell, NA = sea or not data
gridmask.r<-crop(gridmask.r,e.dem) # IMPORTANT: crop to SAME geographical extent of DEM and inv.lsratio raster

#####################################################################

angle<-270
# For each wind direction and coastal effect map
for (angle in seq(0,350,10))
{
  
  # Load coast index maps and extract values for centre of 5km cells
  in.file<-paste(dir_lsratio,"invratio_",direction,"deg.tif",sep="")
  print(in.file)
  inv.lsratio<-raster(in.file) # load raster of lsratio

  # 1. Calculate Lref values for each 5km cell using mean L of all 100m land cells
  lref5km.r<-LalltoLref(inv.lsratio,gridmask.r)
  
  # 2. Convert lref back to 100m cell raster
  lref100.r<-disaggregate(lref5km.r,50)
  #plot(lref100.r,main="5km lsratio at 100m cells - lref100.r")
  #if (min(getValues(lref100.r),na.rm=TRUE)<=0){warning("Lref values of 0 or less in lref100.r")}
  
  lref100.r2<-mask(lref100.r,dem) # but is missing some 100m land cells where no overlapping 5km cell
  #plot(lref100.r2,main="lref100.r2")
  
  # set missing lref cells (currently NA in lref100.r2) to '0'
  lref100.r3<-overlay(lref100.r2,dem, fun=function(x,y){ifelse(is.na(x)&!is.na(y),0,x)})
  #plot(lref100.r3,main="lref100.r3")
  
  # 3. For each missing 100m cell set value to 'nearest' 100m  cell Lref
  # get xy and val of missing 100m cells
  xy<-xyFromCell(lref100.r3,1:ncell(lref100.r3))  
  vals<-getValues(lref100.r3)
  xyvals<-cbind(xy,vals)
  sel<-which(xyvals[,3]==0)
  missingxyv<-xyvals[sel,1:3]

  # assign value from nearest cell with Lref(ie cell covered by 5km grid)
  sel<-which(!is.na(xyvals[,3]) & xyvals[,3]!=0)
  lrefvals<-xyvals[sel,1:3]
  x<-missingxyv[,1]
  y<-missingxyv[,2]
  refx<-lrefvals[,1]
  refy<-lrefvals[,2]
  
  sqdif<-function(x,refx){return(abs(x-refx)^2)}
  
  #sel<-dim(distance)[2])

  for (i in 1: length(x)){
    distance<-sqrt(sapply(x[i],sqdif,refx) + sapply(y[i],sqdif,refy))
    nearest<-apply(distance,2,min)
    sel<-which(distance==min(distance))
    missingxyv[,3]<-lrefvals[sel,3]
  }
 
  
  
  
  
  
  dist<-function(x,y,refx,refy){
    dist<-sqrt( abs(x-refx)^2 + abs(y-refy)^2 ) 
    return(dist)
  }  
    
result<-matrix(ncol=length(x),nrow=length(refx)) 
d<-vapply()





    sel<-which(dist==min(dist))
  if (length(sel)==0){warning(paste("Warning - no min found in nearLref: x=",x," y=",y," sel=",sel," min(dist)=",min(dist),sep=""))}
  if (length(sel)>1){Lref<-mean(lrefvals[sel,3],na.rm=TRUE)
  } else {Lref<-lrefvals[sel,3] }
  missingxyv[i,3]<-Lref
  }# end for loop
  


  # 4. Update missing values (0) in lref100 to missingxyv values  
  # missingxyv[,3]<-2
  lref100.r4<-rasterize(missingxyv[,1:2], lref100.r3, missingxyv[,3], fun=max, update=TRUE)
  plot(lref100.r4)
  # test - problematic area where angle=270? - OK
  #e<-extent(240000,260000,130000,150000)
  #plot(crop(lref100.r4,e))

} # end angle loop'
  
  
  #####################################################################
  #CUT OUTS
  #####################################################################
  # then take the raster value with lowest distance to point AND non-NA value in the raster
  sampled = apply(X = xy, MARGIN = 1,
            FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
  
  
  
  Lref<-rep(NA,ncell(dem))
  for (i in 1:length(Lref)){
    if (!is.na(extract(dem,i))){dist<-sqrt( abs(x-lref5km[,1])^2 + abs(y-lref5km[,2])^2 ) 
                                sel<-which(dist==min(dist))
                                if (length(sel)==0){warning(paste("Warning - no min found in nearLref: x=",x," y=",y," sel=",sel," min(dist)=",min(dist),sep=""))}
                                if (length(sel)>1){Lref[i]<-mean(lref5km[sel,3],na.rm=TRUE)
                                } else {Lref[i]<-lref5km[sel,3] }

    }else{Lref[i]<-NA}
  }
  
  
  
  nearLref<-function(x,y,lref)
  { 
    dist<-sqrt( abs(x-lref[,1])^2 + abs(y-lref[,2])^2 ) 
    sel<-which(dist==min(dist))
    if (length(sel)==0){warning(paste("Warning - no min found in nearLref: x=",x," y=",y," sel=",sel," min(dist)=",min(dist),sep=""))}
    if (length(sel)>1){Lref<-mean(lref[sel,3],na.rm=TRUE)
    } else {Lref<-lref5km[sel,3] }
    #print(Lref)
    return(Lref)
  }
  
  Lref<-rep(NA,ncell(dem))
  for (i in 1:length(Lref)){
    if (!is.na(extract(dem,i))){Lref[i]<-nearLref(xFromCell(dem,i),yFromCell(dem,i),lref5km.r)
   }else{Lref[i]<-NA}
  }
  
  
  
  # for each mising cell set Lref value to that of nearest  valid 100m 
  #cell100xy<-xyFromCell(dem,1:ncell(dem)) # more selective eg not NA to reduce  loop
  #sel<-which(!is.na(getValues(dem)))
  #cell100xy<-cell100xy[sel,1:2] # = cell numbers not NA
  #refcell<-rep(NA,dim(cell100xy)[1])
  #Lref<-rep(NA,dim(cell100xy)[1])
  
  
  #Lref<-ifelse(!is.na(getValues(dem)), nearLref(xFromCell(dem),yFromCell(dem),lref5km.r)  , NA )
  
  for (i in 1:length(Lref)){
    dist<-sqrt( abs(cell100xy[i,1]-lref5km[,1])^2 + abs(cell100xy[i,2]-lref5km[,2])^2 ) 
    sel<-which(dist==min(dist))
    if (length(sel)>1){rpts<-rpts+1}
    refcell[i]<-lref5km[sel,3]
    #print(refcell)
  }
  print(rpts)
  
  lref5km.r<-setValues(gridmask.r,lref5km[,3]) # create 5km cell raster
  plot(lref5km.r,main="5km lsratio cells")
   
  # convert lref back to 100m cell raster
  lref100.r<-disaggregate(lref5km.r,50)
  plot(lref100.r,main="5km lsratio at 100m cells")
  lref100.r2<-mask(lref100.r,dem)
  lref100.r2<-mask(lref100.r,dem, inverse=TRUE,updatevalue=lref100.r)
  
  ldif.r<-lref100.r-inv.lsratio
  plot(ldif.r,main="Lref-Lcell")
  
} # end for each wind direction


Lref<-resample()
Lref<-mask(Lref,hrtemp) # set cells without historic temp values to NA 

# BY BLOCKS:
# NB dataset too big to work out inverse land -sea ratio in entirity, so does for each 10km square in turn (selected a
numcols<-( ((xmax(lref100.r3)-xmin(lref100.r3)) /50000 )-1) 
numrows<-( ((ymax(lref100.r3)-ymin(lref100.r3)) /50000 )-1)
for (b.col in 0:numcols ){
  for (b.row in 0:numrows ){
    xmn<-xmin(lref100.r3)+(b.col*50000)-50000
    xmx<-xmn+50000
    ymn<-(b.row*50000)+ymin(lref100.r3)-50000
    ymx<-ymn+5000
    e.cblock<-extent(c(xmn,xmx,ymn,ymx))
    centre.block<-crop(lref100.r3,e.cblock)
    
    
    get values of r block 
    correct<-function(value,constant){
      newval<-value*mean(constant)
      print(constant)
      return(newval)
    }
    correctL<-function(lcell,lrefgrid){
      # assign nearby cell Lref value to missing value cell
      dist<-sqrt( abs(x-nearby[,1])^2 + abs(y-nearby[,2])^2 ) 
      sel<-which(dist==min(dist))
      if (length(sel)==0){warning(paste("Warning - no min found in nearLref: x=",x," y=",y," sel=",sel," min(dist)=",min(dist),sep=""))}
      if (length(sel)>1){Lref<-mean(nearby[sel,3],na.rm=TRUE)
      } else {Lref<-lrefvals[sel,3] }
      missingxyv[i,3]<-Lref
    }
    
    
    
    v.block<-getValues(c.block)
    correct.vals<-ifelse(v.block<5,correct(v.block,c.block),v.block)
    print(correct.vals)
    
    
    
    
    xy<-xyFromCell(r.block,1:ncell(r.block))  
    vals<-getValues(r.block)
    xyvals<-cbind(xy,vals)
    sel<-which(xyvals[,3]==0)
    missingxyv<-xyvals[sel,1:3]
      Yes - calculate nearest neighbour 
      No - next block
  
    
    
for (i in 1:dim(missingxyv)[1]){
  if (i%%10==0){print(i)}
  x<-missingxyv[i,1]
  y<-missingxyv[i,2]
  xcol<-(x-min(lref100.r3))/res(lref100.r3)
  yrow<-(y-min(lref100.r3))/res(lref100.r3)
  for (step in 1:50){ # check up to 50 cells away
    # get values for whole rows
    test.x<-x-step:x+step
  }    

}
    
    
    

for (i in 1:dim(missingxyv)[1]){
  if (i%%10==0){print(i)}
  x<-missingxyv[i,1]
  y<-missingxyv[i,2]
  # define nearby valid cells within 5km radius
  xmn<-max(x-5000,xmin(lref100.r3))
  xmax<-min(x+5000,xmax(lref100.r3))
  ymn<-max(y-5000,xmin(lref100.r3))
  ymx<-min(y+5000,xmax(lref100.r3))
  e.nearby<-extent(xmn,xmx,ymn,ymx)
  nearby.r<-crop(lref100.r3,e.nearby)
  nearxy<-xyFromCell(nearby.r,1:ncell(lref100.r3))  
  nearvals<-getValues(nearby.r)
  nearby<-cbind(nearxy,nearvals)
  sel<-which(!is.na(nearby[,3]) & nearby[,3]!=0)
  nearby<-nearby[sel,1:3]
  # assign nearby cell Lref value to missing value cell
  dist<-sqrt( abs(x-nearby[,1])^2 + abs(y-nearby[,2])^2 ) 
  sel<-which(dist==min(dist))
  if (length(sel)==0){warning(paste("Warning - no min found in nearLref: x=",x," y=",y," sel=",sel," min(dist)=",min(dist),sep=""))}
  if (length(sel)>1){Lref<-mean(nearby[sel,3],na.rm=TRUE)
  } else {Lref<-lrefvals[sel,3] }
  missingxyv[i,3]<-Lref
}# end for loop

WORKS!!
sqdif<-function(x,refx){return(abs(x-refx)^2)}
result<-sqrt(sapply(x,sqdif,refx) + sapply(y,sqdif,refy))



# Create ever larger fileter and apply with focal
# for use where sea=0, missing = NA
correct.Lval<-function(r){
rad<-200
while(length(is.na(getValues(r)))>=1){
  nearby.r<-raster(ncols=rad/100, nrows=rad/100, xmn=0)  
  nearby<-focalWeight(nearby.r,rad,type=c("circle"))
  new<-focal(r,nearby,NAonly=TRUE)
  print(paste("rad= ",rad," filter: ",dim(nearby),sep=""))
  rad<-rad+100
  if (rad>1000) break
  }
if (length(is.na(getValues(r)))>=1){print("still missing cells")}
return(r)
}

# c. Calc SeaST-5kmref difference 