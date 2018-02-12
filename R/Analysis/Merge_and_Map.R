print("Calling set up...")
source("/home/ISAD/jm622/rscripts/setup_v4_carson.R") # loads & runs setup file
#library(insol)

######################################################
# Merge and print temperature for defined times
#cells<-c(718,719,745,746)
cells<-c(718,719,720,744,745,746,769,770,771)
#cells<-c(898,899,900,911,912,913,921,922,923)
numcells<-length(cells)
# Set year and load year's data'
year<-1983; 
month<-6  ; day<-25  
DOY<-daydoy(year,month,day)
datelabel<-paste(year,month,day,sep="_")
print(datelabel)
#######################################################
par(mfrow=c(3,4))

# 1. Print hourly Temperature maps for 24hr 
# save block of 24hr data
map.24h<-brick()

for (hr in 0:23){
blocks<-vector("list",numcells)
i<-1
for (ukcpcell in cells){
  # Define dem.block
  print(paste("UKCP cell= ", ukcpcell,sep=""))
  print("Get Cell Coordinates")
  gridmask.r<-land5km.r # land5km.r defined in setup
  vals<-getValues(gridmask.r)
  xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
  sel<-which(vals==1)
  landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell
  print(dim(landcells))
  x<-landcells[ukcpcell,1]
  y<-landcells[ukcpcell,2]
  e.block<-extent(x-2500,x+2500,y-2500,y+2500)
  dem.block<-crop(demuk,e.block) 
  
  # load this year's temperature data at 100m and hourly resolution
  infile<-paste(dir_finalt,"block-",ukcpcell,"-",year,".R",sep="")
  print(infile)
  load(infile) # = tmodel.r
  numdays<-dim(tmodel.r)[3]/24
  print(numdays)
  
  # Extract 24h of temperature data
  #start<-1+(doy-1)*24
  #end<-doy*24
  #print (paste(start,"  ",end))
  #t.24h<-tmodel.r[,,start:end]
  
  #print (paste(start,"  ",end))
  t.hr<-tmodel.r[,,(hr+(doy-1)*24)]
  
  blocks[[i]]<-raster(as.matrix(t.hr),template=dem.block)
  i<-i+1
}

### SAVE and PLOT merged blocks
map.24h<-addLayer(map.24h,do.call(merge, blocks))
plot(map.r,main=paste(datelabel," ",hr,"h",sep=""))
} # end 24 hr

writeRaster(map.24h,file=paste(dir_tmaps,"tmap24h_camel_",datelabel,".tif",sep=""))


#######################################################
# 3. Plot elements of model or input data
#year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day ;# values for t

for (hr in 0:23){
  
  blocks<-vector("list",numcells)
  i<-1
  for (ukcpcell in cells){
    # Define dem.block
    print(paste("UKCP cell= ", ukcpcell,sep=""))
    print("Get Cell Coordinates")
    gridmask.r<-land5km.r # land5km.r defined in setup
    vals<-getValues(gridmask.r)
    xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
    sel<-which(vals==1)
    landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell
    print(dim(landcells))
    x<-landcells[ukcpcell,1]
    y<-landcells[ukcpcell,2]
    e.block<-extent(x-2500,x+2500,y-2500,y+2500)
    dem.block<-crop(demuk,e.block) 
    
    # load and extract RHref data
    rh.filein<-paste(dir_rh5km,"RH_100m_",year,"_",month,"_",day,".R",sep="")
    load(rh.filein) # rh.day
    rhref.block<-crop(raster(rh.day[,,hr+1],template=dem),dem.block)
    rhref.block<-crop(raster(rh.day[,,hr+1],template=dem),dem.block)
    
    blocks[[i]]<-raster(as.matrix(rhref.block),template=dem.block)
    i<-i+1
  }
  
  ### PLOT merged blocks
  map.r<-do.call(merge, blocks)
  titletext<-paste("RHref ",hr,"hr on DOY= ",doy+90,sep="")
  plot(map.r,main=titletext)
  
} # end 24 hr










# Load and reproject CAL data to OS (downscaling by resampling takes place in hourly interval)
cal.filein<-paste(dir_calday,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
print(paste("CAL file in: ",cal.filein,sep=""))
cal.day<-brick(cal.filein)
cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")
cal.block<-tps.resample(crop(raster(cal.day,layer=hr+1),dem.buffer),dem.block)




#######################################################




    blocks<-vector("list",9)
    i<-1
    for (n in(1:length(cells))) { 
      infile<-paste(stats[x],cells[n],"-",year,".tif",sep="")
      print(infile)
      block.r<-raster(infile)
      #plot(block.r,main=cells[n])
      #if(ukcpcell==745) block.r<-calc(block.r,fun=function(x){x=x+10})
      #if(ukcpcell==744) block.r<-calc(block.r,fun=function(x){x=x-10})
      blocks[[i]]<-block.r
      i<-i+1
    }
    map.r<- do.call(merge, blocks)
    
    
    # Print
    brks <- seq(0, 14, by=1) 
    plot(map.r, col=rev(heat.colors(14)),axes=FALSE,box=FALSE)
    plot(map.r, col=rev(terrain.colors(length(seq(500, 1200, by = 50))-1)),axes=FALSE,box=FALSE, breaks=seq(500, 1200, by = 50),legend=F )
    