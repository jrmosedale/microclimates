args <-commandArgs(trailingOnly = TRUE)
print(args)
day <- as.integer(args[1])
month<-as.integer(args[2])
year<-as.integer(args[3] )
print(paste(day,"/",month,"/",year,sep=""))

print("Calling set up...")
source("/home/ISAD/jm622/rscripts/setup_v4_carson.R") # loads & runs setup file

######################################################
# Merge and print temperature for defined times
#cells<-c(718,719,720)
cells<-c(718,719,720,744,745,746,769,770,771)
#cells<-c(898,899,900,911,912,913,921,922,923)
numcells<-length(cells)
# Set year and load year's data'
#year<-2012; month<-6  ; day<-24  
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
    dem.block<-crop(dem,e.block) 
    
    # load this year's temperature data at 100m and hourly resolution
    infile<-paste(dir_finalt,"block-",ukcpcell,"-",year,".R",sep="")
    print(infile)
    load(infile) # = tmodel.r
    numdays<-dim(tmodel.r)[3]/24
    print(numdays)
    
    # Extract 24h of temperature data
    # THIS DEPENDS ON AMOUNT OF DATA IN TMODEL.R
    #t.hr<-tmodel.r[,,(1+hr+((DOY-1)*24))]
    t.hr<-tmodel.r[,,(hr+1)]
    blocks[[i]]<-raster(as.matrix(t.hr),template=dem.block)
    i<-i+1
  }
  ### SAVE and PLOT merged blocks
  map.r<-do.call(merge, blocks)
  plot(map.r,main=paste(datelabel," ",hr,"h",sep=""))
  map.24h<-addLayer(map.24h,map.r)
  
} # end 24 hr

writeRaster(map.24h,file=paste(dir_tmaps,"tmap24h_test_",datelabel,".tif",sep=""),overwrite=TRUE)
