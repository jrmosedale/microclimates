# Calculate grid cell coordinates from dem using set grid cell size
resolution <-100 # of dem original
ukdem.file<-paste(dir_dem,"demoriginal.tif",sep="")
demuk<-raster(ukdem.file)

# Core area of interest boundary
e.dem<-extent(c( 120000,420000,0,180000 ))
gridcell<- 20000 # in metres
e.use<-extent(c( 120000,420000,0,180000 ))

# Calculate if land cell at 100m

land100m.r<-calc(dem,function(x) ifelse(is.na(x),NA,1))

land100m.r<-flip(land100m.r,direction='x')
landgridcell.r<-aggregate(land100m.r,fact=(gridcell/100),fun=sum,expand=TRUE)
landgridcell.r<-flip(landgridcell.r,direction='x')

landgridcell.r<-calc(landgridcell.r,function(x) ifelse(is.na(x),NA,1))
plot(landgridcell.r)

vals<-getValues(landgridcell.r)
xy<-xyFromCell(landgridcell.r,1:ncell(landgridcell.r))
sel<-which(vals==1)

cells.df<-as.data.frame(xy[sel,1:2])   # = coordinates for middle of each ukcp09 cell
cells.df$xmin<-cells.df$x-(gridcell/2)
cells.df$xmax<-cells.df$x+(gridcell/2)
cells.df$ymin<-cells.df$y-(gridcell/2)
cells.df$ymax<-cells.df$y+(gridcell/2)


# OR SIMPLY EXTEND RASTER TO SOUTH ASSUMING ALL SEA CELLS (NA)

# Extend raster r to extent e
add_seacells<-function(r,e){
  newdem<-extend(r,e)
}


newdem<-extend(dem,extent(c( 120000,420000,0,180000 )))
land100m.r<-calc(newdem,function(x) ifelse(is.na(x),NA,1))

# Calculate number of 100m land cells in gridcell 
landgridcell.r<-aggregate(land100m.r,fact=(gridcell/100),fun=sum,expand=TRUE)
landgridcell.r<-calc(landgridcell.r,function(x) ifelse(is.na(x),NA,1))
plot(landgridcell.r)

vals<-getValues(landgridcell.r)
xy<-xyFromCell(landgridcell.r,1:ncell(landgridcell.r))
sel<-which(vals==1)

cells.df<-as.data.frame(xy[sel,1:2])   # = coordinates for middle of each ukcp09 cell
cells.df$xmin<-cells.df$x-(gridcell/2)
cells.df$xmax<-cells.df$x+(gridcell/2)
cells.df$ymin<-cells.df$y-(gridcell/2)
cells.df$ymax<-cells.df$y+(gridcell/2)

# sort so 1st gridcell in SW corner
cells.df<-cells.df[order(cells.df$x, cells.df$y),]

# TEST
for (n in 1:nrow(cells.df)){
  plot(crop(newdem,extent(cells.df$xmin[n],cells.df$xmax[n],cells.df$ymin[n],cells.df$ymax[n])),main=paste("20km cell ",n))
}