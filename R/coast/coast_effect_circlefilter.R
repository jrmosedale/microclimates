
dir_out<-"~/Documents/Exeter/Data2015/Coast_effect/"
#dir_out<-"C:/Data2015/Coast_effect"

########################################################################################
#dem<-raster("C:/Data2015/DEM100/demoriginal.tif")
demuk<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif", crs=("+init=epsg:27700"))
# Create sw dem of interest (no buffer)
#e.dem<-extent(c(70000,420000,0,180000 )) # includes scilly isles
e.dem<-extent(c( 120000,420000,0,180000 )) # excludes scilly isles
dem<-crop(demuk,e.dem)

# Create buffered UK dem
demukb<-extend(demuk,c(200,200),value=NA) # extend 20km of "sea" around dem
# 1.Create buffered south west dem raster required
# Sea cells = NA, land cells = elevation
buffer<-20000 # include 20km buffer area around zone of interest
#e.dem<-extent(c((70000-buffer),(420000+buffer),(0-buffer),(180000+buffer) )) # includes scilly isles
e.dem<-extent(c( (120000-buffer),(420000+buffer),(0-buffer),(180000+buffer) )) # excludes scilly isles
demsw<-crop(demukb,e.dem)
plot(demsw,main="DEM-sw")
e.dem <-extent(demsw)
# PROBLEM: ideally require land outline/raster for whole UK/N.France 

# 2. Create simplified land/sea raster - USED to calculate coastal index.
# Cell values: Sea=1, Land=0 no NA
dem.vals<-getValues(demsw)
dem.vals<-ifelse(is.na(dem.vals),1,0)
landsea.r<-demsw
landsea.r<-setValues(landsea.r,dem.vals)
plot(landsea.r,main="landsea.r")

# 3. Create raster of cells within 20km of coast from orig dem - those for which coastal index required
# a Create simple land/sea raster. Cell values: LAND cells = NA, Sea cells=0
dem.vals<-getValues(demsw)
dem.vals<-ifelse(is.na(dem.vals),-999,NA)
seaonly.r<-setValues(demsw,dem.vals)
plot(seaonly.r)

# b Calculate coastal distance raster (val=metres from sea cell (!=NA))
coastdist.r<-distance(seaonly.r) # this takes a few mins - sea cells (prev NA) recorded as 0
#vals<-getValues(coastdist.r)
#vals<-ifelse(vals==0,NA,vals)
#coastdist.r<-setValues(coastdist.r,vals)
plot(coastdist.r)
# write(sea.dist,file=paste(dir_out,"disttosea.r",sep=""))

# c Create raster of cells within 20km of coast - set other cells to NA 
coastal.r<-coastdist.r
coastal.r[coastal.r==0]<-NA
m <- c(0 , 20000, 0,  20000, 100000, 1)
rclas <- matrix(m, ncol=3, byrow=TRUE)
coastal.r <- reclassify(coastal.r,rclas)
plot(coastal.r)
# crop to only those cells within original dem (not buffered region for which coast effects not calculated)

########################################################################################

# Create fileter for using with focal function
# Create matrix to hold circlulaar area of rad r

circle.cells<-function(res, radius) {  # circular filter weighted according to distance from centre, direction etc
  m<-matrix(0,nrow=1+(2*radius/res),ncol=1+(2*radius/res))
  dimnames(m)[[1]]<-seq(-radius,radius,by=res)
  dimnames(m)[[2]]<-seq(-radius,radius,by=res)
  
  for (row in 1:nrow(m)){
    for (col in 1:ncol(m)){
      dist<-sqrt((as.numeric(dimnames(m)[[1]])[row])^2 +
                   (as.numeric(dimnames(m)[[1]])[col])^2)
      if(dist<=radius) {m[row,col]<-1} # change here method of calculating index
    }
  }
  circle.pts<-sum(m)
  return(circle.pts)
}
########################################################################################

fill.filter<-function(m) {  # circular filter weighted according to distance from centre, direction etc
  area<-pi*radius^2 # find function to calculate number of cells within area?
    for (row in 1:nrow(m)){
        for (col in 1:ncol(m)){
            dist<-sqrt((as.numeric(dimnames(m)[[1]])[row])^2 +
                         (as.numeric(dimnames(m)[[1]])[col])^2)
            if(dist<=radius) {m[row,col]<-(1-(dist/radius))/area} # change here method of calculating index
        }
    }
    return(m)
}

# Alternative method of creating filter
r<-raster(ncols=10,nrows=10,xmn=0,res=res)
filter2<-focalWeight(r,radius,"circle")



radius<-5000
res<-100

filter<-matrix(0,nrow=1+(2*radius/res),ncol=1+(2*radius/res))
dimnames(filter)[[1]]<-seq(-radius,radius,by=res)
dimnames(filter)[[2]]<-seq(-radius,radius,by=res)

filter<-fill.filter(filter) # circular filter =1 at centre, 0 at edge of radius
filter.r<-raster(filter)
plot(filter.r)

#####################################################################

# Apply filter to raster where sea=1, land=0
e<-extent(120000,150000,0,64000)
tmp.r<-crop(landsea.r,e)
plot(tmp.r)

new.r<-focal(tmp.r,filter) # where sea=1, land=0 
new.r<-focal(landsea.r,filter)
plot(new.r); new.r

# set sea cells back to NA and inland cells to 1 (currently 0)
final.r<-mask(new.r, demsw)
plot(final.r)

