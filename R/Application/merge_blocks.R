# Re-load and stitch together block results if required using Mosaic etc 
print("Calling set up...")
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file

cells<-c(925,926,927,928,929,930,931,932,934,935)
jd<-JDdmy(6,7,1992)

print("Get Cell Coordinates")
gridmask.r<-land5km.r # land5km.r defined in setup
vals<-getValues(gridmask.r)
xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
sel<-which(vals==1)
landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell


# 1. Model Temperature
par(mfrow=c(2,2))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map.s<-stack()

for (hr in 0:23)  { 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  map.r<-raster()
  for (n in(1:length(cells))) { # create day stack of hourly data for whole of area
    file.24h<-paste(dir_finalt,"block-",cells[n],"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(file.24h)
    t.r<-raster(file.24h,band=hr+1)
    projection(t.r)<-CRS("+init=epsg:27700")
    #plot(t.r,main=ukcpcell[n])
    if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  }
  plot(map.r,main=paste("Model T: ",timelabel))
  map.s<-stack(map.s,map.r)
}
#writeRaster(map.s,file=paste(dir_tmaps,"Lizardmap_",day,"_",month,"_",year,".tif",sep=""),overwrite=TRUE,format="GTiff")



# 2. Ref Temperature
par(mfrow=c(1,1))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map2.s<-stack()

file.24h<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_100m.r", sep="")
print(file.24h)
load(file.24h)

for (hr in 0:23)  { 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  map.r<-raster()
  for (n in(1:length(cells))) { # create day stack of hourly data for whole of area
    print("Get Cell Coordinates")

    # Crop dem to fit cell and buffered cell
    x<-landcells[cells[n],1]
    y<-landcells[cells[n],2]
    e.block<-extent((x-2500),(x+2500),(y-2500),(y+2500))
    #dem.block<-crop(demuk,e.block) ; print (dem.block)
    t.r<-crop(raster(t100m.day[,,hr+1],template=dem),e.block)
    projection(t.r)<-CRS("+init=epsg:27700")
    #plot(t.r,main=ukcpcell[n])
    if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  }
  plot(map.r,main=paste("Ref T: ",timelabel))
       map2.s<-stack(map2.s,map.r)
}

# 3. Rel Humidity
par(mfrow=c(1,1))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map2.s<-stack()

file.24h<-paste(dir_rh5km,"RH_100m_",DMYjd(jd)$year,"_",DMYjd(jd)$month,"_",DMYjd(jd)$day,".R",sep="")
print(file.24h)
load(file.24h)

for (hr in 0:23)  { 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  map.r<-raster()
  for (n in(1:length(cells))) { # create day stack of hourly data for whole of area
    print("Get Cell Coordinates")
    
    # Crop dem to fit cell and buffered cell
    x<-landcells[cells[n],1]
    y<-landcells[cells[n],2]
    e.block<-extent((x-2500),(x+2500),(y-2500),(y+2500))
    #dem.block<-crop(demuk,e.block) ; print (dem.block)
    t.r<-crop(raster(rh.day[,,hr+1],template=dem),e.block)
    projection(t.r)<-CRS("+init=epsg:27700")
    #plot(t.r,main=ukcpcell[n])
    if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  }
  plot(map.r,main=paste("Rel Humidity: ",timelabel))
  map2.s<-stack(map2.s,map.r)
}


# 4. DEM
par(mfrow=c(1,1))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map4.s<-stack()
map.r<-raster()
for (n in(1:length(cells))) { # create day stack of hourly data for whole of area
  print("Get Cell Coordinates")
  
  # Crop dem to fit cell and buffered cell
  x<-landcells[cells[n],1]
  y<-landcells[cells[n],2]
  e.block<-extent((x-2500),(x+2500),(y-2500),(y+2500))
 t.r<-crop(demuk,e.block) ; print (dem.block)
  projection(t.r)<-CRS("+init=epsg:27700")
  #plot(t.r,main=ukcpcell[n])
  if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
}
plot(map.r,main="DEM")
map2.s<-stack(map2.s,map.r)

cells<-c(c(925,926,927,928,929,930,931,932,934,935))
# A. RadEffect
par(mfrow=c(2,2))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map.s<-stack()

for (hr in 0:23)  { 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  map.r<-raster()
  for (n in(1:length(cells))) { # create day stack of hourly data for whole of area
    file.24h<-paste(dir_finalt,"radef-",cells[n],"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(file.24h)
    t.r<-raster(file.24h,band=hr+1)
    projection(t.r)<-CRS("+init=epsg:27700")
    #plot(t.r,main=ukcpcell[n])
    if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  }
  plot(map.r,main=paste("RadEff: ",timelabel))
  map.s<-stack(map.s,map.r)
}

#B. LatEffect
par(mfrow=c(2,2))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map.s<-stack()

for (hr in 0:23)  { 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  map.r<-raster()
  for (n in(1:length(cells))) { # create day stack of hourly data for whole of area
    file.24h<-paste(dir_finalt,"latef-",cells[n],"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(file.24h)
    t.r<-raster(file.24h,band=hr+1)
    projection(t.r)<-CRS("+init=epsg:27700")
    #plot(t.r,main=ukcpcell[n])
    if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  }
  plot(map.r,main=paste("LatEff: ",timelabel))
  map.s<-stack(map.s,map.r)
}

#C. Coast effect
par(mfrow=c(2,2))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map.s<-stack()

for (hr in 0:23)  { 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  map.r<-raster()
  for (n in(1:length(cells))) { # create day stack of hourly data for whole of area
    file.24h<-paste(dir_finalt,"coastef-",cells[n],"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(file.24h)
    t.r<-raster(file.24h,band=hr+1)
    projection(t.r)<-CRS("+init=epsg:27700")
    #plot(t.r,main=ukcpcell[n])
    if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  }
  plot(map.r,main=paste("Coast Eff: ",timelabel))
  map.s<-stack(map.s,map.r)
}


#D. Anom effect
par(mfrow=c(2,2))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map.s<-stack()

for (hr in 0:23)  { 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  map.r<-raster()
  for (n in(1:length(cells))) { # create day stack of hourly data for whole of area
    file.24h<-paste(dir_finalt,"anom-",cells[n],"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(file.24h)
    t.r<-raster(file.24h,band=hr+1)
    projection(t.r)<-CRS("+init=epsg:27700")
    #plot(t.r,main=ukcpcell[n])
    if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  }
  plot(map.r,main=paste("Anom: ",timelabel))
  map.s<-stack(map.s,map.r)
}




# Plot elevation effect over 24hrs for Lizard - Confirms grid effect on downscaled temperatue 
hr<-7
ukcpcell<-c(916:919,924:932,934:935)

# load stack of 5km temperature data 
t5km.s<-stack()
filein<-paste(dir_hrtemp,"HrTemp_2000-05-10.r", sep="")
load(filein)
print(filein)
for (n in 1:24){
  t5km.r<-raster(t5km.day[,,n],template=grid5km.r)
  t5km.s<-addLayer(t5km.s,t5km.r) 
}

# 1. Calculate map of interpolated tref values using tps
ptm <- proc.time()[1:3] 
tref.map<-raster()
for (n in 1:length(ukcpcell)){
  x<-landcells[ukcpcell[n],1]
  y<-landcells[ukcpcell[n],2]
  e.block<-extent(x-2500,x+2500,y-2500,y+2500)
  dem.block<-crop(demuk,e.block)
  #tref.r<-ref5km.to.block100m(dem.block,t5km.s[[1+hr]])
  tref.r<-tps.resample(t5km.s[[n]],dem.block) 
  if (n==1) tref.map<-tref.r else tref.map<-merge(tref.map,tref.r)
}
print(proc.time()[1:3]-ptm)
plot(tref.map)



#Calculate interpolated tref + elev effect and plot
elevdif.r<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # created by elevation.dif.map function
effect.r<- (elevdif.r/100)*0.66
effect.r<-crop(effect.r,tref.map)
plot(effect.r,main="elevation anomaly")
t.r<-tref.map+effect.r
plot(t.r,main=paste("Tref+elev effect at hr ",hr,sep=""))
  


################################
# 2. CHALLANGE FID ALTERATIVE TO BELOW _ USING tps to create smooth surface of ref temperatures 
# Assign temperature values from 5km cells to 100m cells in block
# Function returns 100m raster from 5km raster for block
# Compare percentland_map_function & elevdif_map_function

tps5km.to.block100m<-function(dem.buffer,t5km.r)
{
  t5km.buffer<-tps.resample(t5km.r,dem.buffer)
  ref.block<-mask(ref.block,dem.block)
  #plot(ref.block,main="ref block")
  return(ref.block)
} # end function 


ukcpcell<-c(916:919,924:932,934:935)
jd<-JDdmy(10,5,2000)
projection(map.r)<-CRS("+init=epsg:27700")

par(mfrow=c(1,1))
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day
map.s<-stack()

gridmask.r<-land5km.r # land5km.r defined in setup
vals<-getValues(gridmask.r)
xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
sel<-which(vals==1)
landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell
print(dim(landcells))
hr<- 1

ptm <- proc.time()[1:3] 
map.r<-raster()
for (n in(1:length(ukcpcell)))
{ 
  x<-landcells[ukcpcell[n],1]
  y<-landcells[ukcpcell[n],2]
  e.block<-extent(x-2500,x+2500,y-2500,y+2500)
  dem.block<-crop(demuk,e.block)# create day stack of hourly data for whole of area
  t.r<-tps.resample(t5km.r,dem.block)
  if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  
}
print(proc.time()[1:3]-ptm)

plot(map.r,main=hr)


# 3. Calculate using resample bilinear
ptm <- proc.time()[1:3] 
map.r<-raster()
for (n in(1:length(ukcpcell)))
{ 
  x<-landcells[ukcpcell[n],1]
  y<-landcells[ukcpcell[n],2]
  e.block<-extent(x-2500,x+2500,y-2500,y+2500)
  dem.block<-crop(demuk,e.block)# create day stack of hourly data for whole of area
  t.r<-resample(t5km.r,dem.block,method='bilinear')
  if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  
}
map.r<-mask(map.r,crop(dem,map.r))
print(proc.time()[1:3]-ptm)
plot(map.r,main="bilinear")

# 4. Calculate using resample bilinear
library(PEIP)
ptm <- proc.time()[1:3] 
map.r<-raster()

for (n in(1:length(ukcpcell)))
{ 
  x<-landcells[ukcpcell[n],1]
  y<-landcells[ukcpcell[n],2]
  e.block<-extent(x-2500,x+2500,y-2500,y+2500)
  dem.block<-crop(demuk,e.block); plot(dem.block)
 # incells<-Which(!is.na(t5km.r), cells = TRUE) 
  incells<-c(1:ncell(t5km.r))
  xin<-xyFromCell(t5km.r,incells)[,1]
  yin<-xyFromCell(t5km.r,incells)[,2]
  zin<-getValues(t5km.r,format="matrix")
  outcells<-Which(!is.na(dem.block), cells=TRUE)
  xout<-xyFromCell(dem.block,outcells)[,1]
  yout<-xyFromCell(dem.block,outcells)[,2]
  indata<-list(xin,yin,zin)
  outdata<-list(xout,yout)
  t.r<-interp.surface.grid(indata,outdata )
  
  #mat<-extract(t5km.r,cellnnumbers=cellnumbers,format="matrix")
#  t.r<-interp2grid(mat,xout,yout,xin,yin,type=1)
  
  
  if (n==1) map.r<-t.r else map.r<-merge(map.r,t.r)
  
}
map.r<-mask(map.r,crop(dem,map.r))
print(proc.time()[1:3]-ptm)
plot(map.r,main="bilinear")


# TRy using alternative fields functions for interpolation


#### Create interpolation model and compare for every 6 hour of one day

# plot hourly time series 
#  day<-
par(mfrow=c(3,2))

tps<-list();
tps.hr<-list()
for (hr in seq(1,24,6)){
  print(hr)
  hr.new<-hr-1
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  t5km.r<-t5km.s[[hr]]; plot (t5km.r,main="hr")
  proc.time()->ptm
  xy <- data.frame(xyFromCell(t5km.r,1:ncell(t5km.r)))
  v <- getValues(t5km.r)
  tps.new <- Tps(xy, v) # fit tps model
  print (proc.time()-ptm)
  #plot(tps.new)
  #surface(tps.new)
  #print(summary(tps.new))
  #tps<-c(tps,tps.tab)
  #tps.hr<-c(tps.hr,hr.new)
  t.block<-interpolate(dem.block,tps.new)
  print (proc.time()-ptm)
  
}
model10052000<-tps


# Compare interpolation using models from different times of day 
output<-raster(dem.block)
ptm<-proc.time()
t100.r<-interpolate(t5km.s[[1]],tps.new)
print(proc.time()-ptm)
