par(mfrow=c(2,3))
cells<-c(925,926,927,928,929,930,931,932,934,935)
day<-DMYjd(jd)$day;month<-DMYjd(jd)$month;year<-DMYjd(jd)$year


# 1. TESTING RAD LWR AND CAL
# Load day files
filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
dnr.24h<-stack(filein1) ; projection(dnr.24h)<-CRS("+init=epsg:4326")
sis.24h<-stack(filein2) ; projection(sis.24h)<-CRS("+init=epsg:4326")
dnr.24h<-projectRaster(dnr.24h,crs="+init=epsg:27700")
sis.24h<-projectRaster(sis.24h,crs="+init=epsg:27700")
# Load daily files of hourly values of temperature, RH, lwr
t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_100m.r", sep="")
load(t.filein) #t100m.day
# load 5km RH data
rh.filein<-paste(dir_rh5km,"RH_100m_",year,"_",month,"_",day,".R",sep="")
load(rh.filein) # rh.day
# Load and reproject CAL data to OS (downscaling by resampling takes place in hourly interval)
cal.filein<-paste(dir_calday,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
print(paste("CAL file in: ",cal.filein,sep=""))
cal.day<-brick(cal.filein)
cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")


for (hr in 6:10)
{
tradmap<-raster()
calmap<-raster()
lwrmap<-raster()

for (n in 1: length(cells)) {
print(cells[n]) 
# Extract hourly data from day stack
dnr.r<-raster(dnr.24h,layer=hr+1); #plot(dnr.r)
sis.r<-raster(sis.24h,layer=hr+1); #plot(sis.r)

# calc dem.block/buffer
x<-landcells[cells[n],1]
y<-landcells[cells[n],2]
e.block<-extent(x-2500,x+2500,y-2500,y+2500)
dem.block<-crop(demuk,e.block) ;# plot(dem.block)
e.buffer<-extent(x-2500-buffer,x+2500+buffer,y-2500-buffer,y+2500+buffer)
dem.buffer<-crop(demuk,e.buffer)
# Slope & aspect cropped to dem.buffer for use in radiation_downscale
slope.buffer<-crop(slope.r,dem.buffer)
aspect.buffer<-crop(aspect.r,dem.buffer)

rad.results<-radiation_downscale_stack(day,month,year,hr,sis.r,dnr.r,dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=FALSE)
direct.block<-rad.results[[1]]*W.to.MJhr
diffuse.block<-rad.results[[2]]*W.to.MJhr
total.block<-rad.results[[3]]*W.to.MJhr # = total incoming shortwave radiation - IGNORING ALBEDO EFFECT

### CALCULATE LWR (effective resolution 5km) ###
tref.block<-crop(raster(t100m.day[,,hr+1],template=dem),dem.block)
rhref.block<-crop(raster(rh.day[,,hr+1],template=dem),dem.block)
cal.block<-tps.resample(crop(raster(cal.day,layer=hr+1),dem.buffer),dem.block)
lwr.block<-calc_lwr_block(jd,cal.block,rhref.block,tref.block,plotlwr=FALSE,writefile=FALSE)

#plot(total.block,main=(paste("Total Rad at ",hr,"h00",sep="")))
#plot(tref.block,main=paste("Tref ",hr))
#plot(rhref.block,main=paste("RHref ",hr))
#plot(cal.block,main=paste("CAL ",hr))
#plot(lwr.block,main=paste("LWR ",hr))


if (n==1) tradmap<-total.block else tradmap<-merge(tradmap,total.block)
if (n==1) lwrmap<-lwr.block else lwrmap<-merge(lwrmap,lwr.block)
if (n==1) calmap<-cal.block else calmap<-merge(calmap,cal.block)
} # cell
plot(tradmap,main=paste("Total Rad at ",hr))
plot(calmap,main=paste("CAL at ",hr))
plot(lwrmap,main=paste("LWR at ",hr))
}# hr


#####################################



# 2. TESTING WIND
# Load day files
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))
interval<-10
for (hr in 6:10)
{
  wstrmap<-raster()
  wdirmap<-raster()
  invwmap<-raster()

    for (n in 1: length(cells)) {
    print(cells[n]) 
  
    # calc dem.block/buffer
    x<-landcells[cells[n],1]
    y<-landcells[cells[n],2]
    e.block<-extent(x-2500,x+2500,y-2500,y+2500)
    dem.block<-crop(demuk,e.block) ;# plot(dem.block)
    e.buffer<-extent(x-2500-buffer,x+2500+buffer,y-2500-buffer,y+2500+buffer)
    dem.buffer<-crop(demuk,e.buffer)
    # Slope & aspect cropped to dem.buffer for use in radiation_downscale
    slope.buffer<-crop(slope.r,dem.buffer)
    aspect.buffer<-crop(aspect.r,dem.buffer)
    shelter.block<-block.sheltermap(dem.block,dir_shelter,interval)
    
    wind.results<-wind.tpsdownscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval,print.res=FALSE,write.files=FALSE )
    wdir.block<-wind.results[[1]]
    wstr.block<-wind.results[[2]]
    invwstr.block<-wind.results[[3]]     
 
    if (n==1) wstrmap<-wstr.block else wstrmap<-merge(wstrmap,wstr.block)
    if (n==1) wdirmap<-wdir.block else lwrmap<-merge(wdirmap,wdir.block)
    if (n==1) invwmap<-invwstr.block else invwmap<-merge(invwmap,invwstr.block)
  } # cell
  
  plot(wstrmap,main=paste("Wind Str at ",hr))
  plot(wdirmap,main=paste("Wind dir at ",hr))
  plot(invwmap,main=paste("Inv Wind Str at ",hr))
}# hr



#####################################


# 3. Test all RafEffect factors incl albedo, 1 & 2
# Load day files
filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
dnr.24h<-stack(filein1) ; projection(dnr.24h)<-CRS("+init=epsg:4326")
sis.24h<-stack(filein2) ; projection(sis.24h)<-CRS("+init=epsg:4326")
dnr.24h<-projectRaster(dnr.24h,crs="+init=epsg:27700")
sis.24h<-projectRaster(sis.24h,crs="+init=epsg:27700")
# Load daily files of hourly values of temperature, RH, lwr
t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_100m.r", sep="")
load(t.filein) #t100m.day
# load 5km RH data
rh.filein<-paste(dir_rh5km,"RH_100m_",year,"_",month,"_",day,".R",sep="")
load(rh.filein) # rh.day
# Load and reproject CAL data to OS (downscaling by resampling takes place in hourly interval)
cal.filein<-paste(dir_calday,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
print(paste("CAL file in: ",cal.filein,sep=""))
cal.day<-brick(cal.filein)
cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")

elevdif.r<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # created by elevation.dif.map function
projection(elevdif.r)<-CRS("+init=epsg:27700")

load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))
interval<-10

albedomap.r<-raster(paste(dir_albedo,"albedo_mean_dembuf.tif",sep=""))
projection(albedomap.r)<-CRS("+init=epsg:27700")

#sheltermap<-block.sheltermap(dem,dir_shelter,interval)

for (hr in 6:10)
{
  tradmap<-raster()
  calmap<-raster()
  lwrmap<-raster()
  wstrmap<-raster()
  wdirmap<-raster()
  invwmap<-raster()
  RADEFmap<-raster()
  ELEVEF<-raster()
  tmap<-raster()
  
  for (n in 1: length(cells)) {
    print(cells[n]) 
    # Extract hourly data from day stack
    dnr.r<-raster(dnr.24h,layer=hr+1); #plot(dnr.r)
    sis.r<-raster(sis.24h,layer=hr+1); #plot(sis.r)
    
    # calc dem.block/buffer
    x<-landcells[cells[n],1]
    y<-landcells[cells[n],2]
    e.block<-extent(x-2500,x+2500,y-2500,y+2500)
    dem.block<-crop(demuk,e.block) ;# plot(dem.block)
    e.buffer<-extent(x-2500-buffer,x+2500+buffer,y-2500-buffer,y+2500+buffer)
    dem.buffer<-crop(demuk,e.buffer)
    # Slope & aspect cropped to dem.buffer for use in radiation_downscale
    slope.buffer<-crop(slope.r,dem.buffer)
    aspect.buffer<-crop(aspect.r,dem.buffer)
    shelter.block<-block.sheltermap(dem.block,dir_shelter,interval)
    albedo.block<-crop(albedomap.r,dem.block)
    
    wind.results<-wind.tpsdownscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval,print.res=FALSE,write.files=FALSE )
    wdir.block<-wind.results[[1]]
    wstr.block<-wind.results[[2]]
    invwstr.block<-wind.results[[3]]     
    
    rad.results<-radiation_downscale_stack(day,month,year,hr,sis.r,dnr.r,dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=FALSE)
    direct.block<-rad.results[[1]]*W.to.MJhr
    diffuse.block<-rad.results[[2]]*W.to.MJhr
    total.block<-rad.results[[3]]*W.to.MJhr # = total incoming shortwave radiation - IGNORING ALBEDO EFFECT
    
    elevdif.block<-crop(elevdif.r,dem.block)
    elev.tdif<- (elevdif.block/100)*0.66 
    
    ### CALCULATE LWR (effective resolution 5km) ###
    tref.block<-crop(raster(t100m.day[,,hr+1],template=dem),dem.block)
    rhref.block<-crop(raster(rh.day[,,hr+1],template=dem),dem.block)
    cal.block<-tps.resample(crop(raster(cal.day,layer=hr+1),dem.buffer),dem.block)
    lwr.block<-calc_lwr_block(jd,cal.block,rhref.block,tref.block,plotlwr=FALSE,writefile=FALSE)
    
    #plot(total.block,main=(paste("Total Rad at ",hr,"h00",sep="")))
    #plot(tref.block,main=paste("Tref ",hr))
    #plot(rhref.block,main=paste("RHref ",hr))
    #plot(cal.block,main=paste("CAL ",hr))
    #plot(lwr.block,main=paste("LWR ",hr))
    
    # CALCULATE RAD EFFECT
    RadEffect<- params$estimates[3]*total.block^2.5 +  
      params$estimates[14]*total.block*albedo + 
      params$estimates[4]*lwr.block + 
      params$estimates[5]*albedo+
      params$estimates[16]*lwr.block*wstr.block + 
      params$estimates[15]*total.block^2.5*wstr.block
    ElevEffect<-params$estimates[1]+elev.tdif + 
      params$estimates[6]*wstr.block + 
      params$estimates[7]*invwstr.block
    
    anom.block<-RadEffect+ElevEffect
    t.block<-tref.block+anom.block
    
    if (n==1) tradmap<-total.block else tradmap<-merge(tradmap,total.block)
    if (n==1) lwrmap<-lwr.block else lwrmap<-merge(lwrmap,lwr.block)
    if (n==1) calmap<-cal.block else calmap<-merge(calmap,cal.block)
    
    if (n==1) wstrmap<-wstr.block else wstrmap<-merge(wstrmap,wstr.block)
    if (n==1) wdirmap<-wdir.block else wdirmap<-merge(wdirmap,wdir.block)
    if (n==1) invwmap<-invwstr.block else invwmap<-merge(invwmap,invwstr.block)
    
    if (n==1) RADEFmap<-RadEffect else RADEFmap<-merge(RADEFmap,RadEffect)
    if (n==1) ELEVEFmap<-ElevEffect else ELEVEFmap<-merge(ELEVEFmap,ElevEffect)
    if (n==1) tmap<-t.block else tmap<-merge(tmap,t.block)
    
    
  } # cell
  plot(tradmap,main=paste("Total Rad at ",hr))
  plot(calmap,main=paste("CAL at ",hr))
  plot(lwrmap,main=paste("LWR at ",hr))
  
  plot(wstrmap,main=paste("Wind Str at ",hr))
  plot(wdirmap,main=paste("Wind dir at ",hr))
  #plot(invwmap,main=paste("Inv Wind Str at ",hr))
  
  plot(RADEFmap,main=paste("RAD effect at ",hr))
  plot(ELEVEFmap,main=paste("ELEV effect at ",hr))
  plot(tmap,main=paste("Model T at ",hr))
  
} #hr


################################################


# 4. Test LATENT
# Load day files
filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
dnr.24h<-stack(filein1) ; projection(dnr.24h)<-CRS("+init=epsg:4326")
sis.24h<-stack(filein2) ; projection(sis.24h)<-CRS("+init=epsg:4326")
dnr.24h<-projectRaster(dnr.24h,crs="+init=epsg:27700")
sis.24h<-projectRaster(sis.24h,crs="+init=epsg:27700")
# Load daily files of hourly values of temperature, RH, lwr
t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_100m.r", sep="")
load(t.filein) #t100m.day
# load 5km RH data
rh.filein<-paste(dir_rh5km,"RH_100m_",year,"_",month,"_",day,".R",sep="")
load(rh.filein) # rh.day
# Load and reproject CAL data to OS (downscaling by resampling takes place in hourly interval)
cal.filein<-paste(dir_calday,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
print(paste("CAL file in: ",cal.filein,sep=""))
cal.day<-brick(cal.filein)
cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")

elevdif.r<-raster(paste(dir_terrain,"eref-edem_100m.tif",sep="")) # created by elevation.dif.map function
projection(elevdif.r)<-CRS("+init=epsg:27700")

load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))
interval<-10

albedomap.r<-raster(paste(dir_albedo,"albedo_mean_dembuf.tif",sep=""))
projection(albedomap.r)<-CRS("+init=epsg:27700")

#sheltermap<-block.sheltermap(dem,dir_shelter,interval)

for (hr in 6:10)
{
  tradmap<-raster()
  calmap<-raster()
  lwrmap<-raster()
  wstrmap<-raster()
  wdirmap<-raster()
  invwmap<-raster()
  RADEFmap<-raster()
  ELEVEF<-raster()
  LATEFmap<-raster()
  tmap<-raster()
  
  for (n in 1: length(cells)) {
    print(cells[n]) 
    # Extract hourly data from day stack
    dnr.r<-raster(dnr.24h,layer=hr+1); #plot(dnr.r)
    sis.r<-raster(sis.24h,layer=hr+1); #plot(sis.r)
    
    # calc dem.block/buffer
    x<-landcells[cells[n],1]
    y<-landcells[cells[n],2]
    e.block<-extent(x-2500,x+2500,y-2500,y+2500)
    dem.block<-crop(demuk,e.block) ;# plot(dem.block)
    e.buffer<-extent(x-2500-buffer,x+2500+buffer,y-2500-buffer,y+2500+buffer)
    dem.buffer<-crop(demuk,e.buffer)
    # Slope & aspect cropped to dem.buffer for use in radiation_downscale
    slope.buffer<-crop(slope.r,dem.buffer)
    aspect.buffer<-crop(aspect.r,dem.buffer)
    shelter.block<-block.sheltermap(dem.block,dir_shelter,interval)
    albedo.block<-crop(albedomap.r,dem.block)
    
    wind.results<-wind.tpsdownscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval,print.res=FALSE,write.files=FALSE )
    wdir.block<-wind.results[[1]]
    wstr.block<-wind.results[[2]]
    invwstr.block<-wind.results[[3]]  
    refwstr.block<-wind.results[[4]]  
    
    rad.results<-radiation_downscale_stack(day,month,year,hr,sis.r,dnr.r,dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=FALSE)
    direct.block<-rad.results[[1]]*W.to.MJhr
    diffuse.block<-rad.results[[2]]*W.to.MJhr
    total.block<-rad.results[[3]]*W.to.MJhr # = total incoming shortwave radiation - IGNORING ALBEDO EFFECT
    reftotal.block<-rad.results[[4]]*W.to.MJhr
    
    elevdif.block<-crop(elevdif.r,dem.block)
    elev.tdif<- (elevdif.block/100)*0.66 
    
    ### CALCULATE LWR (effective resolution 5km) ###
    tref.block<-crop(raster(t100m.day[,,hr+1],template=dem),dem.block)
    rhref.block<-crop(raster(rh.day[,,hr+1],template=dem),dem.block)
    cal.block<-tps.resample(crop(raster(cal.day,layer=hr+1),dem.buffer),dem.block)
    lwr.block<-calc_lwr_block(jd,cal.block,rhref.block,tref.block,plotlwr=FALSE,writefile=FALSE)
    
    ### CALCULATE LATENT HEAT FACTORS - SIMPLIFIED###
    if (hr>0) prevtref.block<-crop(raster(t100m.day[,,hr],template=dem),dem.block) else prevtref.block<-tref.block
    prevt.block<-prevtref.block
    
    # 2. Calculate for 5km reference (ignore terrain effects)
    refnetr.block<-(reftotal.block*(1-albedo))-lwr.block # use downscaled lwr to simplify
    # Calculate sea level Pressure for 5km cell at 100m res Prog: prepare.pressurre   - CHECK!!! MORE!
    pref.block<-downscale.pressure(dem.block,p.ncfile,jd,write.file=FALSE) # = Pressure at sea level
    
    # Calculate for 100m incl terrain effects
    p100.block<-correct.pressure(pref.block,tref.block,elevdif.block)   # but what of where no tref 5km cell - find nearest??
    tproxy.block<-tref.block # Simplify tproxy to tref
    netr.block<-(total.block*(1-albedo.block))-lwr.block
    rh.block<-rh.change(tref.block,tproxy.block,rhref.block) 
  
    # 3. Calculate Evapotranspiration
    CRE.5km<-CRE(tref.block,refnetr.block,rhref.block,pref.block,dn,refwstr.block)
    CRE.100m<-CRE(tproxy.block,netr.block,rh.block,p100.block,dn,wstr.block)
    evapdif<-CRE.100m-CRE.5km
    
    # Calculate difference in water condensation - water.conden function
    wc.5km<-Water.conden(prevtref.block,tref.block,rhref.block)
    wc.100m<-Water.conden(prevt.block,tproxy.block,rh.block)
    condendif<-wc.100m-wc.5km 
    
    #plot(total.block,main=(paste("Total Rad at ",hr,"h00",sep="")))
    #plot(tref.block,main=paste("Tref ",hr))
    #plot(rhref.block,main=paste("RHref ",hr))
    #plot(cal.block,main=paste("CAL ",hr))
    #plot(lwr.block,main=paste("LWR ",hr))
    plot(evapdif,main=paste("Evapdif ",hr))
    plot(condendif,main=paste("Condendif ",hr))
    plot(CRE.5km,main="CRE5km")
    plot(CRE.100m,main="CRE100m")
    
    # CALCULATE RAD EFFECT
    RadEffect<- params$estimates[3]*total.block^2.5 +  
      params$estimates[14]*total.block*albedo + 
      params$estimates[4]*lwr.block + 
      params$estimates[5]*albedo+
      params$estimates[16]*lwr.block*wstr.block + 
      params$estimates[15]*total.block^2.5*wstr.block
    ElevEffect<-params$estimates[1]+elev.tdif + 
      params$estimates[6]*wstr.block + 
      params$estimates[7]*invwstr.block
    LatentEffect<-params$estimates[2]*wc.100m + 
      params$estimates[10]*evapdif + 
      params$estimates[11]*condendif
    
    anom.block<-RadEffect+ElevEffect+LatentEffect
    t.block<-tref.block+anom.block
    
    if (n==1) tradmap<-total.block else tradmap<-merge(tradmap,total.block)
    if (n==1) lwrmap<-lwr.block else lwrmap<-merge(lwrmap,lwr.block)
    if (n==1) calmap<-cal.block else calmap<-merge(calmap,cal.block)
    
    if (n==1) wstrmap<-wstr.block else wstrmap<-merge(wstrmap,wstr.block)
    if (n==1) wdirmap<-wdir.block else wdirmap<-merge(wdirmap,wdir.block)
    if (n==1) invwmap<-invwstr.block else invwmap<-merge(invwmap,invwstr.block)
    
    if (n==1) RADEFmap<-RadEffect else RADEFmap<-merge(RADEFmap,RadEffect)
    if (n==1) ELEVEFmap<-ElevEffect else ELEVEFmap<-merge(ELEVEFmap,ElevEffect)
    if (n==1) LATEFmap<-LatentEffect else LATEFmap<-merge(LATEFmap,LatentEffect)
    
    if (n==1) tmap<-t.block else tmap<-merge(tmap,t.block)
    
    
  } # cell
 # plot(tradmap,main=paste("Total Rad at ",hr))
  #plot(calmap,main=paste("CAL at ",hr))
  #plot(lwrmap,main=paste("LWR at ",hr))
  
  #plot(wstrmap,main=paste("Wind Str at ",hr))
  #plot(wdirmap,main=paste("Wind dir at ",hr))
  #plot(invwmap,main=paste("Inv Wind Str at ",hr))
  
  plot(RADEFmap,main=paste("RAD effect at ",hr))
  plot(ELEVEFmap,main=paste("ELEV effect at ",hr))
  plot(LATEFmap,main=paste("LAT effect at ",hr))
  
  plot(tmap,main=paste("Model T at ",hr))
  

} #hr

#####################################
