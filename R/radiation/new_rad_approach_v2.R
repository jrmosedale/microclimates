# Set time range for which data will be analysed
start.year<-1992
start.month<-7
start.day<-30
hr<-0
end.year<-1992
end.month<-7
end.day<-15
start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
jd<-start.jd

##########################################################################################
#  BLOCK LOOP 
##########################################################################################
#ukcpcell<-930 # for testing only lizRD - 936-939
for(ukcpcell in c(931,932,933)){


# Get ukcp09 grid cells coordinates
in.file<-paste(dir_grids,"ukcpmask.grd",sep="")
print(in.file)
gridmask.r<-raster(in.file) #  1=valid cell, NA = sea or not data
gridmask.r<-crop(gridmask.r,e.dem) 
vals<-values(gridmask.r)
xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
sel<-which(vals==1)
landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell

### load WIND data into memory
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))

# Link to sea surface pressure file - already unpacked
p.ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")

# block.width<-5000 = ASSUMPTION
buffer<-20000 # set for max required for any individual program (coast effect???)

# Extract data for BLOCK & BUFFER REGION
x<-landcells[ukcpcell,1]
y<-landcells[ukcpcell,2]
e.block<-extent(x-2500,x+2500,y-2500,y+2500)
dem.block<-crop(demuk,e.block)
e.buffer<-extent(x-2500-buffer,x+2500+buffer,y-2500-buffer,y+2500+buffer)
dem.buffer<-crop(demuk,e.buffer)
plot(dem.buffer,main=ukcpcell);plot(dem.block,main=ukcpcell)

# Perform operations for each block but fixed for all timeperiods - ie CROPPING Terrain datasets
# Wind - load shelter matrix for block - Prog: wind_downscale_blocks
shelter.block<-block.sheltermap(dem.block,dir_shelter)

# Slope & aspect cropped to dem.buffer for use in radiation_downscale
slope.buffer<-crop(raster(paste(dir_terrain,"slope.tif",sep="")),dem.buffer)
aspect.buffer<-crop(raster(paste(dir_terrain,"aspect.tif",sep="")),dem.buffer)

# Calculate elevation difference from 5km mean
elevdif.filein<-paste(dir_terrain,"eref-edem_100m.tif",sep="") # created by elevation.dif.map function
elevdif.block<-crop(raster(elevdif.filein),dem.block)
#elevdif.block<-dem.block-cellStats(dem.block, stat='mean', na.rm=TRUE)

# Load and extract ldif for each wind direction for block - used with hourly wind direction to calc coastal effect
for (direction in seq(0,350,interval)) {
  ldif.filein<-paste(dir_ldif,"ldif_",direction,"deg_from_percland_in_",radius/1000,"km.tif",sep="")
  ldif.layer<-crop(raster(ldif.filein),dem.block)
  if (direction==0) ldif.stack<-stack(ldif.layer) else ldif.stack<-stack(ldif.stack,ldif.layer)
} # end for direction

# Extract ALBEDO for area
albedo.block<-raster(array(0.2,c(dim(dem.block)[1],dim(dem.block)[2])),template=dem.block)

# Extract FLOW ACC for area
infile<-paste(dir_flowacc,"flowacc_multi.tif",sep="")
flowacc<-raster(infile,res=100)
minval<-cellStats(flowacc,min)  
flowacc<-flowacc/minval  # divide by min val (area of single cell?)
flow.block<-crop(flowacc,dem.block)

##########################################################################################
# DAY LOOP 
##########################################################################################
jd<-start.jd
year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day ;# values for t

# Create blank raster stack to hold results for one day
final.t.stack<-stack()

# Load daily RADIATION files
filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
dnr.24h<-stack(filein1)
sis.24h<-stack(filein2)

# DOWNSCALE RAD
rad.24h<-rad_downscale_100m(dem.buffer,dnr.24h,sis24h)

# Load daily TEMPERATURE FILES
#prev day temperature file and store last hour of temperature data from previous day before loading new files
prevt.filein<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),".r", sep="")
load(prevt.filein) #t5km.day
t5km.prevhr<-raster(t5km.day[,,23],template=grid5km.r) # sets prev hr to 23h00 from previous day to jd
prevtref.block<-ref5km.to.block100m(dem.block,t5km.prevhr) # extracts at 100m for block
# Load daily files of hourly values of temperture, RH, lwr
t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),".r", sep="")
load(t.filein) #t5km.day

# load 5km RH data
rh.filein<-paste(dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R",sep="")
load(rh.filein) # rh.day

# load lwr file
lwr.filein<-paste(dir_lwr,"lwr_",year,"_",month,"_",day,".R",sep="")
load(lwr.filein) # lwr.day

# load sst data for day - resample to 100m for block and mask with dem.buffer to set land cells to NA
infile.sst<- paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")
sst5km.r<-raster(infile.sst)  
sst.buffer<-resample(sst5km.r,dem.buffer)
sst.buffer<-mask(sst.buffer,dem.buffer, inverse=TRUE) 

### Load PARAMETERS ###
parameters<-"~/Documents/Exeter/Data2015/parameters/testparams.csv"
params<-read.csv(parameters)
par(mfrow=c(2,2))

##########################################################################################
# HOUR LOOP 
##########################################################################################
for (hr in 7:13){
  # Construct time label for plots 
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  print(timelabel)
  
  ### 5km TEMPERATURE Calculations ###
  # Extract 5km reference temperature data  for block at 100m resolution - uses: ref5km.to.block100m
  print("Temperature")
  tref5km.r<-raster(t5km.day[,,hr+1],template=grid5km.r)
  tref.block<-ref5km.to.block100m(dem.block,tref5km.r)
  
  # WIND    
  print("Wind")
  wind.results<-wind.tpsdownscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval=10,print.res=FALSE,write.files=TRUE )
  wdir.block<-wind.results[[1]]
  wstr.block<-wind.results[[2]]
  invwstr.block<-wind.results[[3]]   
  #plot(wdir.block,main=paste("Wind dir ",timelabel,sep=""))
  
  # RADIATION CALCULATIONS
  print("Radiation")
  dnr.buffer<-raster(rad.24h[[1]],hr+1)
  sis.buffer<-raster(rad.24h[[2]],hr+1)
  plot(crop(sis.buffer,dem.block),main=paste("SIS ",timelabel,sep=""))
  # Calculate Rcorr use: calculate_Rcorr
  R.block<-calculate_Rcorr(jd,hr,sis.buffer,dnr.buffer,dem.buffer,slope.buffer,aspect.buffer,dem.block)
  #plot(R.block,main="Rc")
  
  # LWR Calcualtions
  lwr5km.r<-raster(lwr.day[,,hr+1],template=grid5km.r)
  lwr.block<-ref5km.to.block100m(dem.block,lwr5km.r)
  
  # COASTAL EFFECT
  print("Coastal")
  upwindsst.block<-upwind.sst.block(wdir.block,sst.buffer,dem.buffer,dem.block)
  sst.dif<-sst.tref(upwindsst.block,tref.block)  
  #plot (sst.dif,main=paste("SST dif ",timelabel,sep=""))
  # Calculate ldif according to wind direction
  ldif.block<-calc.ldif.block(ldif.stack,wdir.block,interval)   
  #plot(ldif.block,main=paste("L dif ",timelabel,sep=""))
  
  ### FLOW ACCUMULATION EFFECT ###
  # Test if thermal invesrsion conditions (tic) ie low wind and high lwr emission
  tic.lwr<-1.45 # 400Wm2 1.4MJm2hr ?
  tic.wstr<-1
  tic<- overlay(lwr.block,wstr.block,fun=function(x,y){ifelse(x>tic.lwr & y<tic.wstr,1,0)})  
  
  
  ##########################################################################################
  # CALCULATE  TEMPERATURE MAP 
  ##########################################################################################
  
  # APPLY Rcorr to tref
  # Possible apply Rcorr to increase/decrease from prev timestep rather than actual temp values???
  trad.block<-tref.block*Rcorr.block
  plot(trad.block,main=paste("Trad at ",timelabel,sep=""))
  
  elev.tdif<-elev.tdif<- (elevdif.block/100)*0.66  # altitude effect=1.98 per 300m - see Wikipedia!
  
  # APPPLY other factors
  anom.block<-(params$estimates[1]+elev.tdif+   # difference in temperature due to altitude
                 params$estimates[4]*lwr.block+  # net longwave radiation
                 params$estimates[5]*albedo.block+        # albedo
                 params$estimates[6]*wstr.block+    # wind speed one metre above ground
                 params$estimates[7]*invwstr.block+    #  = 1/(sqrt(Site.u1)+1) 
                 params$estimates[8]*ldif.block+     # lsrdif = Culdrose land sea ratio - site land sea ratio
                 params$estimates[9]*sst.dif+    # tempdif = SST - Culdrose temperature
                 params$estimates[12]*flow.block+     # flow accumulation
                 params$estimates[13]*tic+      # temperature inversion conditions (binary)
                 #params$estimates[14]*total.block*albedo+
                 #params$estimates[15]*total.block*wstr.block+
                 params$estimates[16]*lwr.block*wstr.block+
                params$estimates[17]*invwstr.block*ldif.block+
                params$estimates[18]*wstr.block*sst.dif+
                params$estimates[19]*ldif.block*sst.dif+
                params$estimates[20]*flow.block*tic
  )
  plot(anom.block,main="Other Anom")
  t.block<-trad.block+anom.block
  plot(t.block,main=paste("Temperature at ",timelabel,sep=""))
  #plot(dem.block,main="dem")
  
  # OUTPUTS: SAVE DOWNSCALED TEMP FOR DEM.BLOCK for every Hour 
  final.t.stack<-stack(final.t.stack,t.block)
  } # end hr loop

  # WRITE OUTPUT to Daily file of Hourly values
  print(paste(dir_finalt,"block_",ukcpcell,"_on_",day,"_",month,"_",year,".tif",sep=""))
  writeRaster(final.t.stack,file=paste(dir_finalt,"block_",ukcpcell,"_on_",day,"_",month,"_",year,".tif",sep=""), format="GTiff",overwrite=TRUE)

} # loop block














