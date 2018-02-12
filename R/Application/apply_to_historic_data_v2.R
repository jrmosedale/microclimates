# Progs:
#wind_downscale_blocks
#radprog_blocks
#upwind_sst_blocks
#calc_ldif_block
#prepare_pressure
#latentheat.block

# Set time range for which data will be analysed
start.year<-2000
start.month<-5
start.day<-10
hr<-0
end.year<-2000
end.month<-5
end.day<-11
start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
jd<-start.jd
### Load PARAMETERS ###
parameters<-"~/Documents/Exeter/Data2015/parameters/testparams.csv"
params<-read.csv(parameters)


##########################################################################################
# Define BLOCKS for analysis in parallel = 5km blocks matching UKCP09 grid cells 
# Progs:  wind_downscale_blocks
#         radprog_blocks
##########################################################################################
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
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))

# Link to sea surface pressure file - already unpacked
p.ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")

# block.width<-5000 = ASSUMPTION
buffer<-20000 # set for max required for any individual program (coast effect???)


###### LOOP for every block... #######
ukcpcell<-933 # for testing only lizRD - 936-939

#for (ukcpcell in(1:length(landcells))){
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
  albedo<-0.2
  
  # Extract FLOW ACC for area
  infile<-paste(dir_flowacc,"flowacc_multi.tif",sep="")
  flowacc<-raster(infile,res=100)
  minval<-cellStats(flowacc,min)  
  flowacc<-flowacc/minval  # divide by min val (area of single cell?)
  flow.block<-crop(flowacc,dem.block)

  ####### LOOP FOR EACH TIME PERIOD from START to END DATE - how to loop using jd in hour steps, starting at 12.00??  ########
  for (jd in start.jd:end.jd ) 
  {
    year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day ;# values for t
    
    # Load day stack of direct and diffuse RADIATION
    filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
    filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),sep="")
    dnr.24h<-stack(filein1)
    sis.24h<-stack(filein2)
    
    # Create blank raster stack to hold results for one day
    final.t.stack<-stack()
    
    # load 5km RH data
    rh.filein<-paste(dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R",sep="")
    load(rh.filein) # rh.day
    
    # Load prev day temperature file and store last hour of temperature data from previous day before loading new files
    prevt.filein<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),".r", sep="")
    load(prevt.filein) #t5km.day
    t5km.prevhr<-raster(t5km.day[,,23],template=grid5km.r) # sets prev hr to 23h00 from previous day to jd
    prevtref.block<-ref5km.to.block100m(dem.block,t5km.prevhr) # extracts at 100m for block
       
    # Load daily files of hourly values of temperture, RH, lwr
    t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),".r", sep="")
    load(t.filein) #t5km.day
    
    # load lwr file
    lwr.filein<-paste(dir_lwr,"lwr_",year,"_",month,"_",day,".R",sep="")
    load(lwr.filein) # lwr.day
    
    # load sst data for day - resample to 100m for block and mask with dem.buffer to set land cells to NA
    infile.sst<- paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")
    sst5km.r<-raster(infile.sst)  
    sst.buffer<-resample(sst5km.r,dem.buffer)
    sst.buffer<-mask(sst.buffer,dem.buffer, inverse=TRUE) 
    
    for (hr in 0:23)  #loop for hours
    {
      # Construct time label for plots 
      timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
      print(timelabel)
      
      ### REF TEMPERATURES ###
      # Extract 5km reference temperature data  for block at 100m resolution - uses: ref5km.to.block100m
      print("Temperature")
      tref5km.r<-raster(t5km.day[,,hr+1],template=grid5km.r)
      tref.block<-ref5km.to.block100m(dem.block,tref5km.r)
      
      ### CALC RADIATION EFFECT ###
      # Downscale wind to direction , strength and inverse wind strength for block - Prog: Wind_downscale_blocks    
      # Output: = writes raster for wind str, wind dir and inv wind str (true=print results)
      print("Wind")
      wind.results<-wind.tpsdownscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval=10,print.res=FALSE,write.files=TRUE )
      wdir.block<-wind.results[[1]]
      wstr.block<-wind.results[[2]]
      invwstr.block<-wind.results[[3]]     
      #invwstr.dif<-invwstr.ref-invwstr.block
      
      # Radiation downscale for block - Prog: radprog_blocks - #output  = list(direct.r,diffuse.r,total.r)  
      # Outputs converted to MJ/m2 from W/m2
      print("Radiation")
      # Extract hourly data from day stack
      dnr.r<-raster(dnr.24h,layer=hr+1); #plot(dnr.r)
      sis.r<-raster(sis.24h,layer=hr+1); #plot(sis.r)
      rad.results<-radiation_downscale_stack(day,month,year,hr,sis.r,dnr.r,dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=FALSE)
      direct.block<-rad.results[[1]]*W.to.MJhr
      diffuse.block<-rad.results[[2]]*W.to.MJhr
      total.block<-rad.results[[3]]*W.to.MJhr # = total incoming shortwave radiation - IGNORING ALBEDO EFFECT
      
      # Extract lwr for hour and block at 100m - effective res of 5km
      lwr5km.r<-raster(lwr.day[,,hr+1],template=grid5km.r)
      lwr.block<-ref5km.to.block100m(dem.block,lwr5km.r)
      
      # NOT REQUIRED - CAL for comparison purposes
      #cal.infile<-paste(dir_calimp,"CALimp_5km_",year,"_",month,"_",day,".R",sep="")
      #load(cal.infile) # =CALimp.day
      #CAL.r<-raster(CALimp.day[,,hr+1],template=grid5km.r)
      #CAL.5km<-crop(CAL.r,dem.block)  

      ### CALC COASTAL EFFECT ###
      print("Coastal")
      # Calculate upwind.sst & sst-tref for block - Prog: upwind_sst_blocks
      # Calculate upwind sea temperature - record NA if no sea cell upwind within 20km
      upwindsst.block<-upwind.sst.block(wdir.block,sst.buffer,dem.buffer,dem.block)
      sst.dif<-sst.tref(upwindsst.block,tref.block)  
      #plot (sst.dif)
      
      # Calculate ldif according to wind direction
      ldif.block<-calc.ldif.block(ldif.stack,wdir.block,interval)    
      
      ### LATENT HEAT EFFECT CALCS ###
      print("Latent")
      # record whether day or night - used in latent heat calculations
      if (cellStats(total.block,max)>0) dn<-1 else dn<-0
      
      # Extract Relative Humidity for hr and for block 
      rh5km.r<-raster(rh.day[,,hr+1],template=grid5km.r)
      rhref.block<-ref5km.to.block100m(dem.block,rh5km.r)
      
      # Calculate sea level Pressure for 5km cell at 100m res Prog: prepare.pressurre   - CHECK!!! MORE!
      pref.block<-downscale.pressure(dem.block,p.ncfile,jd,write.file=FALSE) # = Pressure at sea level
      # p.sw<-downscale.pressure(ncfile,grid5km.r,t,write.file=TRUE)
      # p.blk<-crop(p.sw,dem.block)
      # Correct for elevation variation within 5km cell 
      p100.block<-correct.pressure(pref.block,tref.block,elevdif.block)   # but what of where no tref 5km cell - find nearest??
      
      # Calculate latent heat difference Prog latentheat.block.functions  
      
      # CHECK THIS Calculate proxy temperature at 100m from tref (5km) and previous anomaly between 5km and 100m
      if (jd==start.jd & hr==0) tproxy.block<-prevtref.block else tproxy.block<-tref.block+prevtanom.block
      # If first reading set previous temp at 100m to previous tref recorded for 5km 
      if (jd==start.jd & hr==0) prevt.block<-prevtref.block
      
      # Calculate relative humdidity at 100m res
      rh.block<-rh.change(tref.block,tproxy.block,rhref.block) 
      
      # Calcute difference in Evapotranspiration - CRE function CRE<-function(Temp,Rn,RH,P,dn,u2)
      # 1. Calculate net radiation for block from incoming shortwave, albedo and outgoing longwave radiation
      netr.block<-(total.block*(1-albedo))-lwr.block
      # 2. Calculate mean Rn and Wind speed for 5km cell
      wstr.mean<-matrix(cellStats(wstr.block, stat='mean', na.rm=TRUE),nrow=nrow(wstr.block),ncol=ncol(wstr.block) )
      wstr.mean.block<-raster(wstr.mean,template=wstr.block)
      wstr.mean.block<-mask(wstr.mean.block,wstr.block)
      netr.mean<-matrix(cellStats(netr.block, stat='mean', na.rm=TRUE),nrow=nrow(netr.block),ncol=ncol(netr.block) )
      netr.mean.block<-raster(netr.mean,template=netr.block)
      netr.mean.block<-mask(netr.mean.block,netr.block)   
      # 3. Calculate Evapotranspiration
      CRE.5km<-CRE(tref.block,netr.mean.block,rhref.block,pref.block,dn,wstr.mean.block)
      CRE.100m<-CRE(tproxy.block,netr.block,rh.block,p100.block,dn,wstr.block)
      evapdif<-CRE.100m-CRE.5km
      
      # Calculate difference in water condensation - water.conden function
      wc.5km<-Water.conden(prevtref.block,tref.block,rhref.block)
      wc.100m<-Water.conden(prevt.block,tproxy.block,rh.block)
      condendif<-wc.100m-wc.5km  
      
      ### ALTITUDE EFFECT ON TEMPERATURE uses ref-dem ###
      elev.tdif<- (elevdif.block/100)*0.66  # altitude effect=1.98 per 300m - see Wikipedia!
      
      ### FLOW ACCUMULATION EFFECT ###
      # Test if thermal invesrsion conditions (tic) ie low wind and high lwr emission
      tic.lwr<-1.45 # 400Wm2 1.4MJm2hr ?
      tic.wstr<-1
      tic<- overlay(lwr.block,wstr.block,fun=function(x,y){ifelse(x>tic.lwr & y<tic.wstr,1,0)})  
      
      ####################################################################################
      # Apply model paramters to calculate temperature anomaly
      ####################################################################################
      #t.block<-tref.block+elev.effect+(rad.effect*2)+(latent.effect*0.3)+(coast.effect/2)
      print("Calculating 100m temperatures...")
      anom.block<-(params$estimates[1]+elev.tdif+   # difference in temperature due to altitude
               #params$estimates[2]*wc.100m+   # effect of water condensation due to temporal change in temperature at site
               params$estimates[3]*total.block^2+      # net short wave radiation
               params$estimates[4]*lwr.block+  # net longwave radiation
               params$estimates[5]*albedo+        # albedo
               params$estimates[6]*wstr.block+    # wind speed one metre above ground
               params$estimates[7]*invwstr.block+    #  = 1/(sqrt(Site.u1)+1) 
               params$estimates[8]*ldif.block+     # lsrdif = Culdrose land sea ratio - site land sea ratio
               params$estimates[9]*sst.dif+    # tempdif = SST - Culdrose temperature
               #params$estimates[10]*evapdif+   # evaporation at site - evaporation at culdrose
               #params$estimates[11]*condendif+ # condensation at site - condensation at Culdrose
               params$estimates[12]*flow.block+     # flow accumulation
               params$estimates[13]*tic+      # temperature inversion conditions (binary)
               params$estimates[14]*total.block*albedo+
               params$estimates[15]*total.block*wstr.block+
               params$estimates[16]*lwr.block*wstr.block+
               params$estimates[17]*invwstr.block*ldif.block+
               params$estimates[18]*wstr.block*sst.dif+
               params$estimates[19]*ldif.block*sst.dif+
               params$estimates[20]*flow.block*tic
               )
      
      
      
      
      
      #### Calculate new temperature raster at 100m resolution and anomaly from T5km ####
      t.block<-tref.block+anom.block
      tanom.block<-t.block-tref.block # Check correct order/sign - used in next hour timestep
      
      ####################################################################################
      #### Plot results ####
      par(mfrow=c(3,3))
      plot(t.block,main=paste("T100m ",ukcpcell," at ",timelabel,sep=""))
      plot(t.block-prevt.block,main="T hourly change")
      plot(tanom.block,main="Tanom")
      #plot(tref.block,main="Tref")
      #plot(crop(aspect.buffer,dem.block),main=paste("Aspect",timelabel,sep=""))
      plot(wstr.block,main="Wstr")
      #plot(wdir.block,main="Wdir")
      #plot(p100.block,main=paste("P 100m ",timelabel,sep=""))
      #plot(rhref.block,main=paste("RH ref ",timelabel,sep=""))
      #plot(rh.5km,main=paste("RH 5km ",timelabel,sep=""))
      #plot(rh.block,main=paste("RH 100m ",timelabel,sep=""))
      #plot(CRE.5km,main=paste("Evap Ref ",timelabel,sep=""))
      #plot(CRE.100m,main=paste("Evap 100m  ",sep=""))      
      #plot(wc.5km,main=paste("Cond Ref ",timelabel,sep=""))
      #plot(wc.100m,main=paste("Cond 100m ",sep=""))
      plot(total.block,main=paste("Total Rad ",timelabel,sep=""))
      #plot(direct.block,main=paste("Direct ",sep=""))
      #plot(diffuse.block,main=paste("Diffuse ",sep=""))
      #plot(lwr.block,main=paste("LWR ",sep=""))
      #plot(CAL.5km,main="CAL")
      #plot(flow.block*tic,main="flow*tic")
      plot(sst.dif,main="SST dif")
      
      # OUTPUTS: SAVE DOWNSCALED TEMP FOR DEM.BLOCK for every Hour 
      final.t.stack<-stack(final.t.stack,t.block)
      #if (ukcpcell%%10==0){plot(dem.block, main=paste("UKCP09 cell= ",ukcpcell,sep=""))}
      
      # Delete downscaled block files except final hourly temperature files
      
      # Save tref and tanomaly.100m for next time step
      prevtref.block<-tref.block
      prevt.block<-t.block
      prevtanom.block<-tanom.block
      
    }# end loop for every hour
    
    # WRITE OUTPUT to Daily file of Hourly values
    names(final.t.stack)<-0:23
    writeRaster(final.t.stack,file=paste(dir_finalt,"block_",ukcpcell,"_on_",day,"_",month,"_",year,".tif",sep=""), format="GTiff",overwrite=TRUE)
    print(paste(dir_finalt,"block_",ukcpcell,"_on_",day,"_",month,"_",year,".tif",sep=""))
    # load by: test<-stack(paste(dir_finalt,"name.tif",sep=""))
    
  } # end loop for every day
  
  
#} # end for every block

##########################################################################################
# Re-load and stitch together block results if required using Mosaic etc
for (t in start.jd:end.jd by ??? ) 
{
  year<-DMYjd(t)$year; month<-DMYjd(t)$month ; day<-DMYjd(t)$day
  result.lyr<-raster(demsw)
    
  for (ukcpcell in(1:length(landcells)))
  { # create day stack of hourly data for whole of area
    day.block<-paste(paste(dir_newt,"block_",ukcpcell,"_on_",day,"_",month,"_",year,".tif",sep=""))
    block.stack<-stack(in.block)
    sw.stack<-moasic(sw.stack,block.stack)
  }
  writeRaster(newt.stack,file=paste(dir_newtsw,"SW_on_",day,"_",month,"_",year,".tif",sep=""),overwrite=TRUE)
}




# plot hourly time series 
#  day<-
par(mfrow=c(3,4))

for (hr in 0:23){
  timelabel<-paste(day,"/",month,"/",year," ",hr,"h00",sep="")
  tfinal<-stack(paste(dir_finalt,"block_",ukcpcell,"_on_",day,"_",month,"_",year,".tif",sep=""))
  plot(tfinal,hr+1,main=timelabel) 
  #temp.plot(raster(tfinal,layer=hr+1),hr,day,month,year)
}



##########################################################################################
# FUNCTIONS USED IN ABOVE
##########################################################################################
# LATENTHEAT.BLOCK
# This function calculates crop reference evapotranspiration (using the Penman-Monteith equation).
# Inputs and outputs are as rasters excpet for dn which may be raster or a single value
# Details of algorithm here: http://www.fao.org/docrep/x0490e/x0490e00.htm
# Input variables:
# Temp is the temperature at the site in degrees C. You will need to use the anomoly between the
#5 km grid and the 100m cell in the previous time-step to estimate this (as the CRE function is
# used to derive the local temperature)
# Net radiation  - importantly this is in MJ m-2 hour-1 and may require conversion of units from the
# satellite derived estimates (typical value ~0.2)
# RH - relative humidity expressed as a percentage (typical value 80%)
# P - Atmposheric rressure - in millibars (typical value ~ 1000)
# currently not assumed to vary by location within a 5km grid cell, but we could do an altitude correction:
# http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric pressure (p)
# dn  # a binary variable specifying whether it is day (1)  or   night (0) needed as the equation for Soil Heat Flux
# changes depending on whether it is night or day - can be a Raster or single value
# u2 is wind speed at at 2 m height [m s-1] (typical value ~5)
# Output variable:
#Crop reference evapotranspiration [mm m-2 hr-1]
CRE<-function(Temp,Rn,RH,P,dn,u2){
  e0<-0.6108*exp(17.27*Temp/(Temp+237.3)) # saturated vapour pressure
  ea<-e0*(RH/100) # actual vapour pressure
  delta<-4098*(0.6108*exp(17.27*Temp/(Temp+237.3)))/((Temp+237.3)^2)  # slope vapour pressure curve
  if(class(dn)=="RasterLayer") G<-overlay(Rn,dn,fun = function(x, y) ifelse(y>0, 0.1*Rn, 0.5*Rn)) else {
    if (dn>0) G<-0.1*Rn else G<-0.5*Rn } # Soil heat flux
  gamma<-0.000665*(P/10)  # psychrometric constant
  ET0<-(0.408*delta*(Rn-G)+gamma*37/(Temp+273)*u2*(e0-ea))/(delta+gamma*(1+0.34*u2))
  ET0<-calc(ET0,fun=function(x){ifelse(x>0,x,0)})
  ET0
}

# This function calculates the change in relative humidity between two locations as a result of the temperature change
# Input variables:
# t1 is the reference temperature in degrees C - i.e. the value for each 5 km grid cell
# t2 is the is the temperature at the site in degrees C. You will need to use the anomoly between the
#5 km grid and the 100m cell in the previous time-step to estimate this
# rh is the reference relative humdity, expressed as a percentage - i.e. the value for each 5 km grid cell
# Output variable:
# the relative humdity at the site, expressed as a percentage - i.e. the value for each 100m grid cell
rh.change<-function(t1,t2,rh){
  # absolute humidity at t1
  e0.1<-0.6108*exp(17.27*t1/(t1+237.3))
  e<-e0.1*(rh/100)
  a.hum.1<-(2165*e)/(t1+273.16) # grams per metre cubed
  # rel humidity at t2
  s.e<-(a.hum.1*(t2+273.16))/2165
  s.e0<-0.6108*exp(17.27*t2/(t2+237.3))
  rhs<-(s.e/s.e0)*100
  rhs
}

# This function calculates the amount of water that can be expected to condense,
# as either a result of a change in temperature from one place to anotehr or through time
# Input variables:all as RASTERS
# t1 is the reference temperature in degrees C - i.e. the value for each 5 km grid cell
# t2 is the is the temperature at the site in degrees C. You will need to use the anomoly between the
# 5km grid and the 100m cell in the previous time-step to estimate this
# rh is the reference relative humdity, expressed as a percentage - i.e. the value for each 5km grid cell
# Output variable:
# the amount of water condensed [mm m-2 hr-1]. Zero if relative humidity is less than 100% - as RASTER
Water.conden<-function(t1,t2,rh){
  # absolute humidity at t2, rh=100
  e0.100<-0.6108*exp(17.27*t2/(t2+237.3))
  e100<-e0.100*(100/100)
  a.hum.100<-(2165*e100)/(t2+273.16) # grams per metre cubed
  # absolute humidity at rh2
  rh2<-rh.change(t1,t2,rh)
  e0.2<-0.6108*exp(17.27*t2/(t2+237.3))
  e2<-e0.2*(rh2/100)
  a.hum.2<-(2165*e2)/(t2+273.16) # grams per metre cubed
  a.hum.2<-a.hum.2-a.hum.100
  #a.hum.2<-ifelse(rh2>100,a.hum.2,0)
  a.hum.2<-overlay(a.hum.2,rh2,fun=function(x,y) ifelse(y>100,x,0) )
  a.hum.2
}
##########################################################################################
# PREPARE_PRESSURE
# Prepare Sea level Pressure data - extract data at 100m resolution for time jd
# Requires JD functions
# Works for whole area if block=dembuf or similar or for small block
# Write raster for area grid

downscale.pressure<-function(block,p.ncfile,jd,write.file=FALSE)
{
  # unzip .gz file if necessary
  #gzfile<-paste(dir_pressure025,"pp_0.25deg_reg_v11.0.nc.gz",sep="")
  #ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")
  #gunzip(filename=gzfile, destname=ncfile, overwrite=TRUE)
  #ncdf_test<-nc_open(ncfile)
  
  # Select time bands 
  # Base jd value = 1/1/1950
  jd.base<-JDdmy(1,1,1950)
  jd.band<-jd-jd.base
  
  # Load raster from ncdf file 
  pressure.r<-raster(p.ncfile,band=jd.band)
  #plot(pressure.r)
  
  # Reproject cropped version to OSGB projection 
  pressure.r<-projectRaster(crop(pressure.r,c(-7,0,49,52)),crs="+init=epsg:27700")
  #plot(pressure.r)
  
  # Resample to 5km
  p5km.r<-resample(pressure.r,grid5km.r)
  p5km.r<-mask(p5km.r,grid5km.r)
  
  # The convert to 100m grid using nearest neighbour method 
  p.block<-ref5km.to.block100m(dem.block,p5km.r)
  
  # Resample & crop to 100m resolution of dembuf
  #p100m.r<-resample(pressure.r,block)
  #p100m.r<-mask(p100m.r,block)
  
  # Crop and write pressure raster (all times) for area of interest
  if(write.file){
    plot(p.block,main=paste("Sea level Pressure 100m cell for ",DMYjd(jd)$day,"/",DMYjd(jd)$month, "/",DMYjd(jd)$year,sep=""))
    outfile<-paste(dir_pressure,"pressure.tif",sep="")  
    writeRaster(file=outfile,p.block,overwrite=TRUE )
  }
  return(p.block)
} # end function

# Correct pressure for elevation according to hypsometric formula
# Inputs: sea level pressure, tmp at location, elevation above sea level
# References:  http://keisan.casio.com/exec/system/1224579725, http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric%20pressure
correct.pressure<-function(sl.pressure,tmp,elevation)
{
  kelvins<-273+tmp
  pressure<-sl.pressure*((kelvins-(0.0065*elevation)) / kelvins )^5.25588
  return(pressure)
}
##########################################################################################
# Calculates single raster  ldif.block
# Input: ldif.stack of ldif for each wind direction
#        wind direction for block
#        direction interval for which ldif calculated

calc.ldif.block<-function(ldif.stack,wdir.block,interval=10){
  ldif.block<-raster
  wdir.layer<-round(wdir.block/interval) # creates raster holding layer in ldif.stack for wind direction
  wdir.layer<-calc(wdir.layer,fun=function(x){ifelse(x==0,36,x)}) # corrects layer 0 to 36 (angle 0 to 360)
  ldif.block<-stackSelect(ldif.stack,wdir.layer)
  return(ldif.block)
} # end function
##########################################################################################

# Functions to calculate Upwind SST and SST-Tref for each 100m cell in block given wind direction and sst map
# Input:  100m dem grid -  block and buffered regions
#         raster of sst for buffer attime  t
#         raster of wind direction for block at time t
# Output: raster of SST and SST-Tref at t

#####################################################################
# FUNCTIONS 
# Used by upwind.sst
findsst<-function(x)
{
  is.sea<-ifelse(x==0,NA,1)
  nearest.sea<-match(1,is.sea)
  if (is.na(nearest.sea)){
    seacell<-NA 
  } else {seacell<-x[nearest.sea]}  
  return(seacell)
}
#####################################################################
# UPWIND>SST BLOCK 
# Finds nearest upwind sea cell within 'buffer' (20km) for a block of cells for a SINGLE wind direction 
# If no sea cell found returns value of...0 sea cells returned = NA
# Otherwise records upwind sst held in sst.r
# Input:  sstbuffer.r - sea surface temperatures for buffer region
#         gridbuffer.r (including buffer) - could be dem
#         gridblock.r - defines area block of interest within buffer - could be dem
#         direction assumed to be constant across area
#         distance = max distance a search for nearest sea cell (= buffer of 20km)
# called by: upwind.sst.block
upwind.sst<-function(sstbuffer.r,gridbuffer.r,gridblock.r,direction,distance=10000)
{
  x<-dim(sstbuffer.r)[1]
  y<-dim(sstbuffer.r)[2]
  step<-distance/res(sstbuffer.r)[1] # max number of cells from focal cell to be searched
  
  # create matrix holding cell numbers for sea cells and 0 for landcells
  vals<-getValues(gridbuffer.r,format="matrix") # land/sea 1/NA values
  cells<-getValues(sstbuffer.r,format="matrix")
  seacells<-ifelse(is.na(vals),cells,0) # matrix of -999 if land or sea temperature if sea - could be changed to sst values as long as land =0
  #plot(raster(seacells,template=gridbuffer.r))
  
  store<-array(0,dim=c(dim(seacells)[1]-(2*step),dim(seacells)[2]-(2*step),step+1)) # 3d array to hold values for non-buffered region for every 'step'
  store[,,1]<-seacells[(step+1):(dim(seacells)[1]-step),(step+1):(dim(seacells)[2]-step)]
  #plot(raster(store[,,1],template=gridblock.r))
  
  for (i in 1:step)
  {
    xshift<-round(i*sin(direction*(pi/180)),0) 
    yshift<-round(i*cos(direction*(pi/180)),0)
    #print(paste("i: ",i, " xshift: ",xshift," yshift: ",yshift,sep=""))
    yshift<-yshift*(-1)
    store[,,(i+1)]<-seacells[(step+1+yshift):(dim(seacells)[1]-step+yshift),(step+1+xshift):(dim(seacells)[2]-step+xshift)]
    
  } # end
  
  # use first/last to find nearest sea cell??
  storev<-array(store,dim=c((dim(store)[1]*dim(store)[2]),step+1))
  
  upwind.sst<-apply(storev,1,findsst) # if no sea within buffer then = NA?
  upwind.sst<-matrix(upwind.sst,nrow=nrow(store),ncol=ncol(store)) 
  
  upwind.sst.r<-raster(upwind.sst,template=gridblock.r)
  upwind.sst.r<-raster::mask(upwind.sst.r,gridblock.r)
  #plot(upwind.sst.r,main=paste("nearest seacell where direction= ",direction,sep=""))
  
  return(upwind.sst.r)
}# end function

#####

# Function to return upwind.sst for block given w.dir for block
upwind.sst.block<-function(wdir.block,sst.buffer,dem.buffer,dem.block, plotresult=FALSE)
{
  # define output raster
  sst.vals<-array(NA,dim=c(nrow(wdir.block),ncol(wdir.block)))
  # Round wind direction to nearest degree
  wdir.block<-round(wdir.block,0)
  # Calc number of unique wind dir - create raster stack of upwind.sst for each unique w.direction
  wdir.vals<-unique(round(getValues(wdir.block)[which(!is.na(getValues(wdir.block)))])) # every unique wind.val excluding NA and rounding to nearest degree
  
  # Create stack of upwind sst for block for every unique wind direction udirng time t
  upwind.vals<-array(0,dim=c(nrow(wdir.block),ncol(wdir.block),length(wdir.vals)))  
  for (i in 1:length(wdir.vals)){
    upwind.r<-upwind.sst(sst.buffer,dem.buffer,dem.block,wdir.vals[i],distance=20000)
    upwind.vals[,,i]<-getValues(upwind.r, format="matrix")
  }
  # assign sst according to w.dir of cell
  wdir.class<-cbind(wdir.vals,1:length(wdir.vals)) # create 2col vector for reclassifying
  wdir.class.r<-reclassify(wdir.block,wdir.class)  #reclass so that values = layer of stack holding sst values 
  class.vals<-getValues(wdir.class.r, format="matrix")
  
  for (i in 1:length(wdir.vals)){
    sel<-which(class.vals==i, arr.ind=TRUE)
    upwind.i<-upwind.vals[,,i]
    sst.vals[sel]<-upwind.i[sel]
    #plot(raster(sst.vals,template=wdir.block))
  }
  result<-raster(sst.vals,template=wdir.block)
  if (plotresult) plot(result)
  return(result)
} # end function

#######
# Function to calculate sst-tref 
sst.tref<-function(sst.r,tref.r)
{
  if (compareRaster(sst.r,tref.r)!=TRUE){warning("!!sst-tref rasters not comparable!!")}
  sst.tref.r<-overlay(sst.r,tref.r,fun=function(x,y){ifelse(is.na(x)&!is.na(y),0,y-x)})
  return(sst.tref.r)
} # end function

#####################################################################
# WIND_DOWNSCALE FUNCTIONS
#####################################################################
# REQUIRES: Function JDdmy for computing Julian data to work out number of days after 1st Jan 1960 (1st wind data observation),
# so that correct element of array can be extracted

# Based on hour, day, month and year, extracts the required value from the array of values stored by wind_downscale1.R
# Inputs:
# hr: the hour (0-23) 0 = midnight
# day: the day of the month (0-31)
# month: the month of the year (numeric: 0-12)
# year: any year from 1960 to 2014
# Output: the element of the array stored by wind_downscale1.R that corresponds to either that hour, or the latest
# period immediatly before that hour (data only available 6-hourly)
array.val<-function(hr,day,month,yr)
{
  jd.base=JDdmy(1,1,1960)
  jd<-JDdmy(day,month,yr)
  dval<-(jd-jd.base)*4
  hval<-floor(hr/6)
  val<-dval+hval+1
  val
}

####################################################
# A. Prepare and load data for all time periods - Requires dem.block and e.block
####################################################
# 1. Load shelter maps for block using interval to which wind direction is rounded - eg 10 = every 10 degrees
# using interval, e.block, dem.block
block.sheltermap<-function(dem.block,dir_shelter,interval=10){
  dem.m<-getValues(dem.block,format="matrix")
  shelter<-array(NA, dim=c((360/interval),nrow(dem.m),ncol(dem.m)))
  for (i in 1:(360/interval)) {
    dir<-i*interval
    in.file<-paste(dir_shelter,"Shelter_",sprintf("%03d",dir,sep=""),"_deg.tif",sep="")
    print(in.file)
    wcoef.r<-raster(in.file) 
    wcoef.r<-crop(wcoef.r,dem.block)
    shelter[i,,]<-getValues(wcoef.r,format="matrix") # fills shelter[i,1:end,1] to shelter[i,1:10,end]
    print(paste("i= ",i," dir= ",dir))
  }   
  return(shelter)
} # end function


# 2. Load wind data (single file for whole time period - created by wind_downscale1)
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v
#load(file=paste(dir_wind,"wind_u.r",sep=""))
#load(file=paste(dir_wind,"wind_v.r",sep=""))


####################################################
# # # This Function downscales the wind to 100m cells for a specific hour
####################################################
# Specify hour,day, month and year for which data are required
#year=2010;month=6;day<-10;hr<-11
# Stages:
# (1) get wind values for a given hour, day, month and year
# (2) convert to 100m resolution raster OSGB grid reference using THIN SPLINET INTERPOLATION - see sepaate function
# (3) adjust based on altitude
# (4) adjust based on shelter coefficient

wind.tpsdownscale<-function(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter,interval=10,print.results=TRUE,write.files=FALSE)
{
  # Can be quite slow. Allows you to keep tabs on progress by printing hour, day, month & year
  tp<-paste("year=",year," month=",month," day=",day," hour=",hr,sep="")
  print(tp)
  dir_windstrength<-paste(dir_wind,"strength/",sep="")
  dir_winddirection<- paste(dir_wind,"direction/",sep="")
  dir_windinvstr<- paste(dir_wind,"invstr/",sep="")
  
  #############
  # Stage 1: get wind values for a given day month and year
  #############
  # As original data are 4x daily, but data are required for each hour,
  # this bit reads in the data for the periods immediatly before after for which there are data and calculates
  # weighted mean
  av1<-array.val(hr,day,month,year)
  av2<-av1+1
  rem<-hr/6-floor(hr/6)
  uwind1<-wind_u[,,av1]
  uwind2<-wind_u[,,av2]
  vwind1<-wind_v[,,av1]
  vwind2<-wind_v[,,av2]
  uwind<-(1-rem)*uwind1+rem*uwind2
  vwind<-(1-rem)*vwind1+rem*vwind2
  
  #############
  # Stage 2: convert to 100m resolution raster OSGB grid reference using tps resampling
  #############
  # Convert to raster (original lat long format and resolution - CHECK LAT/LON MAX/MIN
  uwind.r<-raster(uwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
  vwind.r<-raster(vwind,xmn=-7.5,xmx=0,ymn=47.5,ymx=52.5)
  # Reproject in OSGB projection
  crs(uwind.r)<-latlong
  crs(vwind.r)<-latlong
  u_osgb<-projectRaster(uwind.r,crs=ukgrid)
  v_osgb<-projectRaster(vwind.r,crs=ukgrid)
  # tps resampling to 100m using FUNCTION tps.resample
  u_100<-tps.resample(u_osgb,dem.block)
  v_100<-tps.resample(v_osgb,dem.block)
  
  #############
  # Stage 3: adjust based on altitude of terrain
  # NB Height adjustment based on wind spped values at different pressures downloaded  Earth System Research Lab
  # Typical heights at different pressures calculated from Allen et al 1998 http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric pressure (p)
  # Quadratic function fitted - NB this works well for heights up to ~1800m. IT won't work above ~2000m
  # Function was first derived by comparing values at different pressures (heights) over the course of a year (2014)
  #############
  # adjust wind speeds by height of dem
  # convert to matrices
  uwind.m<-getValues(u_100,format="matrix")
  vwind.m<-getValues(v_100,format="matrix")
  dem.m<-getValues(dem.block,format="matrix")
  # adjust wind by height
  ustr<-sqrt(uwind.m^2) # wind strength
  vstr<-sqrt(vwind.m^2) # wind strength
  udir<-ifelse(uwind.m>0,1,-1) # positive or negative
  vdir<-ifelse(vwind.m>0,1,-1) # positive or negative
  u.adj<-ustr*((-0.000000108025)*dem.m^2+0.000408692*dem.m+0.956139) # NB don't worry about warnings. Calculation assigns NAs to the sea
  v.adj<-vstr*((-0.000000108025)*dem.m^2+0.000408692*dem.m+0.956139)  # NB don't worry about warnings. Calculation assigns NAs to the sea
  # adjust values to correspond to wind speed 1m above theground
  # rescaling factor first derived by comparing values ot Culdrose wind data using wind_downscale3
  # however, in line wiht what you'd expect from: http://www.fao.org/docrep/x0490e/x0490e07.htm#wind profile relationship
  u.adj<-u.adj*0.373686439
  v.adj<-v.adj*0.373686439
  u.m<-u.adj*udir
  v.m<-v.adj*vdir
  u.r<-raster(u.m,template=dem.block) #only for checking
  v.r<-raster(v.m,template=dem.block)# only for checking
  # plot(u.r,main="u.r");plot(v.r,main="v.r")
  
  #############
  # Stage 4: height adjustments done using shelter coefficient maps based on topography and wind direction
  ############# 
  # Calculate Wind Direction
  dir.m <- (180/pi)*(atan2(u.m,v.m))  # NB this is direction in which wind blows to
  dir.m<-ifelse(dir.m<=180,dir.m+180,dir.m-180) # NB this direction from which wind originates (360 deg)
  dir.r<-raster(dir.m,template=dem.block)
  
  # Calculate Wind Strength from u.m and v.m components
  str.m<-matrix(NA,nrow=(nrow(dem.m)),ncol=(ncol(dem.m)) )# matrix for storing all values
  str.m<-sqrt(u.m^2+v.m^2)
  wstr.r<-raster(str.m,template=dem.block)# only for checking
  
  #  plot(wstr.r,main="wstr") ; plot(dir.r,main="dir")
  
  # Uses Rounded wind direction of each cell to select correct shelter coefficient
  dir.m<-round(dir.m/interval)*interval # round to interval used for shelter maps
  dir.m<-ifelse(dir.m==0,360,dir.m) # converts 0 to 360 degree direction 
  # Applies shelter coefficient to calculate wind strength if land cell (else NA)
  mxrws<-nrow(dem.m)
  mxcls<-ncol(dem.m)
  for (rws in 1:mxrws) {
    for (cls in 1:mxcls) {    
      str.m[rws,cls]<-str.m[rws,cls]*shelter[(dir.m[rws,cls]/interval),rws,cls] 
    }
  }  
  
  #############
  # Stage 5: Format outputs
  #############
  # Calculate inverse wind strength
  invstr.m<-1/(sqrt(str.m+1))
  
  # Convert to raster and save tif files 
  wstr.r<-raster(str.m,template=dem.block)
  wdir.r<-raster(dir.m,template=dem.block)
  invwstr.r<-raster(invstr.m,template=dem.block)
  if (write.files==TRUE){
    fileout.1<-paste(dir_windstrength,"strength_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    fileout.2<-paste(dir_winddirection,"direction_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    fileout.3<-paste(dir_windinvstr,"invstr_",year,"_",month,"_",day,"_",hr,".tif",sep="")
    #print(fileout.1); print(fileout.2); print(fileout.3)
    writeRaster(wstr.r,file=fileout.1,overwrite=TRUE)
    writeRaster(wdir.r,file=fileout.2,overwrite=TRUE)
    writeRaster(invwstr.r,file=fileout.3,overwrite=TRUE)
  }
  #Create raster stack and print if requested
  if (print.results==TRUE){
    shelterblock<-raster(shelter[(mean(dir.m,na.rm=TRUE)/interval),,],template=dem.block)
    result.stack<-stack(u.r,v.r,shelterblock,wstr.r,wdir.r,invwstr.r)
    names(result.stack)<-c("u.r","v.r","shelter block","wind strength","wind direction","inv wind str")
    par(mfrow=c(2,3))
    plot(result.stack)  
  } # end if
  
  wind.results<-c(wdir.r,wstr.r,invwstr.r)
  return(wind.results)
  
} # end function wind.downscale



####################################################
# RADPROG FUNCTIONS
####################################################
# Converts a raster to a matrix for use with solar index functions NOT USED - USE getValues(r,format="matrix")
use.raster<-function(r)
{
  xr<-dim(r)[1]
  xc<-dim(r)[2]
  m<-array(getValues(r),dim=c(xc,xr))
  m<-t(m) # transpose
  m
}

# OS GB grid ref to Lat and Long
OSGBtolatlong<-function(x,y)
{
  pt = data.frame(x,y)
  coordinates(pt)=~x+y
  proj4string(pt)=CRS("+init=epsg:27700")
  latlong<-spTransform(pt,CRS("+init=epsg:4326"))
  ll<-as.data.frame(latlong)
  ll
}

# Needed for solar index function
solartime <- function(localtime,Long,Julian,merid=0,dst=0)
{
  Bn <- 2 * 3.141 * (Julian - 81) / 364
  eot <- 9.87 * sin(2 * Bn) - 7.53 * cos(Bn) - 1.5 * sin(Bn)
  solartime <- localtime + (4 / 60) * (2 * 3.141 * (merid - Long) / 360) + (eot / 60) - dst
  solartime
}

solalt <- function(localtime,Lat,Long,Julian,merid=0,dst=0)
{
  stime<-solartime(localtime,Long,Julian,merid,dst)
  tt <- 0.261799 * (stime - 12)
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((Julian - 171) / 365.25))
  Sinh = sin(declin) * sin(Lat * pi / 180) + cos(declin) * cos(Lat * 3.141 / 180) * cos(tt)
  solalt = (180 * atan(Sinh / sqrt(1 - Sinh * Sinh))) / pi
  solalt
}

solazi <- function(localtime,Lat,Long,Julian,merid=0,dst=0)
{
  stime<-solartime(localtime,Long,Julian,merid,dst)
  tt = 0.261799 * (stime - 12)
  declin = (pi * 23.5 / 180) * cos(2 * pi * ((Julian - 171) / 365.25))
  Sinh = sin(declin) * sin(Lat * pi / 180) + cos(declin) * cos(Lat * pi / 180) * cos(tt)
  hh = (atan(Sinh / sqrt(1 - Sinh * Sinh)))
  Sinazi = cos(declin) * sin(tt) / cos(hh)
  cosazi = (sin(Lat * pi / 180) * cos(declin) * cos(tt) - cos(pi * Lat / 180) * sin(declin)) / sqrt((cos(declin) *
                                                                                                       sin(tt)) ^ 2 + (sin(pi * Lat / 180) * cos(declin) * cos(tt) - cos(pi * Lat / 180) * sin(declin)) ^ 2)
  solazi = 180 + (180 * atan(Sinazi / sqrt(1 - Sinazi * Sinazi))) / pi
  if (cosazi < 0) {
    if (Sinazi < 0) {
      solazi = 180 - solazi
    } else {
      solazi = 540 - solazi
    }
  }
  solazi
}

horizonangle <- function(dtm,azimuth,res=100)
{
  dtm<-(dtm*5)/res
  azimuth<-azimuth-90
  azi <- azimuth * (pi/180)
  horizon <- array(0,dim(dtm))
  dtm3 <- array(0,dim(dtm)+200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)] <- dtm
  for (step in 1:10) {
    horizon[1:x,1:y] <- pmax(horizon[1:x,1:y], (dtm3[(101+sin(azi)*step^2):(x+100+sin(azi)*step^2),(101+cos(azi)*step^2):(y+100+cos(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(5*step^2))
  }
  horizon
}

solarindex <- function(slope,aspect,localtime,Lat,Long,Julian,dtm=array(0,dim=c(1,1)),res=100,merid=0,dst=0,shadow=TRUE)
{
  saltitude<-array(solalt(localtime,Lat,Long,Julian,merid,dst),dim(dtm))
  alt <- saltitude * (pi/180)
  zen <- pi/2 - alt
  sazimuth<-array(solazi(localtime,Lat,Long,Julian,merid,dst),dim(dtm))
  azi <- sazimuth * (pi/180)
  sl <- slope * (pi/180)
  asp <- aspect * (pi/180)
  shadowmask <- array(1,dim(dtm))
  horangle<-horizonangle(dtm,sazimuth)
  #plot(raster(horangle,template=dem.buffer))
  if(shadow) {
    shadowmask[horizonangle(dtm,sazimuth)>tan(alt)] <- 0
  }
  index <- array(0,dim(dtm))
  index <- cos(zen) * cos(sl) + sin(zen) * sin(sl) * cos(azi - asp)
  index[index<0] <- 0
  index <- index * shadowmask
  #plot(raster(index,template=dem.buffer),main="Solar Index")
  index
}

skyview <-function(dtm,steps=36)
{
  sky <- array(1,dim(dtm))
  for (s in 1:steps) {
    sky <- sky-atan(horizonangle(dtm,s*360/steps))/((pi/2)*steps)
  }
  sky
}

####################################################################################
# Direct normal radiation: to downscaled direct radiation
####################################################################################

radiation_downscale_stack<-function(day,month,year,hr,sis.r,dnr.r, dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=TRUE)
{ 
  # Calculate lat and long of centre of grid
  ll<-OSGBtolatlong(xmin(dem.buffer)+(0.5*(xmax(dem.buffer)-xmin(dem.buffer))) , ymin(dem.buffer)+(0.5*(ymax(dem.buffer)-ymin(dem.buffer))) )
  lat<-as.numeric(ll[2])
  long<-as.numeric(ll[1])
  # Calculate if day or night from sunvector
  jd.h<-JDdmy(day,month,year)+((-12+hr)/24 )# Julian day with hour fraction
  if (sunvector(jd.h,lat,long,0)[3]>=0) # if DAYTIME 
  {
    # resample to 100m and extent of dem.buffer
    #dnr.buffer<-raster::resample(dnr.r,dem.buffer) # change this to different interpolation method?
    #sis.buffer<-raster::resample(sis.r,dem.buffer)
    # crop to 3* dem.buffer then tps model and interpolate
    xdim<-(xmax(dem.buffer)-xmin(dem.buffer))
    ydim<-(ymax(dem.buffer)-ymin(dem.buffer))
    e.tps<-extent(xmin(dem.buffer)-xdim,xmax(dem.buffer)+xdim,ymin(dem.buffer)-ydim,ymax(dem.buffer)+ydim)# Run setup programs for creating constant raster maps etc
    
    sis.r<-crop(sis.r,e.tps)
    sis.buffer<-tps.resample(sis.r,dem.buffer,FALSE)
    dnr.r<-crop(dnr.r,e.tps)
    dnr.buffer<-tps.resample(dnr.r,dem.buffer,FALSE)
    
    #plot(sis.buffer); plot(dnr.buffer)
    
    # work out Julian day and time - ncdf time in hrs since 1/1/1983 0:00
    #jul.base<-JDdoy(1,1983)
   # hrs<-ncvar_get(ncdf_dnr,"time")
    #days<-floor(hrs/24)
   # jul.day<-as.numeric(jul.base+days) #JD at 12:00 on ncdf day
   # h<-as.numeric(hrs%%24)
   # if (h!=hr) {print("WARNING h^= hr!!!")}     
    
    # creates raster template for storing direct and diffuse radiation values with cell values of -9
    direct.r<-dem.buffer*0-9
    diffuse.r<-dem.buffer*0-9
    total.r<-dem.buffer*0-9
    
    # convert values to matrices for use with solar index function
    m.dem<-getValues(dem.buffer,format="matrix") # use.raster function above
    m.slope<-getValues(slope.buffer,format="matrix")
    m.aspect<-getValues(aspect.buffer,format="matrix")
    
    # converts NAs to zeros - NB: buffer of boundary blocks = NA->0
    sel<-which(is.na(m.dem)==T); m.dem[sel]<-0
    sel<-which(is.na(m.slope)==T); m.slope[sel]<-0
    sel<-which(is.na(m.aspect)==T); m.aspect[sel]<-0
    
    # Calculate solar index for dem 
    si<-solarindex(slope=m.slope,aspect=m.aspect,localtime=hr,
                   Lat=lat,Long=long,Julian=jd,dtm=m.dem)
   # Calculate single value solar index for flat slope
    slope.flat<-0;aspect.flat<-0
    si.flat<-solarindex(slope=slope.flat,aspect=aspect.flat,localtime=hr,
                        Lat=lat,Long=long,Julian=jd,shadow=F)[1,1]
    #plot(raster(si,template=dem.buffer),main="Solar index")
   
    # Direct normal radiation:
    dnr.m<-getValues(dnr.buffer,format="matrix")
    # Direct radiation: flat
    dir.m<-dnr.m*si.flat # i.e. SID=Cos(sza)*DNR
    # Diffuse radiation flat - 
    sis.m<-getValues(sis.buffer,format="matrix")
    dif.m<-sis.m-dir.m
    #plot(raster(dir.m,template=dem.buffer),main="Direct - flat")
    #plot(raster(dif.m,template=dem.buffer),main="Diffuse - flat")
   
    # Calculate matrix coords that define central 'block'
    xmn<-1+(xmin(dem.block)-xmin(dem.buffer))/res(dem.block)[1]
    xmx<-dim(dem.buffer)[1]-(xmax(dem.buffer)-xmax(dem.block))/res(dem.block)[1]
    ymn<-1+(ymin(dem.block)-ymin(dem.buffer))/res(dem.block)[1]
    ymx<-dim(dem.buffer)[2]-(ymax(dem.buffer)-ymax(dem.block))/res(dem.block)[1]
    
    # downscaled direct radiation
    direct<-dnr.m*si;   #  direct<-dir.m*si ?????
    direct<-direct[xmn:xmx,ymn:ymx] # Extract block
    #plot(raster(direct,template=dem.buffer),main="Old Direct")
    # downscaled diffuse radiation
    sv<-skyview(m.dem)
    diffuse<-dif.m*sv
    diffuse<-diffuse[xmn:xmx,ymn:ymx] # Extract block
    # add mask back in so that SIs are only produced for land
    mask<-getValues(dem.block,format="matrix")*0
    direct<-direct+mask
    diffuse<-diffuse+mask
    total<-direct+diffuse
    
    # Extract central Block and Convert to rasters
    direct.r<-raster(direct,template=dem.block)
    diffuse.r<-raster(diffuse,template=dem.block)
    total.r<-raster(total,template=dem.block) 
    
  } else {         # IF NIGHTIME         
    direct.r<-dem.block*0
    diffuse.r<-dem.block*0
    total.r<-dem.block*0
  }
  
  #Create raster stack and print if requested
  if (print.results){
    result.stack<-stack(direct.r,diffuse.r,total.r)
    names(result.stack)<-c("direct","diffuse","total")
    par(mfrow=c(2,2))
    plot(result.stack)  
  } # end if
  
  results<-c(direct.r,diffuse.r,total.r)
  return(results)
} # end function


###################################################################################
#FUNCTION FOR RESAMPLING RASTER using thin plate spline
#library(fields) # required for thin plate spline
tps.resample<-function(input.r,output.r,msk=TRUE){
  xy <- data.frame(xyFromCell(input.r,1:ncell(input.r)))
  v <- getValues(input.r)
  tps <- Tps(xy, v) # fit tps model (Don't worry about warning)
  result<- raster(output.r) # create blank raster at resolution of output.r
  
  # use model to predict values for all locations
  result<- interpolate(result,tps)
  if (msk==TRUE) result<-mask(result,output.r)
  #plot(result,main="Thin-plate spline - min")
  
  return(result)
}# end function

  

