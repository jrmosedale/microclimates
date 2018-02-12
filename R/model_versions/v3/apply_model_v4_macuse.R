# Apply model to historic data
# v1 - loads most data at daily NOT hourly (v1) step 
# v2 - modification to Latent heat calc - use of 100m res reference data wihtout terrain effects - calc in radiation and wind downsclaing
# Conducts RH and LWR calculations within daily time step
##########################################################################################
# Define INPUTS - dates, 5km cell, parameters
# Most input = daily stack of hourly data
##########################################################################################
args <-commandArgs(trailingOnly = TRUE)
print(args)
start.day <- as.integer(args[1])
start.month<-as.integer(args[2])
start.year<-as.integer(args[3] )
print(paste("Start: ",start.day,"/",start.month,"/",start.year,sep=""))
end.day<-as.integer(args[4] )
end.month<-as.integer(args[5] )
end.year<-as.integer(args[6] )
print(paste("End: ",end.day,"/",end.month,"/",end.year,sep=""))
ukcpcell<-as.integer(args[7])
print(paste("UKCPCELL= ",ukcpcell,sep=""))

##########################################################################################
# Prepare INPUTS for ANY cell
##########################################################################################

print("Calling set up...")
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file
print("Calling functions...")
source("/home/ISAD/jm622/rscripts/apply_model_v2_functions.R") # loads & runs functions
print("Calling inputs...")

### Set MODEL PARAMETERS from Maclean et al 2016 ###
params<-rep(0,12)
intercept<-0
params[1]<-  0.8  # total shortwave rad^2.5
params[2]<- -0.4 # longwave
params[3]<-  3.11   # albedo*total
params[4]<- -0.4 # windspeed*total^2.5
params[5]<- -0.45 # Coast inv windspeed*ldif 
params[6]<-  0.05 # Coast windspeed*ldif 
params[5]<-  0.15 # Coast ldif * sst-tref
params[8]<- -5.8 #evapdif
params[9]<-  0.064 # condendif
params[10]<--0.1 # coldair effect

print("Parameters: ")
print(params)

print(paste("ADIABTIC LAPSE RATE= ",adiabatic.lapserate))
print(paste("BUFFER= ",buffer))
print(paste("UKCPcell= ", ukcpcell) )

##########################################################################################
# Prepare data required for ANY cell
##########################################################################################
### Load WIND data into memory
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))

# Link to sea surface pressure file 
p.ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")
nc_open(p.ncfile)

print("Get Cell Coordinates")
#in.file<-paste(dir_grids,"ukcpmask.grd",sep="")
#gridmask.r<-raster(in.file) #  1=valid cell, NA = sea or not data
#gridmask.r<-crop(gridmask.r,dem) 
gridmask.r<-land5km.r # land5km.r defined in setup
vals<-getValues(gridmask.r)
xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
sel<-which(vals==1)
landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell
print(dim(landcells))

##########################################################################################
# Set time/date/cell parameters
##########################################################################################
cells<-c(718,719,720,744,745,746,769,770,771)

# Frosty day
start.day <- 4; start.month<-4; start.year<-1992
print(paste("Start: ",start.day,"/",start.month,"/",start.year,sep=""))
end.day<-5; end.month<-4; end.year<-1992

# Hot day
start.day <- 24; start.month<-6 ; start.year<- 2012; hr<-0
end.day<-25; end.month<-6; end.year<-2012


print(paste("End: ",end.day,"/",end.month,"/",end.year,sep=""))
start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
jd<-start.jd
print(start.jd);print(end.jd)

print(paste("UKCPCELL= ",ukcpcell,sep=""))

# Set cell and time range for testing
# start.day <- 2; start.month<-7 ; start.year<- 1992; hr<-0
# end.day<-3; end.month<-7; end.year<-1992
# ukcpcell<-926
# lizardcells<-c(916,917,918,919,924,925,926,927,928,929,930,931,932,934,935)
# dir_finalt<-"~/Documents/Exeter/Data2015/Temperature_100m/"

# 

plot.var<-FALSE
plot.effects<-TRUE

##########################################################################################
# APPLY MODEL BY CELL ID to each day and hour 
########################################################################################## 
for (ukcpcell in cells){

#source("/home/ISAD/jm622/rscripts/apply_model_v2_inputs_macuse.R") # defines inputs,parameters, etc
source("~/Documents/Exeter/Rprojects/Project2015/model_versions/v2/apply_model_v2_inputs_macuse.R") # defines inputs,parameters, etc

# Load previous hr data for first time step
#tproxy.block<-prevtref.block 
#prevt.block<-prevtref.block

for (jd in start.jd:end.jd) 
  {
    #ptm <- proc.time()[1:3] 
    year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day ;# values for t
    
    # Load day stack of direct and diffuse RADIATION and reproject to UK OS coordinates
    filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
    filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
    dnr.24h<-stack(filein1) ; projection(dnr.24h)<-CRS("+init=epsg:4326")
    sis.24h<-stack(filein2) ; projection(sis.24h)<-CRS("+init=epsg:4326")
    dnr.24h<-projectRaster(dnr.24h,crs="+init=epsg:27700")
    sis.24h<-projectRaster(sis.24h,crs="+init=epsg:27700")
 
    # Create blank raster stack to test results for one day
    tmodel.day<-stack() # for TESTING ONLY
    radef.day<-stack() # for TESTING ONLY
    latef.day<-stack() # for TESTING ONLY
    coastef.day<-stack() # for TESTING ONLY
    elevef.day<-stack() # for TESTING ONLY
    invwstr.day<-stack() # for TESTING ONLY
    anom.day<-stack() # for TESTING ONLY
    ldif.day<-stack() # for TESTING ONLY
    sstdif.day<-stack() # for TESTING ONLY
    
    # Load prev day temperature file and store last hour of data for block
    if (jd!=start.jd){
    prevt.filein<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_100m.r", sep="")
    load(prevt.filein) #t100m.day
    prevtref.block<-crop(raster(t100m.day[,,23],template=dem),dem.block) # sets prev tref to 23h00 from previous day to jd
    }
    # Load next day's  temperature file and store first hour (for block)cropped later for block??)
    if (jd!=end.jd){
      nextt.filein<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_100m.r", sep="")
    load(nextt.filein) #t100m.day
    nexttref.block<-crop(raster(t100m.day[,,1],template=dem),dem.block)
    }
    # Load daily files of hourly values of temperature, RH, lwr
    t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_100m.r", sep="")
    load(t.filein) #t100m.day
    
    # Use today's data for prev and next if first or last day
    if (jd==start.jd) prevtref.block<-crop(raster(t100m.day[,,1],template=dem),dem.block)
    if (jd==end.jd) nexttref.block<-crop(raster(t100m.day[,,23],template=dem),dem.block)
      
    # load 5km RH data
    rh.filein<-paste(dir_rh5km,"RH_100m_",year,"_",month,"_",day,".R",sep="")
    load(rh.filein) # rh.day
    
    # Load and reproject CAL data to OS (downscaling by resampling takes place in hourly interval)
    cal.filein<-paste(dir_calday,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
    print(paste("CAL file in: ",cal.filein,sep=""))
    cal.day<-brick(cal.filein)
    cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")
    
    # load sst data for day - resample to 100m for block and mask with dem.buffer to set land cells to NA
    infile.sst<- paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")
    sst5km.r<-raster(infile.sst)  
    sea.buffer<-calc(dem.buffer,function(x) ifelse(is.na(x),0,NA))
    sst.buffer<-tps.resample(crop(sst5km.r,sea.buffer),sea.buffer)
    
    ### Calculate ALTITUDE effect on T   ###
    elev.tdif<- (elevdif.block)*adiabatic.lapserate  # lapserate C per m difference
    
    # Set index to track number of hours of data calculated
    hr.result<-1
    print("End of day processes ")  
    #print(proc.time()[1:3] - ptm[1:3])
    
    for (hr in 0:23)  #loop for hours
    {
      # Construct time label for plots 
      timelabel<-paste(DMYjd(jd)$day,"/",DMYjd(jd)$month,"/",DMYjd(jd)$year," ",hr,"h00",sep="")
      print(timelabel)
      
      ### EXTRACT HOURLY BLOCK DATA ###
      # Extract reference temperature data for block and hour
      print("Temperature")
      tref.block<-crop(raster(t100m.day[,,hr+1],template=dem),dem.block)
      # Extract Relative Humidity for hr and for block 
      rhref.block<-crop(raster(rh.day[,,hr+1],template=dem),dem.block)
       # Extract CAL, crop to dem.buffer then resample to dem.block
      cal.block<-tps.resample(crop(raster(cal.day,layer=hr+1),dem.buffer),dem.block)
      cal.block<-calc(cal.block,fun=function(x)ifelse(x<0,0,x))

      ### DOWNSCALE WIND to direction , strength and inverse wind strength for block - Prog: Wind_downscale_blocks    
      # Output: = writes raster for wind str, wind dir and inv wind str (true=print results)      print("Wind")
      wind.results<-wind.tpsdownscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval,print.res=FALSE,write.files=FALSE )
      wdir.block<-wind.results[[1]]
      wstr.block<-wind.results[[2]]
      invwstr.block<-wind.results[[3]]  
      refwstr.block<-wind.results[[4]]  
      
      ### CALC COASTAL EFFECT ###
      print("Coastal")
      # Calculate upwind.sst & sst-tref for block - Prog: upwind_sst_blocks
      # Calculate upwind sea temperature - record NA if no sea cell upwind within 20km
      upwindsst.block<-upwind.sst.block(wdir.block,sst.buffer,dem.buffer,dem.block)
      sst.dif<-sst.tref(upwindsst.block,tref.block)  
      #plot (sst.dif)
      
      # Calculate ldif according to wind direction
      ldif.block<-calc.ldif.block(ldif.stack,wdir.block,interval)    
      
      ### CALC NET RADIATION EFFECT ###
      # Radiation downscale for block -  #output= list(direct.r,diffuse.r,total.r)  
      # Outputs converted to MJ/m2 from W/m2
      print("Start Rad processes ") 
      # Extract hourly data from day stack
      dnr.r<-raster(dnr.24h,layer=hr+1); #plot(dnr.r)
      sis.r<-raster(sis.24h,layer=hr+1); #plot(sis.r)
      rad.results<-radiation_downscale_stack(day,month,year,hr,sis.r,dnr.r,dem.buffer,dem.block,slope.buffer,aspect.buffer,print.results=FALSE)
      direct.block<-rad.results[[1]]*W.to.MJhr
      diffuse.block<-rad.results[[2]]*W.to.MJhr
      total.block<-rad.results[[3]]*W.to.MJhr # = total incoming shortwave radiation - IGNORING ALBEDO EFFECT
      reftotal.block<-rad.results[[4]]*W.to.MJhr
      
      # CALCULATE LWR (effective resolution 5km) #
      print("Start LWR processes ") 
      lwr.block<-calc_lwr_block(jd,cal.block,rhref.block,tref.block,plotlwr=FALSE,writefile=FALSE)
      
      # CALCULATE Rnet reference and refined blocks -  used in latent heat calculations
      Rnetref.block<-(reftotal.block*(1-albedo))-lwr.block # use downscaled lwr to simplify
      Rnet.block<-(total.block*(1-albedo))-lwr.block
      
      ### LATENT HEAT EFFECT CALCS ###
      print("Latent")
      # record whether day or night - used in latent heat calculations
      if (cellStats(total.block,max)>0) dn<-1 else dn<-0
      
      # 1 CHECK THIS Calculate proxy temperature at 100m from tref (5km) and previous anomaly between 5km and 100m
      if (jd==start.jd & hr==0) tproxy.block<-prevtref.block else tproxy.block<-tref.block+prevtanom.block
      # If first reading set previous temp at 100m to previous tref recorded for 5km 
      if (jd==start.jd & hr==0) prevt.block<-prevtref.block
      
      # 2. Calculate sea level Pressure at reference and corrected for terrain Prog: prepare.pressurre   - CHECK!!! MORE!
      pref.block<-downscale.pressure(dem.block,p.ncfile,jd,write.file=FALSE) # = Pressure at sea level
      p100.block<-correct.pressure(pref.block,tref.block,elevdif.block)   # but what of where no tref 5km cell - find nearest??
      
      # 3. Calculate RH  for 100m incl terrain effects
      rh.block<-rh.change(tref.block,tproxy.block,rhref.block) 
      
      # 4. Calculate Evapotranspiration
      CRE.5km<-CRE(tref.block,Rnetref.block,rhref.block,pref.block,dn,refwstr.block)
      CRE.100m<-CRE(tproxy.block,Rnet.block,rh.block,p100.block,dn,wstr.block)
      evapdif<-CRE.100m-CRE.5km
      
      # 5. Calculate difference in water condensation - water.conden function
      wc.5km<-Water.conden(prevtref.block,tref.block,rhref.block)
      wc.100m<-Water.conden(prevt.block,tproxy.block,rh.block)
      condendif<-wc.100m-wc.5km 
      
      ### FLOW ACCUMULATION EFFECT - NIGHT ONLY ###
      if (dn==0){
        # Test if thermal inversion conditions (tic) ie low wind and high lwr emission
        tic.lwr<-0.2 # 400Wm2 1.4MJm2hr ?
        tic.wstr<-0.5
        tic<- overlay(lwr.block,wstr.block,fun=function(x,y){ifelse(x>tic.lwr & y<tic.wstr,1,0)})  
      }
      if (dn==1) {tic<-raster(matrix(rep(0,100),nrow=NROW(dem.block),ncol=NCOL(dem.block)),template=dem.block)}
      coldair.block<-tic*(altdif.block*adiabatic.lapserate)*twi.block
      
      ####################################################################################
      # Apply model paramters to calculate temperature anomaly
      ####################################################################################
      #t.block<-tref.block+elev.effect+(rad.effect*2)+(latent.effect*0.3)+(coast.effect/2)
      print("Calculating 100m temperatures...")
      RadEffect<- params[1]*calc(total.block,function(x){ifelse(x>0,x^2.5,0)}) + 
                  params[2]*lwr.block + 
                  params[3]*total.block*albedo + 
                  params[4]*calc(total.block,function(x){ifelse(x>0,x^2.5,0)})*wstr.block
      CoastEffect<-params[5]*invwstr.block * ldif.block +
                   params[6]*wstr.block* ldif.block +
                   params[7]*ldif.block*sst.dif 
      LatentEffect<-params[8]*evapdif + 
                    params[9]*condendif
      ElevEffect<-elev.tdif # 
      FlowEffect<-params[10]*coldair.block

      anom.block<-intercept + RadEffect+CoastEffect+LatentEffect+ElevEffect+FlowEffect
     
      #### Calculate new temperature raster at 100m resolution and anomaly from T5km ####
      t.block<-tref.block+anom.block # CREATES new t.block from tref and modelled anomaly

      ####################################################################################
      #### Plot results ####
      ####################################################################################
      par(mfrow=c(3,3))
      if (plot.var==TRUE){
          plot(t.block,main=paste("T100m ",ukcpcell," at ",timelabel,sep=""))
          plot(t.block-prevt.block,main="T hourly change")
          plot(anom.block,main="Tanom")
          plot(prevtref.block,main="Prev Tref")
          plot(tref.block,main="Tref")
          plot(rhref.block,main=paste("RH ref ",timelabel,sep=""))
          
          plot(prevt.block,main="PrevT")
          plot(tproxy.block,main="Tproxy")
          plot(rh.block,main=paste("RH 100m ",timelabel,sep=""))
          
          #plot(total.block,main=paste("Total Rad ",timelabel,sep=""))
          #plot(Rnet.block,main=paste("Net Rad",timelabel,sep=""))
          ##plot(direct.block,main=paste("Direct ",sep=""))
          #plot(diffuse.block,main=paste("Diffuse ",sep=""))
          #plot(rh.5km,main=paste("RH 5km ",timelabel,sep=""))
          #plot(rh.block,main=paste("RH 100m ",timelabel,sep=""))
          #plot(CRE.5km,main=paste("Evap Ref ",timelabel,sep=""))
          #plot(evapdif,main=paste("Evap difference  ",sep=""))    
          #plot(wc.5km,main=paste("Cond ref ",timelabel,sep=""))
          #plot(wc.100m,main=paste("Cond corrected ",timelabel,sep=""))
          #plot(condendif,main=paste("Cond difference ",sep=""))
          #plot(crop(aspect.buffer,dem.block),main=paste("Aspect ",timelabel,sep=""))
          #plot(wstr.block,main="Wstr")
          #plot(wdir.block,main="Wdir")
          #plot(p100.block,main=paste("P 100m ",timelabel,sep=""))
          #plot(lwr.block,main=paste("LWR ",sep=""))
          #plot(cal.block,main="CAL")
          #plot(flow.block*tic,main="flow*tic")
          #plot(ldif.block,main="Ldif")
          #plot(sst.dif,main="SST dif")
          #plot(invwstr.block*ldif.block,main="invwstr x ldif")
          #plot(wstr.block*sst.dif,main="wstr x sst.dif")
          #plot(ldif.block*sst.dif,main="ldif x sst .dif")
          #plot(dem.block)
      }
      if (plot.effects==TRUE){
        plot(t.block,main=paste("T100m ",ukcpcell," at ",timelabel,sep=""))
        plot(t.block-prevt.block,main="T hourly change")
        plot(anom.block,main="Tanom")
        
        #plot(tref.block,main="Tref block")
        #plot(anom.block,main="Tanom") 
        
        plot(RadEffect,main="Radiation Effect ")
        plot(Rnet.block,main="Rnet")
        plot(cal.block,main="CAL")
        plot(LatentEffect,main="Latent Heat Effect")
        
        plot(CoastEffect,main=paste("Coast Effect wdir= ",cellStats(wdir.block,stat="mean",na.rm=TRUE) ,sep=""))
        
        #plot(ElevEffect,main="Elev Effect ")
        plot(FlowEffect,main="Flow Effect")
      }  
    
      # OUTPUTS: SAVE DOWNSCALED TEMP FOR DEM.BLOCK for every Hour 
      tmodel.day<-stack(tmodel.day,t.block)
      invwstr.day<-stack(invwstr.day,invwstr.block)
      radef.day<-stack(radef.day,RadEffect)
      coastef.day<-stack(coastef.day,CoastEffect); 
      elevef.day<-stack(elevef.day,ElevEffect); 
      latef.day<-stack(latef.day,LatentEffect);
      anom.day<-stack(anom.day,anom.block)
      ldif.day<-stack(ldif.day,ldif.block)
      sstdif.day<-stack(sstdif.day,sst.dif)
      
      
      #if (ukcpcell%%10==0){plot(dem.block, main=paste("UKCP09 cell= ",ukcpcell,sep=""))}
      # Save tref and tanomaly.100m for next time step
      prevtref.block<-tref.block
      prevt.block<-t.block
      prevtanom.block<-anom.block
      
      # Print summary statistics of downscaled temperature
      print(paste("Met Off 5km temperature:",cellStats(tref.block,mean),sep="") )
      print(paste("Mean 100m temperature:",cellStats(t.block,mean),sep="") )
      print(paste("Max 100m temperature:",cellStats(t.block,max),sep="") )
      print(paste("Min 100m temperature:",cellStats(t.block,min),sep="")  )                     
      
      # Save hourly result as 3D array
      print(paste("Number of hours of data calculated= ",hr.result,sep=""))
      #tmodel.r[,,hr.result]<-getValues(t.block,format="matrix") # convert to 3d array
      hr.result<-hr.result+1 # advance index ready for next hour
      
    }# end loop for every hour
    
    # FOR TESTING only
    fileout<-paste(dir_finalt,"tmodel-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(tmodel.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"radef-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(radef.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"coastef-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(coastef.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"elevef-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(elevef.day,fileout,format="GTiff",overwrite=TRUE)
    
	  fileout<-paste(dir_finalt,"latef-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(latef.day,fileout,format="GTiff",overwrite=TRUE)
    
	  fileout<-paste(dir_finalt,"anom-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(anom.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"invwstr-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(wstr.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"ldif-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(ldif.day,fileout,format="GTiff",overwrite=TRUE)

    fileout<-paste(dir_finalt,"sstdif-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(sstdif.day,fileout,format="GTiff",overwrite=TRUE)
    
} # end loop for every day

    
# WRITE OUTPUT as 3D array of yearly file of hourly values
#fileout<-paste(dir_finalt,"block-",ukcpcell,"-",year,".R",sep="")

#print(paste("Writing results file ",fileout,sep=""))
#save(tmodel.r,file=fileout)
#print(proc.time()[1:3] - ptm[1:3])

} # END FOR UKCPCELL

##########################################################################################
# Define results array - ASSUME 50x50 row block
numcells<-length(cells)
#max.hr<-(end.jd-start.jd+1)*24
#tmod<-array(NA,c(50,50,max.hr))
#radef<-array(NA,c(50,50,max.hr))
#latef<-array(NA,c(50,50,max.hr))
#coastef<-array(NA,c(50,50,max.hr))

par(mfrow=c(2,3))
map.r<-raster()

# 2. Print hourly Temperature maps for 24hr
for (hr in 0:23){
  blocks0<-vector("list",numcells)
  blocks1<-vector("list",numcells)
  blocks2<-vector("list",numcells)
  blocks3<-vector("list",numcells)
  blocks4<-vector("list",numcells)
  blocks5<-vector("list",numcells)
  blocks6<-vector("list",numcells)
  
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
    
    infile<-paste(dir_finalt,"tmodel-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")  
    print(infile)
    blocks0[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"anom-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks1[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"radef-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks2[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"coastef-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks3[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"ldif-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks4[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"sstdif-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks5[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"latef-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks6[[i]]<-raster(stack(infile),layer=hr+1)
    
    i<-i+1
  }
  
  ### PLOT merged blocks
  map.r<-do.call(merge, blocks0)
  titletext<-paste("Tmodel at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks1)
  titletext<-paste("Anom at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks2)
  titletext<-paste("Rad at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks3)
  titletext<-paste("Coast at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks4)
  titletext<-paste("ldif at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks5)
  titletext<-paste("sst at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
} # end 24 hr
