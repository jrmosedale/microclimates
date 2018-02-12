# Apply model to historic data
# v1 - loads most data at daily NOT hourly (v1) step 
# v2 - modification to Latent heat calc - use of 100m res reference data wihtout terrain effects - calc in radiation and wind downsclaing
# Conducts RH and LWR calculations within daily time step
# Designed to calculate for whole year for ONE 5km block
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
cellnum<-as.integer(args[7])
print(paste("Cells= ",cells,sep=""))

### TEST INPUTS ###
#start.day <- 2; start.month<-7; start.year<-1992
#end.day<-3; end.month<-7; end.year<-1992

##########################################################################################
# Prepare INPUTS for ANY cell
##########################################################################################
print("Calling set up...")
source("/home/ISAD/jm622/rscripts/setup_v20km_carson.R") # loads & runs setup file
print("Calling functions...")
source("/home/ISAD/jm622/rscripts/apply_model_v20km_functions.R") # loads & runs functions

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

print(paste("ADIABTIC LAPSE RATE= ",adiabatic.lapserate))
print(paste("BUFFER= ",buffer))
print(paste("INTERVAL= ",interval))

##########################################################################################
# Set time/date/cell parameters
##########################################################################################

print(paste("Start: ",start.day,"/",start.month,"/",start.year,sep=""))
print(paste("End: ",end.day,"/",end.month,"/",end.year,sep=""))
start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
print(start.jd);print(end.jd)

##########################################################################################
# For every  block - usually single block for main runs ie cells = single value
##########################################################################################
  print("Calling inputs...")
  source("/home/ISAD/jm622/rscripts/apply_model_v20km_inputs.R") # defines inputs,parameters, etc
  #source("~/Documents/Exeter/RProjects/Project2015/model_versions/v10km/apply_model_v10km_inputs.R") 
  
  ptm <- proc.time()

  ### Calculate ALTITUDE effect on T - move elsewhere  ###
  elev.tdif<- (elevdif.block)*adiabatic.lapserate  # lapserate C per m difference
  
  # Define results array -  DIM of Results to match BLOCK SIZE
  max.hr<-(end.jd-start.jd+1)*24
  fileout<-paste(dir_finalt,"block20km-",sprintf("%03d",cellnum,sep=""),"-",year,sep="")
  print(fileout)
  #library(ff)
  #results <- ff(vmode="single",dim=c(ncol(dem.block),nrow(dem.block),max.hr),filename=fileout)
  results<-array(NA,c(ncol(dem.block),nrow(dem.block),max.hr))
  # Set index to track number of hours of data calculated
  hr.index<-0
  
  ##########################################################################################
  # APPLY MODEL to each day and hour 
  ########################################################################################## 
  for (jd in start.jd:end.jd) 
    {
      year<-DMYjd(jd)$year; month<-DMYjd(jd)$month ; day<-DMYjd(jd)$day ;# values for t
      print(paste(day,"/",month,"/",year))
      
      # Load day stack of direct and diffuse RADIATION and reproject to UK OS coordinates
      filein1<-paste(dir_dniday,"DNIhm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
      filein2<-paste(dir_sisday,"SIShm",year,sprintf("%02d",month,sep=""),sprintf("%02d",day,sep=""),".tif",sep="")
      dnr.24h<-stack(filein1) ; projection(dnr.24h)<-CRS("+init=epsg:4326")
      sis.24h<-stack(filein2) ; projection(sis.24h)<-CRS("+init=epsg:4326")
      dnr.24h<-projectRaster(dnr.24h,crs="+init=epsg:27700")
      sis.24h<-projectRaster(sis.24h,crs="+init=epsg:27700")
      
      # CALCULATE 5KM TEMPERATUE DATA AT HOURLY INTERVAL (or reassign for new day)
      if (jd==start.jd) { # load all t files
      t5km.prevday<-hourly_temperatures(jd-1,dir_temp,dir_hrtemp,grid5km.r) 
      t5km.day<-hourly_temperatures(jd,dir_temp,dir_hrtemp,grid5km.r) 
      t5km.nxtday<-hourly_temperatures(jd+1,dir_temp,dir_hrtemp,grid5km.r) 
      }
      if (jd!=start.jd) { 
        t5km.prevday<-t5km.day
        t5km.day<-t5km.nxtday
        t5km.nxtday<-hourly_temperatures(jd+1,dir_temp,dir_hrtemp,grid5km.r) 
      }
      
      # CALCULATE 5km RH data
      rh.day<-rh.hourly(jd,dir_rh,dir_rh5km,t5km.day,t5km.nxtday, grid5km.r,hourplot=FALSE)

      # Load and reproject CAL data to OS (downscaling by resampling takes place in hourly interval)
      cal.filein<-paste(dir_calday,"CALhm",year,sprintf("%02d",month,sep=""),sprintf("%02d.tif",day,sep=""),sep="")
      print(paste("CAL file in: ",cal.filein,sep=""))
      cal.day<-brick(cal.filein)
      cal.day<-projectRaster(cal.day,crs="+init=epsg:27700")
      
      # load sst data for day - resample to 100m for block and mask with dem.buffer to set land cells to NA
      sst.buffer<-sst.time.int(jd,dem.buffer,dir_sstm,dir_ssth,hr=12,msk=TRUE)
      #plot(sst.buffer)

      print("End of day processes ")  
  
      for (hr in 0:23)  #loop for hours
      {
        # Construct time label for plots 
        timelabel<-paste(DMYjd(jd)$day,"/",DMYjd(jd)$month,"/",DMYjd(jd)$year," ",hr,"h00",sep="")
        print(timelabel)
        
        hr.index<-hr.index+1 # advance index ready for next hour
        print(paste("hr.index=",hr.index))
        
        
        ### EXTRACT HOURLY BLOCK DATA ###
        # Extract reference temperature data for block and hour - tps using wider area to ensure some ref twemperature data available
        print("Temperature")
        tref.block<-tps.resample(crop(raster(t5km.day[,,hr+1],template=grid5km.r),dem.buffer),dem.block)
        if (jd==start.jd) prevtref.block<-tref.block
        
        # Interpolate Relative Humidity for hr and for block 
        rhref.block<-tps.resample(crop(raster(rh.day[,,hr+1],template=grid5km.r),dem.buffer),dem.block)
        
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
        # Calculate upwind sea temperature - record NA if no sea cell upwind within 10km
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
  
        # Save tref and tanomaly.100m for next time step
        prevtref.block<-tref.block
        prevt.block<-t.block
        prevtanom.block<-anom.block
        
        # Print summary statistics of downscaled temperature
        print(paste("Met Off 5km temperature:",cellStats(tref.block,mean),sep="") )
        print(paste("Mean 100m temperature:",cellStats(t.block,mean),sep="") )
        print(paste("Max 100m temperature:",cellStats(t.block,max),sep="") )
        print(paste("Min 100m temperature:",cellStats(t.block,min),sep="")  )                     
        
        # Save hourly result 
        print(paste("Number of hours of data calculated= ",hr.index,sep=""))
        results[,,hr.index] <- getValues(t.block,format="matrix")

        # PLOTS for TESTING
        plot(t.block,main=paste("Model T at",timelabel) )
        
      }# end loop for every hour
  } # end loop for every day
  
  # WRITE OUTPUT as 3D array of yearly file of hourly values
  save(results,file=fileout)
  #close(results)
  
  print(proc.time() - ptm)
  

