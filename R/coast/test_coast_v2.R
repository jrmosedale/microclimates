# Apply model to historic data
# v1 - loads most data at daily NOT hourly (v1) step 
# v2 - modification to Latent heat calc - use of 100m res reference data wihtout terrain effects - calc in radiation and wind downsclaing
# Conducts RH and LWR calculations within daily time step
##########################################################################################
# Define INPUTS - dates, 5km cell, parameters
# Most input = daily stack of hourly data
##########################################################################################

##########################################################################################
# Prepare INPUTS for ANY cell
##########################################################################################

# block.width<-5000 = ASSUMPTION
buffer<-20000 # already set in setup??

print("Calling set up...")
source("/home/ISAD/jm622/rscripts/setup_carson.R") # loads & runs setup file
print("Calling functions...")
source("/home/ISAD/jm622/rscripts/apply_model_v2_functions.R") # loads & runs functions
print("Calling inputs...")

### Load or set MODEL PARAMETERS ###
#parameters<-"/home/ISAD/jm622/rscripts/inputs/lizardparams.csv"
parameters<-"~/Documents/Exeter/Data2015/parameters/testparams.csv"
#parameters<-"F:/Data2015/parameters/lizardparams.csv"
params<-read.csv(parameters)
modify_params<-TRUE

if (modify_params==TRUE){   # change parameter values from those in file 
  params$estimates[1]<- 0.528 #intercept
  params$estimates[2]<- 0.137 #wc lapse
  params$estimates[3]<- 0.8  # shortwave rad^2.5
  params$estimates[4]<- 4.3   # longwave
  params$estimates[5]<- -0.43 # albedo
  params$estimates[6]<- -0.015 #windspeed
  params$estimates[7]<- -1.96 #inv wind speed
  params$estimates[8]<- -0.2 #ldif -0.2, -1.5
  params$estimates[9]<- 0.148 #sst-tref
  params$estimates[10]<- -15 #evapdif
  params$estimates[11]<- 0.06 # condendif
  params$estimates[12]<- -0.0012 # flow.acc
  params$estimates[13]<- -0.805 #tic
  params$estimates[14]<- 3.11 # albedo*rad
  params$estimates[15]<- -0.4 #rad*windspeed
  params$estimates[16]<- -0.59 # longwave*windspeed
  params$estimates[17]<- -0.233  #inv windspeed*ldif -2.2
  params$estimates[18]<-0.012 # wind x sst-tref 0.1
  params$estimates[19]<-0.33 # ldif*sst-tref
  params$estimates[20]<-0.0033 #flowacc*tic
}     
print("Parameters: ")
print(params$estimates)

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
cells<-c(682,683,684,715,716,717,741,742,743)
#for (x in cells) cell_location(gridmask.r,cells) # plot map of cell location

start.day <- 25
start.month<-6
start.year<-2012
print(paste("Start: ",start.day,"/",start.month,"/",start.year,sep=""))
end.day<-25
end.month<-6
end.year<-2012
print(paste("End: ",end.day,"/",end.month,"/",end.year,sep=""))

start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
print(start.jd);print(end.jd)

print(paste("UKCPCELL= ",ukcpcell,sep=""))

# Set cell and time range for testing
# start.day <- 2; start.month<-7 ; start.year<- 1992; hr<-0
# end.day<-3; end.month<-7; end.year<-1992
# ukcpcell<-926
# lizardcells<-c(916,917,918,919,924,925,926,927,928,929,930,931,932,934,935)
# dir_finalt<-"~/Documents/Exeter/Data2015/Temperature_100m/"

# jd<-start.jd

plot.var<-FALSE
plot.effects<-FALSE

##########################################################################################
# APPLY MODEL BY CELL ID to each day and hour 
########################################################################################## 

for (ukcpcell in cells){
  
  #source("/home/ISAD/jm622/rscripts/apply_model_v2_inputs_macuse.R") # defines inputs,parameters, etc
  source("~/Documents/Exeter/Rprojects/Project2015/model_versions/v2/apply_model_v2_inputs_macuse.R") # defines inputs,parameters, etc
  
  #plot(dem.block)
  
  # Set index to track number of hours of data calculated
  hr.result<-1
  
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
    tref.day<-stack() # for TESTING ONLY
    tmodel.day<-stack() # for TESTING ONLY
    radef.day<-stack() # for TESTING ONLY
    latef.day<-stack() # for TESTING ONLY
    coastef.day<-stack() # for TESTING ONLY
    elevef.day<-stack() # for TESTING ONLY
    invwstr.day<-stack() # for TESTING ONLY
    anom.day<-stack() # for TESTING ONLY
    ldif.day<-stack() # for TESTING ONLY
    sstdif.day<-stack() # for TESTING ONLY
    wdir.day<-stack() # for TESTING ONLY
    
    # Load prev day temperature file and store last hour of data for block
    prevt.filein<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd-1)$year, "-",sprintf("%02d",DMYjd(jd-1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd-1)$day,sep=""),"_100m.r", sep="")
    load(prevt.filein) #t100m.day
    prevtref.block<-crop(raster(t100m.day[,,23],template=dem),dem.block) # sets prev tref to 23h00 from previous day to jd
    
    # Load next day's  temperature file and store first hour (for block)cropped later for block??)
    nextt.filein<-paste(dir_hrtemp,"HrTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_100m.r", sep="")
    load(nextt.filein) #t100m.day
    #nexttref.r<-raster(t100m.day[,,1],template=dem) # sets next tref to 23h00 from previous day to jd (uncropped for RH calculation)
    nexttref.block<-crop(raster(t100m.day[,,1],template=dem),dem.block)
    
    # Load daily files of hourly values of temperature, RH, lwr
    t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),"_100m.r", sep="")
    load(t.filein) #t100m.day
    
    # load sst data for day - resample to 100m for block and mask with dem.buffer to set land cells to NA
    infile.sst<- paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")
    sst5km.r<-raster(infile.sst)  
    sst.buffer<-resample(sst5km.r,dem.buffer)
    sst.buffer<-mask(sst.buffer,dem.buffer, inverse=TRUE) 
    
    print("End of day processes ")  
    #print(proc.time()[1:3] - ptm[1:3])
    
    for (hr in 0:13)  #loop for hours
    {
      # Construct time label for plots 
      timelabel<-paste(DMYjd(jd)$day,"/",DMYjd(jd)$month,"/",DMYjd(jd)$year," ",hr,"h00",sep="")
      print(timelabel)
      par(mfrow=c(2,3))
      
      
      ### EXTRACT HOURLY BLOCK DATA ###
      # Extract reference temperature data for block and hour
      print("Temperature")
      tref.block<-crop(raster(t100m.day[,,hr+1],template=dem),dem.block)
      plot(tref.block,main="Tref")
     
      ### DOWNSCALE WIND to direction , strength and inverse wind strength for block - Prog: Wind_downscale_blocks    
      # Output: = writes raster for wind str, wind dir and inv wind str (true=print results)      print("Wind")
      wind.results<-wind.tpsdownscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval,print.res=FALSE,write.files=FALSE )
      wdir.block<-wind.results[[1]]
      wstr.block<-wind.results[[2]]
      invwstr.block<-wind.results[[3]]  
      refwstr.block<-wind.results[[4]]  
      #plot(wdir.block,main="wdir")
      
      ### CALC COASTAL EFFECT - NEW!!!!! ###
      print("Coastal")
      # Calculate upwind.sst & sst-tref for block - Prog: upwind_sst_blocks
      # Calculate upwind sea temperature - record NA if no sea cell upwind within 20km
      upwind.results<-upwind.sst.block(wdir.block,sst.buffer,dem.buffer,dem.block)
      upwindsst.block<-upwind.results[[1]]
      plot(upwindsst.block,main="upwind.block")
      seadist.block<-upwind.results[[2]]
      plot(seadist.block,main="seadist.block")
      sst.dif<-sst.tref(upwindsst.block,tref.block,seadist.block)  
      plot (sst.dif,main="sst.dif")
      
      # Calculate ldif according to wind direction
      ldif.block<-calc.ldif.block(ldif.stack,wdir.block,interval)    
      plot(ldif.block,main="ldif.block")
     
      ####################################################################################
      # Apply model paramters to calculate temperature anomaly
      ####################################################################################
      #t.block<-tref.block+elev.effect+(rad.effect*2)+(latent.effect*0.3)+(coast.effect/2)
      print("Calculating 100m temperatures...")
     
      CoastEffect<-params$estimates[8]*ldif.block + 
        params$estimates[9]*sst.dif + # set to zero effect
        params$estimates[17]*invwstr.block*ldif.block+
        params$estimates[18]*wstr.block*sst.dif
      params$estimates[19]*ldif.block*sst.dif
    
      plot (CoastEffect,main="Coast Effect")
      
      ####################################################################################
      #### Plot results ####
      ####################################################################################
      par(mfrow=c(2,3))
      if (plot.var==TRUE){
        plot(t.block,main=paste("T100m ",ukcpcell," at ",timelabel,sep=""))
        #plot(t.block-prevt.block,main="T hourly change")
        #plot(anom.block,main="Tanom")
        plot(tref.block,main="Tref")
        #plot(total.block,main=paste("Total Rad ",timelabel,sep=""))
        ##plot(direct.block,main=paste("Direct ",sep=""))
        #plot(diffuse.block,main=paste("Diffuse ",sep=""))
        #plot(rhref.block,main=paste("RH ref ",timelabel,sep=""))
        #plot(rh.5km,main=paste("RH 5km ",timelabel,sep=""))
        #plot(rh.block,main=paste("RH 100m ",timelabel,sep=""))
        #plot(CRE.5km,main=paste("Evap Ref ",timelabel,sep=""))
        #plot(evapdif,main=paste("Evap difference  ",sep=""))      
        #plot(wc.5km,main=paste("Cond Ref ",timelabel,sep=""))
        #plot(condendif,main=paste("Cond difference ",sep=""))
        #plot(crop(aspect.buffer,dem.block),main=paste("Aspect ",timelabel,sep=""))
        plot(wstr.block,main="Wstr")
        #plot(wdir.block,main="Wdir")
        #plot(p100.block,main=paste("P 100m ",timelabel,sep=""))
        plot(lwr.block,main=paste("LWR ",sep=""))
        #plot(cal.block,main="CAL")
        #plot(flow.block*tic,main="flow*tic")
        plot(ldif.block,main="Ldif")
        plot(sst.dif,main="SST dif")
        #plot(invwstr.block*ldif.block,main="invwstr x ldif")
        #plot(wstr.block*sst.dif,main="wstr x sst.dif")
        #plot(ldif.block*sst.dif,main="ldif x sst .dif")
      }
      if (plot.effects==TRUE){
        plot(t.block,main=paste("T100m ",ukcpcell," at ",timelabel,sep=""))
        #plot(t.block-prevt.block,main="T hourly change")
        #plot(tref.block,main="Tref block")
        plot(anom.block,main="Tanom") 
        
        #plot(wstr.block,main="Wstr") 
        #plot(wdir.block,main="Wind Dir") 
        plot(RadEffect,main="Radiation Effect ")
        plot(LatentEffect,main="Latent Heat Effect")
        plot(params$estimates[2]*wc.100m,main="WaterCond from cooling effect")
        #plot(params$estimates[10]*evapdif,main="Evapdif effect")
        plot(params$estimates[11]*condendif,main="Condendif effect")
        
        plot(CoastEffect,main="Coast Effect ")
        #plot(params$estimates[8]*ldif.block,main="ldif effect") 
        #plot(params$estimates[17]*invwstr.block*ldif.block,main="invwstrxldif effect")
        #plot(params$estimates[18]*wstr.block*sst.dif,main="wstr x sst effect")
        #plot(params$estimates[19]*ldif.block*sst.dif,main="ldif x sstdif effect")
        
        plot(ElevEffect,main="Elev Effect ")
        plot(FlowEffect,main="Flow Effect")
      }  
      
      # OUTPUTS: SAVE DOWNSCALED TEMP FOR DEM.BLOCK for every Hour 
      tref.day<-stack(tref.day,tref.block)
      coastef.day<-stack(coastef.day,CoastEffect)
      invwstr.day<-stack(invwstr.day,invwstr.block)
      ldif.day<-stack(ldif.day,ldif.block)
      sstdif.day<-stack(sstdif.day,sst.dif)
      wdir.day<-stack(wdir.day,wdir.block)
      
      
      
      #if (ukcpcell%%10==0){plot(dem.block, main=paste("UKCP09 cell= ",ukcpcell,sep=""))}
      # Save tref and tanomaly.100m for next time step
      prevtref.block<-tref.block

      hr.result<-hr.result+1 # advance index ready for next hour
      
    }# end loop for every hour
    
    # FOR TESTING only
    fileout<-paste(dir_finalt,"tref-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(tref.day,fileout,format="GTiff",overwrite=TRUE)

    fileout<-paste(dir_finalt,"coastef2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(coastef.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"invwstr2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(wstr.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"ldif2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(ldif.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"sstdif2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(sstdif.day,fileout,format="GTiff",overwrite=TRUE)
    
    fileout<-paste(dir_finalt,"wdir2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    writeRaster(wdir.day,fileout,format="GTiff",overwrite=TRUE)
    
    
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

par(mfrow=c(2,4))
map.r<-raster()

# 2. Print hourly Temperature maps for 24hr
for (hr in 0:23){
  blocks0<-vector("list",numcells)
  blocks1<-vector("list",numcells)
  blocks2<-vector("list",numcells)
  blocks3<-vector("list",numcells)
  blocks4<-vector("list",numcells)
 # blocks5<-vector("list",numcells)
#  blocks6<-vector("list",numcells)
  
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
    
    infile<-fileout<-paste(dir_finalt,"tref-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks0[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"coastef2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks1[[i]]<-raster(stack(infile),layer=hr+1)
    
    fileout<-paste(dir_finalt,"invwstr2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks2[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"wdir2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks3[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"sstdif2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks4[[i]]<-raster(stack(infile),layer=hr+1)
    
    infile<-paste(dir_finalt,"ldif2-",ukcpcell,"-",sprintf("%02d",DMYjd(jd)$day,sep=""),"-",sprintf("%02d",DMYjd(jd)$month,sep=""),"-",year,".tif",sep="")
    print(infile)
    blocks5[[i]]<-raster(stack(infile),layer=hr+1)
    
    i<-i+1
  }
  
  ### PLOT merged blocks
  map.r<-do.call(merge, blocks0)
  titletext<-paste("Tref at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks1)
  titletext<-paste("Coast at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks2)
  titletext<-paste("InvWstr at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks3)
  titletext<-paste("Wdir at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks4)
  titletext<-paste("SST at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
  map.r<-do.call(merge, blocks5)
  titletext<-paste("LDIF at ",hr,"hr",sep="")
  plot(map.r,main=titletext)
  
} # end 24 hr
