# Libraries
library(R.utils) # IMPORTANT - hides raster functions requiring raster::
library(ncdf4)
library(raster)
library(rgdal)
library(sp)
library(chron)
library(insol) # required for julian day functions
library(mgcv) # require package for inputation (of CAL)
library(fields) # required for thin plate spline

# Directories
# Templates - 5km cells
dir_grids<-"~/Documents/Exeter/Data2015/Templates/"

# Terrain
dir_terrain<-"~/Documents/Exeter/Data2015/terrain/"
  
# Temperature
dir_zip<-"~/Documents/Exeter/Data2015/Temp5km/zip/"
dir_temp<-"~/Documents/Exeter/Data2015/Temp5km/extract/"
dir_hrtemp<-"~/Documents/Exeter/Data2015/Temp5km/hourly/" # dir for output files

# Wind
dir_wind<-"~/Documents/Exeter/Data2015/Wind/"
dir_windstrength<-"~/Documents/Exeter/Data2015/Wind/strength/"
dir_winddirection<-"~/Documents/Exeter/Data2015/Wind/direction/"
dir_windinvstr<-"~/Documents/Exeter/Data2015/Wind/invstr/"
dir_shelter<-"~/Documents/Exeter/Data2015/Wind/Shelter/"

# Coast - landsea ratio and %land
dir_percland<-"~/Documents/Exeter/Data2015/CoastEffect/percland/"
dir_lsratio<-"~/Documents/Exeter/Data2015/CoastEffect/lsratio/"
dir_ldif<-"~/Documents/Exeter/Data2015/CoastEffect/ldif/"

# SST
dir_sst<-"~/Documents/Exeter/Data2015/sst/"
#gunzip(paste(dir_sst,"HadISST_sst.nc.gz",sep=""))
dir_sstm<-paste(dir_sst,"monthly/",sep="")
dir_ssth<-paste(dir_sst,"hourly/",sep="")

dir_upwindsea<-"~/Documents/Exeter/Data2015/CoastEffect/upwindmaps/"

# Radiation budget
dir_dnitar<-"~/Documents/Exeter/Data2015/CMSAF-DNI/tar/"
dir_dnigz<-"/Users/jonathanmosedale/Documents/Exeter/Data2015/CMSAF-DNI/ncgz/"
dir_sistar<-"~/Documents/Exeter/Data2015/CMSAF-SIS/tar/"
dir_sisgz<-"/Users/jonathanmosedale/Documents/Exeter/Data2015/CMSAF-SIS/ncgz/"
dir_caltar<-"~/Documents/Exeter/Data2015/CMSAF-CAL/tar/"
dir_calgz<-"/Users/jonathanmosedale/Documents/Exeter/Data2015/CMSAF-CAL/ncgz/"

dir_sis<-"~/Documents/Exeter/Data2015/CMSAF-SIS/extract/"
dir_dni<-"~/Documents/Exeter/Data2015/CMSAF-DNI/extract/"
dir_cal<-"~/Documents/Exeter/Data2015/CMSAF-CAL/extract/"
dir_calimp<-"~/Documents/Exeter/Data2015/CMSAF-CAL/imputed/"
dir_lwr<-"~/Documents/Exeter/Data2015/lwr/"

dir_rad<-"~/Documents/Exeter/Data2015/radiation/"
dir_direct<-"~/Documents/Exeter/Data2015/radiation/rasters/direct/"
dir_diffuse<-"~/Documents/Exeter/Data2015/radiation/rasters/diffuse/"
dir_total<-"~/Documents/Exeter/Data2015/radiation/rasters/total/"

# Sea Pressure
dir_pressure025<-"~/Documents/Exeter/Data2015/SeaPressure/Pdaily025/"
dir_pressure<-"~/Documents/Exeter/Data2015/SeaPressure/"

# Relative Humidity
dir_rh<-"~/Documents/Exeter/Data2015/RelHumidity/"
dir_rh5km<-"~/Documents/Exeter/Data2015/RelHumidity/rh5km/"

# Albedo
dir_albedo<-"~/Documents/Exeter/Data2015/albedo/"

# CONSTANTS
W.to.MJhr<-0.0036 # converting CMSAF satellite rad values (W/m2) to MJ/m2/hour required for calcs
albedo<-0.2

##########################################################################################
# Define Geographical extent of interest 
##########################################################################################
latlong <- "+init=epsg:4326"
ukgrid <- "+init=epsg:27700"
# Define area of interest (no buffer)
#e.dem<-extent(c(70000,420000,10000,180000 )) # includes scilly isles
e.dem<-extent(c( 120000,420000,10000,180000 )) # excludes scilly isles

# Define 100m dem rasters
demuk<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif", crs=(ukgrid))
e.ukexp<-c(0,7e+05,-10000,1200000) # expand to allow 20km buffer to south of area of interest - set to sea (NA)
demuk<-extend(demuk,e.ukexp,values=NA)

dem<-crop(demuk,e.dem)

# define 20km buffered area around region of interest
buffer<-20000
e.buf20km<-extent(xmin(dem)-buffer,xmax(dem)+buffer,ymin(dem)-buffer,ymax(dem)+buffer)# Run setup programs for creating constant raster maps etc
dembuf<-crop(demuk,e.buf20km)

# Define 5km grid rasters
grid5kmuk.r<-raster(paste(dir_grids,"ukhistmask.grd",sep="")) #  1=valid cell, NA = sea or not data
grid5km.r<-crop(grid5kmuk.r,e.dem)
grid5kmbuf.r<-crop(grid5kmuk.r,e.buf20km)

##########################################################################################
# Initial analysis and creation of required files - to be run only ONCE
# Progs:  Wind_sheltermaps_functions
#         percent_land_maps_function
#         inv_lsratio_maps_function
#         ldif_maps_function
#         elevdif_map_function
#         wind_1_readdata
# Calculate dem derived maps etc

# Terrain maps - used by radprog_blocks
slope<-terrain(dembuf, opt='slope', unit='degrees')
aspect<-terrain(dembuf, opt='aspect', unit='degrees')
plot(aspect,main="Aspect")
plot(slope,main="Slope")
writeRaster(slope,file=paste(dir_terrain,"slope.tif",sep=""),overwrite=TRUE)
writeRaster(aspect,file=paste(dir_terrain,"aspect.tif",sep=""),overwrite=TRUE)

# Create wind shelter maps - Prog: Wind_sheltermaps_functions
# Input: requires 20km buffer area. Output tif files saved to dir_shelter
interval<-10 # = division of wind direction in degrees
wind_sheltermaps(demuk,e.dem,interval,dir_shelter)

# Calculate % land in radius of each cell used in historic temp calculations - Prog: percent_land_maps_function (NB: will differ in UKCP09WG analysis)
# Input: requires buffer = radius. Output 5km and 100m grid tif files saved to dir_percland
radius<-10000 #10km or20km
percent_land_maps(dembuf,grid5km.r,e.dem,radius,dir_percland) 

# inv land:sea ratio by wind dir - Prog: inv_lsratio_maps_function
# Input: requires 10km buffer around area of interest
inv.lsmaps(dembuf,dem,interval,dir_lsratio)

# %landref-inv.land.sea by wind.dir (OR Lref-Lcell) - Prog: ldif_maps_function
ldif.maps(interval,radius,dir_percland,dir_lsratio,dir_ldif)

# elevation dif maps using mean of 5km reference cells (central xy can be missing) - Prog: elevdif_map_function
#For each 100m grid cell interpolate vs nearest cells and as function of elevation difference  between 100m elevation and each 5km elevation
elevation.dif.map(dem,grid5km.r,dir_terrain)
  
# Other programs to be run only once (eg unzip, download, basic provessing of data sources)
# Sea Pressure data - daily 0.25 deg - load and crop raster bands to area of interest
gzfile<-paste(dir_pressure025,"pp_0.25deg_reg_v11.0.nc.gz",sep="")
ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")
#gunzip(filename=gzfile, destname=ncfile, overwrite=TRUE)


# Wind_1_readdata - Write single files for wind u and v (all times) - Prog: wind_1_readdata

# Extract Radiation datafiles. Prog: extract_tar_ncdf_function
# WARNING - currently extracts ALL tar files - produces 8750 files per year per factor
extract.tar.to.ncdf(dir_dnitar,dir_dnigz,dir_dni)
extract.tar.to.ncdf(dir_sistar,dir_sisgz,dir_sis)
extract.tar.to.ncdf(dir_caltar,dir_calncgz,dir_calnc)

### load WIND data
# data are arrays of easterly and northerly wind components at spatial resolution of 2.5 degrees and temporal resolution of 4x daily
# data automatically assigned name when written out: wind_u and wind_v
load(file=paste(dir_wind,"wind_u.r",sep=""))
load(file=paste(dir_wind,"wind_v.r",sep=""))




##########################################################################################
# Run processes specific for chosen time period 
# Progs:t5km_to_hourly_blocks
#       sst_downscale_function
#       
##########################################################################################
# Set time range for which data will be analysed
start.year<-1992
start.month<-7
start.day<-1
hr<-0
end.year<-1992
end.month<-7
end.day<-10

### Unzip TEMPERATURE files required - start.jd -1?? - Prog: JDfunctions
start.jd<-JDdmy(start.day,start.month,start.year) 
end.jd<-JDdmy(end.day,end.month,end.year)
unzip_yearfiles(start.jd,end.jd,dir_zip,dir_temp)



### Downscale SST to daily 5km grid - Prog: sst_downscale_function ###

# IMPORTANT: NEEDS TO BE FOR 20KM BUFFER REGION TO ALLOW CALC OF UPWIND SST
in.file<-paste(dir_sst,"HadISST_sst.nc",sep="") # in.file required for sst.spdownsc
ncfile<-nc_open(paste(dir_sst,"HadISST_sst.nc",sep="")) # summary of file variables,  dimensions, attributes

# 1. WRITE monthly files for buffered uk region - normally 10km buffer - CHECK START DATE _ NEED MONTH BEFORE START!!!
# Calculate jd for start of PREVIOUS month to start.jd *
# ?? Uses land mask of grid5km.r
if (start.month==1) { 
  year1<-start.year-1
  month1<-12 } else {
  year1<-start.year
  month1<-start.month-1 }
if (end.month==12) { 
  year2<-end.year+1
  month2<-1 } else {
    year2<-end.year
    month2<-end.month+1 }

sst.spdownsc(year1,month1,year2,month2,grid5kmbuf.r,dir_sstm)  
 
# 2. WRITE daily files 
for (jd in start.jd:end.jd) {
  sst.time.int(jd, dir_sstm, dir_ssth)
} # end day

# 3. delete monthly files
for(y in year1:year2){
  for(m in month1:month2){
    filename<-paste(dir_sstm,"sst_",y,"_",m,".tif",sep="")
    print(fileout)
    file.remove(filename)
 }}
nc_close(ncfile)
###   END  ###

## Downscale daily TEMPERATURE data to hourly  Prog: t5km_to_hrmatrix ###
# WRITES t5km.day dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),".r"
# REQUIRES: elevdif_map FUNCTIONS to fill in values for missing 5km cells with LAND cells at 100m
# includes data for previous day to start.day
hourly_temperatures(start.jd-1,end.jd,dir_temp,dir_hrtemp,grid5km.r) 

#Downscale REL HUMIDITY to hourly 5km Prog: rel.hum.v2
# Writes rh.day as dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R"
rh.hourly(start.jd,end.jd,dir_rh,dir_rh5km,dir_hrtemp,grid5km.r)

# Downscale and impute CLOUD ALBEDO - WARNING - LONG TIME - RUN SEPERATELY??  Prog: calc.cal.hrly
# WRITES calimp.day as dir_calimp,"CALimp_5km_",year,"_",month,"_",day,".R
cal.5km.impute(dir_cal,dir_calimp,start.jd,end.jd,grid5km)

# Calculate LONG WAVE RADITAION at 5km hourly res ffrom CAL, RH T  Prog: longwav_grids
write_lwr_dayfies(start.jd,end.jd,grid5km.r)

# Pressure - already unpacked
p.ncfile<-paste(dir_pressure,"pp_0.25deg_reg_v11.0.nc",sep="")

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

# block.width<-5000 = ASSUMPTION
buffer<-20000 # set for max required for any individual program (coast effect???)

# LOOK for every block...
ukcpcell<-574 # for testing only
for (ukcpcell in(1:length(landcells)))
  {
    x<-landcells[ukcpcell,1]
    y<-landcells[ukcpcell,2]
    e.block<-extent(x-2500,x+2500,y-2500,y+2500)
    dem.block<-crop(demuk,e.block)
    e.buffer<-extent(x-2500-buffer,x+2500+buffer,y-2500-buffer,y+2500+buffer)
    dem.buffer<-crop(demuk,e.buffer)
    plot(dem.block)
    
    # Perform operations for each block but fixed for all timeperiods - ie CROPPING Terrain datasets
    # Wind - load shelter matrix for block - Prog: wind_downscale_blocks
    shelter.block<-block.sheltermap(dem.block,dir_shelter)
    
    # Slope & aspect cropped to dem.buffer for use in radiation_downscale
    slope.buffer<-crop(slope,dem.buffer)
    aspect.buffer<-crop(aspect,dem.buffer)
    
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
    
    # LOOP FOR EACH TIME PERIOD from START to END DATE - how to loop using jd in hour steps, starting at 12.00??
    for (t in start.jd:end.jd by ??? ) 
      {
        year<-DMYjd(t)$year; month<-DMYjd(t)$month ; day<-DMYjd(t)$day ;# hr<-0

        # Load daily files of hourly values of temperture, RH, lwr
        rh.filein<-paste(dir_rh5km,"RH_5km_",year,"_",month,"_",day,".R",sep="")
        load(rh.filein) # rh.day
        
        prevt.filein<-paste(dir_rh5km,"RH_5km_",DMYjd(t-1)$year,"_",DMYjd(t-1)$month,"_",DMYjd(t-1)$day,".R",sep="")
        load(prevt.filein)
        t5km.prevday<-t5km.day # prev day data
        # Store last hour of temperature data from previous day before loading new files
        prev.tref5km.r<-raster(t5km.prevday[,,23],template=grid5km.r)
        
        t.filein<-paste(dir_hrtemp,"HrTemp_", year, "-",sprintf("%02d",month,sep=""),"-", sprintf("%02d",day,sep=""),".r", sep="")
        load(t.filein) #t5km.day
       
        lwr.filein<-paste(dir_lwr,"lwr_",year,"_",month,"_",day,".R",sep="")
        load(lwr.filein) # lwr.day
        
        for (hr in 0:23)  #loop for hours???
        {
            ### REF TEMPERATURES ###
            # Extract 5km reference temperature data  for block at 100m resolution - Prog: ref5km.to.block100m
            tref5km.r<-raster(t5km.day[,,hr+1],template=grid5km.r)
            tref.block<-ref5km.to.block100m(dem.block,tref5km.r)
            
            if (hr>0){prev.tref5km.r<-raster(t5km.day[,,hr],template=grid5km.r)} # else if hr=0 prev.tref5km already set to 23h00 of previous day
            prevtref.block<-ref5km.to.block100m(dem.block,prev.tref5km.r)
            
            ### CALC RADIATION EFFECT ###
            # Downscale wind to direction , strength and inverse wind strength for block - Prog: Wind_downscale_blocks    
            # Output: = writes raster for wind str, wind dir and inv wind str (true=print results)
            wind.results<-wind.downscale(day,month,year,hr,dem.block,wind_u,wind_v,dir_wind,dir_shelter,shelter.block,interval=10,print.res=TRUE,write.files=TRUE )
            wdir.block<-wind.results[[1]]
            wstr.block<-wind.results[[2]]
            invwstr.block<-wind.results[[3]]
            
            #invwstr.dif<-invwstr.ref-invwstr.block
            
            # Ref wind strength = mean for 5km block
            #wstr.anom.block<-wstr.block-cellStats(wstr.block, stat='mean', na.rm=TRUE)
            wstr.mean<-matrix(cellStats(wstr.block, stat='mean', na.rm=TRUE),nrow=nrow(wstr.block),ncol=ncol(wstr.block) )
            wstr.mean.block<-raster(wstr.mean,template=wstr.block)
            wstr.mean.block<-mask(wstr.mean.block,wstr.block)
        
            # Radiation downscale for block - Prog: radprog_blocks - #output  = list(direct.r,diffuse.r,total.r)  
            # Outputs converted to MJ/m2 from W/m2
            rad.results<-radiation_downscale(day,month,year,hr,demuk,dem.buffer,dem.block,slope.buffer,aspect.buffer,dir_rad,print.results=TRUE,write.files=FALSE )
            direct.block<-rad.results[[1]]*W.to.MJhr
            diffuse.block<-rad.results[[2]]*W.to.MJhr
            total.block<-rad.results[[3]]*W.to.MJhr # = total incoming shortwave radiation - IGNORING ALBEDO EFFECT
            
            # Extract lwr for hour and block at 100m - effective res of 5km
            lwr5km.r<-raster(lwr.day[,,hr+1],template=grid5km.r)
            lwr.block<-ref5km.to.block100m(dem.block,lwr5km.r)
            
            # Calculate net radiation for block - from incoming shortwave, albedo and outgoing longwave radiation
            netr.block<-(total.block*(1-albedo))-lwr.block
            
            # Ref netr = mean for 5km block
            netr.mean<-matrix(cellStats(netr.block, stat='mean', na.rm=TRUE),nrow=nrow(netr.block),ncol=ncol(netr.block) )
            netr.mean.block<-raster(netr.mean,template=netr.block)
            netr.mean.block<-mask(netr.mean.block,netr.block)
            
            
            ### CALC COASTAL EFFECT ###
            # If new day then extract 5km SST data for buffer region at 100m 
            if (paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")!=infile.sst){ 
                infile.sst<- paste(dir_ssth,"sst_",year,"_",month,"_",day,"_12h.tif",sep="")
                sst5km.r<-raster(infile.sst)  
                sst.buffer<-resample(sst5km.r,dem.buffer)
             }
         
            # Calculate upwind.sst & sst-tref for block - Prog: upwind_sst_blocks
            # First mask with dem.buffer to set land cells to NA
            sst.buffer<-mask(sst.buffer,dem.buffer) # 
            sst.block<-upwind.sst.block(wdir.block,sst.buffer,dem.buffer,dem.block)
            sst.dif<-sst.tref(sst.block,tref.block)  
            #plot (sstdif.block)
            
            # Calculate ldif according to wind direction
            ldif.block<-calc.ldif.block(ldif.stack,wdir.block,interval)
            
            
            ### LATENT HEAT EFFECT CALCS ###
            # record whether day or night - used in latent heat calculations
            if (cellStats(total.block,max)>0) dn<-1 else dn<-0
            
            # Extract Relative Humidity for hr and for block 
            rh5km.r<-raster(rh.day[,,hr+1],template=grid5km.r)
            rhref.block<-ref5km.to.block100m(dem.block,rh5km.r)
            
            # Calculate sea level Pressure for 5km cell at 100m res Prog: prepare.pressurre   - CHECK!!! MORE!
            pref.block<-downscale.pressure(dem.block,p.ncfile,t,write.file=TRUE) # = Pressure at sea level
            # p.sw<-downscale.pressure(ncfile,grid5km.r,t,write.file=TRUE)
            # p.blk<-crop(p.sw,dem.block)
            # Correct for elevation variation within 5km cell 
            p100.block<-correct.pressure(pref.block,tref.block,elevdif.block)   # but what of where no tref 5km cell - find nearest??
            
            # Calculate latent heat difference Prog latentheat.block.functions        
            # Calculate proxy temperature at 100m from tref (5km) and previous anomaly between 5km and 100m
            if (t==start.jd) tproxy.block<-tref.block else tproxy.block<-tref.block+prev.anom.block
            # If first t set previous temp at 100m to proxytemp ??? 
            if (t==start.jd) prevt.block<-tproxy.block
            
            # Calculate relative humdidity at 100m res
            rh.block<-rh.change(tref.block,tproxy.block,rhref.block) 
            
            # Calcute difference in evapotranspiration - CRE function
            CRE.5km<-CRE(tref.block,netr.mean.block,rhref.block,pref.block,dn,wstr.mean.block)
            CRE.100m<-CRE(tproxy.block,netr.block,rh.block,p100.block,dn,wstr.block)
            evapdif<-CRE.5km-CRE.100m # regress this against temperature in model
            
            # Calculate difference in water condensation - water.conden function
            wc.5km<-Water.conden(prevtref.block,tref.block,rhref.block)
            wc.100m<-Water.conden(prevt.block,tproxy.block,rh.block)
            condendif<-wc.5km-wc.100m   # regress this against temperature in model
       
            
            # Downscale temperature to 100m 
            # t.block<-t.to100m(tref.block,evapdif,condendif,elevdif, )
            
            elev.effect<- (elevdif.block/100)*9.8  # altitude effect
            rad.effect<-netr.block+wstr.block+(wstr.block*netr.block) # radiation effect
            coast.effect<-invwstr.block * ldif.block + (wstr.block*sst.dif) + (ldif.block*sst.dif) # coastal effect
            latent.effect<-evapdif+condendif  # latent heat effect
            
            # OUTPUTS: SAVE DOWNSCALED TEMP FOR DEM.BLOCK for every HOUR 
        
            #if (ukcpcell%%10==0){plot(dem.block, main=paste("UKCP09 cell= ",ukcpcell,sep=""))}
           
            # Delete downscaled block files except final hourly temperature files
            
            # Save tref and tanomaly.100m for next time step
            prevtref.block<-tref.block
            prevt.block<-t.block
            prev.anom.block<-t.block-tref.block # Check correct order/sign
            
        }# end loop for every hour
      
      newt.fileout<-paste(dir_newt,"block_",ukcpcell,"_on_",day,"_",month,"_",year,".tif")
      
    } # end loop for eery day

} # end for every block

##########################################################################################
# Re-load and stitch together block results if required using Mosaic etc


##########################################################################################
# Remove unwanted files specific to time period

##########################################################################################
# Calculate annual risks for time period for every 100m cell

# Fixed phenology values
# Phenology models - to give values by year
# GDD variants & new options (dif start/end, dif res, Hughlin index etc) by  year
# Frost risk - freq, index, ...
# Flowering temperature
# Humidity / Leaf wetness / Disease risk index - OPTIONS
# Yield model? to compare historical predictions with actual?
# Water avail?

# Timeperiods: every cell: historic (30 years); baseline (30x1000); 2020(?), 2050 mid emissions (30x1000)
# For baseline and historic will have to run by cell/block and incorporate with temp downscaling loop to be efficient

# Save and reference as single raster for whole area for each risk:
#                       (i) each historic year risk data - predicted number of frost days, gdd, etc 
#                       (ii) historic risk frequency - average of whole period ?? Why ?? - or the option to calculate for defined time period
#                       (ii) probability (mean) of risk events under different scenarios/timeperiods

##########################################################################################

# Loading of test files
# Input wind dir data - 100m hrly for block 
infile.wind<-paste(dir_winddirection,"direction_",year,"_",month,"_",day,"_",hr,".tif",sep="")
#infile.wind<-paste(dir_winddirection,"direction_2010_1_2010_12.tif",sep="")
wdir.r<-raster(infile.wind)

# Delete unzipped historic 5km temperature files used
for (jd in start.jd:end.jd) {
  max.infile<-paste(dir_temp,"MaxTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="")
  min.infile<-paste(dir_temp,"MinTemp_", DMYjd(jd+1)$year, "-",sprintf("%02d",DMYjd(jd+1)$month,sep=""),"-", sprintf("%02d",DMYjd(jd+1)$day,sep=""),"_ACTUAL.txt", sep="") 
  file.remove(max.infile)
  file.remove(min.infile)
}


# ALTERNATIVE - NOT USING UKCP09 cells - for all cells including sea

block.width<-5000
buffer<-20000 # set for max required for any individual program (coast effect???)

for (demx in (seq(xmin(dem),(xmax(dem)-block.width),block.width)) ){
  for (demy in (seq(ymin(dem),(ymax(dem)-block.width),block.width)) ){
    e.block<-extent(demx,demx+block.width,demy,demy+block.width)
    dem.block<-crop(demuk,e.block)
    e.buffer<-extent(demx-buffer,demx+block.width+buffer,demy-buffer,demy+block.width+buffer)
    dem.buffer<-crop(demuk,e.buffer)
    # LOOP FOR EACH TIME PERIOD from START to END DATE
  
    # CALCULATIONS HERE 
    # OUTPUTS = DOWNSCALED TEMP FOR DEM.BLOCK AT TIME T - SAVE
    if (demy==ymax(dem)-block.width){plot(dem.block)}
    if (demy==ymax(dem)-block.width){plot(dem.buffer)}
  }   
}
# LOAD and MOSAIC individual blocks BY TIME PERIOD

e.test<-extent(c( 175000,180000,15000,20000 )) 
demtest<-crop(demuk,e.test)
plot(demtest)

e.test.buf<-extent(c( 145000,210000,12000,23000 )) 
demtestbuf<-crop(demuk,e.test.buf)
plot(demtestbuf)
