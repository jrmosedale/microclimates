###############################################################################
# UKCP09 ncdf files

###############################################################################

library(ncdf4)
library(raster)

dir_ukcp<-"C:/Data2015/UKCPWG/"

#############################################################
# UKCP job download 

rootlink<-"http://ukclimateprojections-app.metoffice.gov.uk/wps/dl/"

jobs<-c("475414440391872179849uQ","550714440487246037297pG","137514440640759367404jK",
          "952414440728442878183oW","283414441199510635653dG","442014441280932918845sC",
        "479614441457578080566aY","460114441560593752433dM","721814441650139306939nM") # next batch of jobs from 05/10
        
     # to start 07/10 exeter only 

#jobs<-c("108614439773176275491iN","284614439553739535050eH","233514439063789375164oV",
        "576714438688319301186oN","176214438220575502115eH","673314438078091858227eZ",
         "459414437876413956455oA","529414437792178148950mE","157214437336530708216dN",
         "612114437286483733447tY","881714437196666864143cO","344914437110116763614hR",
         "455914436975255494638fP","694714436882382887336vV","268714436474681300976vG",
         "052714436186687635412vB","395414436115794723955qZ","468914435617092874675nX",
          "964314435471234207780uT","843514435234502113340mP","925914435143853672251aG",
         "196314434741817930694lP","014014434641642655445lR","527914434576130631168lI",
        "494614433070302715908tM","657914432970588913779zU","218214432698090514861gY",
         "489614432584946306544pE","428814432194878249123pJ","928114432088094735685bD",
          "761014431972956496062hQ","683014431698754553976kP","061114431214513046922nA",
         "057614431101944236880oQ","218214431016364151968nD","962914430829455803398dJ",
         "481914430339197039976xI","822414430245988544348qE","170814430118616102372jJ",
        "539214430061201676635gE", "611814429979693910245uH","470414429497413721010aP",
         "952214429389714735519yM","837014429306406375249uJ","091114429127005121406dY",
         "462114428665982300759bJ","556614428272715617305qE","008614427668708910238gB",
         "486214427429191981003vQ","967614427010637813939sB","456214426861616057017mO",
         "721514426762401388199zN","424314426028753675044cA","974314425871900423839fP",
         "443014425653908099581vP","122014425241172963715cX","245014425131152385402hK",
         "164614425032441587876aK","216114424312033965330dK","059814422605114870557gG",
         "465414422235467276902uN","867914416336081959259nJ","881814414579750916689wW",
         "998314412817679160744nA","579314411859318377664wR","943514406738640054358kD",
         "985514405890109978659fV", "519214404927671581233oA","572814401467502316372hL",
        "899814399222043126464tK","801714398327992307669aW","786214398054650040710jQ")
print(anyDuplicated(jobs))
# make dir and download zip files  

for (j in 2:length(jobs)) { # for every grid cell
  dir_cell<-paste(dir_ukcp,"job",jobs[j],"/",sep="")
  paste(dir_cell)
  dir.create(dir_cell)
  
  for (i in 1:4) {  # for every output zip file
    link<-paste(rootlink,jobs[j],"/output_",i,".zip",sep="")
    print(link)
    file.out<-paste(dir_cell,"output_",i,".zip",sep="")
    print(file.out)
    #temp<-tempfile()
    download.file(link,file.out)
    #unlink(file.out) - would delete downloaded files
  } # end output loop
  
} # end job loop
    

#############################################################
# Check job and write table of job parameters - using jobs list 
#############################################################

# 1. Create 5km matrix containing Weather Generator grid cell ID values
# NB: historic 5km temp files exclude certain cells (where mid point is not land?)
# Input file from: http://ukclimateprojections-ui.metoffice.gov.uk/ui/docs/grids/wg_5km/index.php
#dir_grids<-"~/Documents/Exeter/Data2015/Templates/"
dir_grids<-"C:/Data2015/Templates/"

file.in<-paste(dir_grids,"grid_box_ids_5km.csv",sep="")
print(file.in)
grid.m<-as.matrix(read.csv(file=file.in,header=FALSE))

# Convert -9999 values to NA
sel<-which(grid.m==-9999)
grid.m[sel]<-NA

# Create raster using easting/northings for bottom left of grid cells (centre=+2500)
grid.5km<-raster(nrows=290,ncols=180,xmn=-200000, ymn=-200000, xmx=700000,ymx=1250000,res=c(5000,5000), crs="+init=epsg:27700")
#grid.5km<-raster(nrows=290,ncols=180,xmn=-197500, ymn=-197500, xmx=697500,ymx=1247500,res=c(5000,5000), crs="+init=epsg:27700")
gridid.r<-raster(grid.m,template=grid.5km)
gridsw.r<-crop(x=gridid.r,y=dem) # crop to geographical extent of DEM raster

#### Same for grid cell file with mask  (file possesses one extra column) #### 
file.in<-paste(dir_grids,"grid_box_ids_with_mask.csv",sep="")
print(file.in)
gridmask.m<-as.matrix(read.csv(file=file.in,header=FALSE))

# Create raster using easting/northings for bottom left of grid cells (centre=+2500)
gridmask.5km<-raster(nrows=290,ncols=181,xmn=-205000, ymn=-200000, xmx=700000,ymx=1250000,res=c(5000,5000), crs="+init=epsg:27700")
#grid.5km<-raster(nrows=290,ncols=180,xmn=-197500, ymn=-197500, xmx=697500,ymx=1247500,res=c(5000,5000), crs="+init=epsg:27700")
grididmask.r<-raster(gridmask.m,template=gridmask.5km)
gridcells.r<-crop(x=grididmask.r,y=dem)

# Vector of 5km grid cells for which historic data
gridcells<-getValues(gridcells.r)
landcells<-gridcells[(which(gridcells!=-9999))]

#################################
# 2. Check and save job metadata
#################################
job.table<-array(NA,dim=c(length(jobs),9))
job.table<-data.frame(job.table)
names(job.table)<-c("job","emissions","timeperiod","location","sampmeth","sims",
                    "years","rdmseed","warnings")


for (j in 1:length(jobs)) { # for every job
  # Read metadata file 
  zipfile<-paste(dir_ukcp,"job",jobs[j],"/output_4.zip",sep="")
  metadata<-read.table(unzip(zipfile,files="metadata.xml")) 
  
  emissions<-substr(metadata[3,1],21,23)
  timeperiod<-substr(metadata[4,1],14,22)
  location<-as.character(strsplit(substr(metadata[6,1],11,22),"<")[[1]])[1]
  sampmeth<-substr(metadata[7,1],17,22)
  sims<-substr(metadata[9,1],24,26)
  years<-substr(metadata[10,1],16,17)
  rdmseed<-substr(metadata[11,1],21,21) 
  warn<-""
  
  # Check parameters correct - stop with warning if not
  
  if (emissions!="a1b") {warn<-paste(warn,"emissions not a1b, ",sep="")}
  if (timeperiod!="2040-2069") {warn<-paste(warn,"timeperiod not 2040-69, ",sep="")}
  if (!(location %in% landcells)) {warn<-paste(warn,"location not listed, ",sep="")}
  if (sampmeth!="random") {warn<-paste(warn,"sampmeth not random, ",sep="")}
  if (sims!="500") {warn<-paste(warn,"sims not 500, ",sep="")}
  if (years!="30") {warn<-paste(warn,"years not 30, ",sep="")}
  if (rdmseed!="F") {warn<-paste(warn,"rdmseed not F, ",sep="")}
  
  # if valid then record parameters 
  job.entry<-array(NA,dim=(9))
  job.entry<-c(jobs[j],emissions,timeperiod,location,sampmeth,sims,years,rdmseed, warn)
  job.table[j,]<-job.entry
  
} # end jobs list

# save jobs.table
write.csv(job.table,file=paste(dir_ukcp,"jobstable.csv",sep=""))
#################################
# 3. Print map of historic 5km cells for which WG data available
#################################
sel<-which(gridcells %in% job.table$location)
gridcells<-gridcells(sel)
gridcells.r<-raster(gridcells,template=gridcells.r)


#############################################################
######  INPUT extracted ncdf files to matrix

#############################################################

for (j in 1:length(jobs)) { # for every grid cell
  dir_cell<-paste(dir_ukcp,"job",jobs[j],"/",sep="")
  paste(dir_cell)
  dir.create(dir_cell)
  
  for (i in 1:4) {  # for every output zip file
    link<-paste(rootlink,jobs[j],"/output_",i,".zip",sep="")


# FOR each SIM
# define yearly data file
# input control data from extracted ncdf file
infile<-paste(dir_ukcp,"N32500E137500M2020Sim500/","r_0001_scen_dly.nc",sep="")
print(infile)
ncdf_infile<-nc_open(infile)
scen.data<-ncvar_get(ncdf_infile) # input as matrix

# Downscale to hrly -> 131,400,000 hrs ~ 1TB per variable per cell?!?
# ~1MB per simulation, 100KB per year?
# Solution: analyse by sim>year>hrly steps  
1.Each day downscale to hourly & 100m - SAVE OUTPUTS temp etc
2. FOr each year calculate3 Phenology &  risk factors ->9,000 element results matrix (per variable per cell)
3. Links between years within sim(carryover between years)
4. Breaks between simulations

5. Results matrixes - stat analysis (by sim/year etc)
6. Summary stats for each cell -> raster maps etc


#############################################################
# input scenario data from extracted ncdf file
infile<-paste(dir_ukcp,"N32500E137500M2020Sim500/","r_0001_scen_dly.nc",sep="")
print(infile)
ncdf_infile<-nc_open(infile)
scenario<-ncvar_get(ncdf_infile) # input as matrix


