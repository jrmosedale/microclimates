# Calculate RAM requirements for holding a large 3d matrix
# assumes double precision ie 8 bytes per element (integer=?)
# see: http://blogs.sas.com/content/iml/2014/04/28/how-much-ram-do-i-need-to-store-that-matrix.html

# Calculate from dem the number of land cells
allcells<-length(m.dem)
sel<-which(!is.na(m.dem))
landcells<-length(sel)
print (paste("Number of land cells in DEM = ",landcells," equalling ",((landcells/allcells)*100),"% of all cells",sep=""))

# Calculate number of hours and days within 30 year time period
days<-365
hours<-24*days

GB<-(landcells*8)/10^9
print(paste("GB for storing land cell matrix for single moment: ",GB,sep=""))
GB<-(hours*landcells*8)/10^9
print(paste("GB for hour land cell matrix: ", GB, sep=""))
