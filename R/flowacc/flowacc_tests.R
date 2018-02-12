# Input Flow Accummulation
dir_flowacc<-"~/Documents/Exeter/Data2015/Flow_acc/"
infile<-paste(dir_flowacc,"flowacc_multi.tif",sep="")
flowacc<-raster(infile,res=100)

par(mfrow=c(1,2))
plot(flowacc)
compareRaster(flowacc,dem)

minval<-cellStats(flowacc,min)
maxval<-cellStats(flowacc,max)

plot(flowacc)
hist(flowacc)

# Best transformation - scale from 0 to 2.5  
flow.log<-log10(flowacc/minval)
plot(flow.log,main="log10flow")
hist(flow.log)


# to get normal distrib 
flow.log<-log10(flowacc-minval)
plot(flow.log,main="log10flow")
hist(flow.log)


# Adjust scale
plot(flowacc)
flow.log<-log10(flowacc-minval)
plot(flow.log,main="log10flow")
hist(flow.log)
plot(flow.log*p.flow,main="flow*param")

plot(flowacc)
flow.ln<-log2(flowacc)
hist(flow.ln)
plot(flow.ln,main="log2flow")
plot(flow.ln*p.flow,main="flow*param")

plot(flowacc)
flow.log1p<-log1p(flowacc)
plot(flow.log1p,main="log1pflow")
hist(flow.log1p)
plot(flow.log1p*p.flow,main="flow*param")


# remove min values
minval<-cellStats(flowacc,min)
flowmin<-flowacc-minval  
plot(flowmin,main="flow-min")
plot(flowmin*p.flow,main="flowmin*minval")  
flowmin.log<-log1p(flowmin)
plot(flowmin.log,main="log")
plot(flowmin.log*p.flow,main="flow.log*p.flow")
hist(flowmin.log)

plot(flowlog)
plot(flowlog*p.flow)

#flowacc<-crop(flowacc,dem)

p.flow<--0.00115536825699

# Write dem
filename<-dir_dem
writeRaster(dem, "~/Documents/Exeter/Data2015/DEM100/dem_used", format="ENVI")