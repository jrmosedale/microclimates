# create 4sided pyramid as test dem 

vals<-matrix(NA,ncol=50,nrow=50)
maxc<-ncol(vals); maxr<-nrow(vals)

for (n in (maxr/2):1){
print(n)
sel<-which(col(vals)==n | col(vals)==(maxc-n+1) | row(vals)==n | row(vals)==(maxr-n+1) )
vals[sel]<-n
#print(vals)
}

dem.test<-raster(vals)
plot(dem.test)



# Calculate matrix coords that define central 'block'
xmn<-1+(xmin(dem.block)-xmin(dem.buffer))/res(dem.block)[1]
xmx<-dim(dem.buffer)[1]-(xmax(dem.buffer)-xmax(dem.block))/res(dem.block)[1]
ymn<-1+(ymin(dem.block)-ymin(dem.buffer))/res(dem.block)[1]
ymx<-dim(dem.buffer)[2]-(ymax(dem.buffer)-ymax(dem.block))/res(dem.block)[1]

slope.block<-terrain(dem.block, opt='slope', unit='degrees')
aspect.block<-terrain(dem.block, opt='aspect', unit='degrees')
plot(slope.block)
plot(aspect.block)
