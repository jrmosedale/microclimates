library(raster)
# read in DEM and ensure divisible by 1? km
dem<-raster("~/Documents/Exeter/Data2015/DEM100/demoriginal.tif")
e.dem<-extent(c(110000,420000,0,180000))
dem<-crop(dem,e.dem)
plot(dem,main="DEM-sw")

########################################
# ALTERNATIVE code without loops for any extent/number of blocks
########################################
# Calculate 10km blocks with extra 1km buffer means
mn1km.r<-aggregate(dem,fact=100,fun=mean,na.rm=TRUE)

# Get cell values,x,y,row and col values
m<-getValues(mn1km.r)
x<-xFromCell(mn1km.r,1:ncell(mn1km.r))
y<-yFromCell(mn1km.r,1:ncell(mn1km.r))
block.x<-colFromCell(mn1km.r,1:ncell(mn1km.r))
block.y<-rowFromCell(mn1km.r,1:ncell(mn1km.r))

# Define store dataframe
store<-array(0,dim=c(length(m),11))
store<-data.frame(store)
names(store)<-c("block.x","block.y","val","c.xmn","c.xmx","c.ymn","c.ymx","b.xmn","b.xmx","b.ymn","b.ymx")


# Assign vals to dataframe
store[,1]<-block.x
store[,2]<-block.y
store[,3]<-m
store[,4]<-x-rep(5000,length(m))
store[,5]<-x+rep(5000,length(m))
store[,6]<-y-rep(5000,length(m))
store[,7]<-y+rep(5000,length(m))
store[,8]<-store[,4]-rep(1000,length(m))
store[,9]<-store[,5]+rep(1000,length(m))
store[,10]<-store[,6]-rep(1000,length(m))
store[,11]<-store[,7]+rep(1000,length(m))

# Remove rows with any NA value (ie for mean)
store<-na.omit(store)

write.csv(store,file="~/Documents/Exeter/Data2015/radiation/blocks.csv",row.names=F)



