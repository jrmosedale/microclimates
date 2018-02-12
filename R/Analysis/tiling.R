root<-"~/Documents/Exeter/Data2015/"


filename<-"~/Documents/Exeter/Data2015/shinydata/rasters/gdd10_gs_mean_1983-2013_cornwall.tif"
test.r<-raster(filename)
plot(test.r)
print(crs(test.r))




The addRasterImage function works by projecting the RasterLayer object to EPSG:3857 and encoding each cell to an RGBA color, to produce a PNG image. That image is then embedded in the map widget.

It’s important that the RasterLayer object is tagged with a proper coordinate reference system. Many raster files contain this information, but some do not. Here is how you’d tag a raster layer object “r” which contains WGS84 data:
  
  crs(r) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")



filename<-"~/Documents/Exeter/Data2015/tiles//Users/jonathanmosedale/Documents/Exeter/Data2015/tiles/gdd10_gs_mean_1983-2013/11/994/1353.png"
tile.r<-raster(filename)



# CReate RGB image
#RGB(r,filename=outfile,)col=palette(n),breaks=c(n+1),zlimcol=,colNA= )
test.r
rgb.r<-RGB(test.r)
plot(rgb.r)
plotRGB(rgb.r)

# Create RGB using specific colour palette
library(colorspace)


outfile<-"~/Documents/Exeter/Data2015/shinydata/rasters/test_gdd10_cornwall.tif"
q<-  quantile(getValues(test.r),c(0,0.05,0.1,0.2,0.3,0.4,0.5, 0.6,0.7, 0.8,0.9,0.95,1),names=FALSE,na.rm=TRUE)
mn<-cellStats(test.r,min)
mx<-cellStats(test.r,max)

p<-  colorBin("YlOrRd",domain=c(round((floor(mn/100)*100),round((ceiling(mx[2])/100)*100)),bins=unique(round(ceiling(q)/100)*100),pretty=TRUE))
p<- colorBin("YlOrRd",)


colrgb.r<-RGB(test.r,col=rev(heat.colors(13)),breaks=q)
colrgb.r<-RGB(test.r,col=p)

colrgb2.r<-RGB(test.r,filename=outfile,col=rev(heat.colors(256)),overwrite=TRUE, colNA='transparent' )
plot(colrgb2.r)
plotRGB(colrgb.r)
writeRaster(colrgb2.r,file=outfile,overwrite=TRUE)
