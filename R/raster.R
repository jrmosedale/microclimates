

r<-raster(xmn=0, xmx=200, ymn=0, ymx=100, ncol=20, nrow=10)
r[] <- 1:ncell(r)
e <- extent(10, 220, 10, 100)
r <- extend(r, e)

ncol(r)
nrow(r)
res(r)


r

e<-extent(10,50,0,100)
r<-crop(r,e)
r