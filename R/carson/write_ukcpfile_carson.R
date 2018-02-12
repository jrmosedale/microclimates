# Code to write file of all UKCP 5km cells to be analysed
# Get ukcp09 grid cells coordinates
in.file<-paste(dir_grids,"ukcpmask.grd",sep="")
print(in.file)
gridmask.r<-raster(in.file) #  1=valid cell, NA = sea or not data
gridmask.r<-crop(gridmask.r,e.dem) 
vals<-values(gridmask.r)
xy<-xyFromCell(gridmask.r,1:ncell(gridmask.r))
sel<-which(vals==1)
landcells<-xy[sel,1:2] # = coordinates for middle of each ukcp09 cell
# Write landcell file

