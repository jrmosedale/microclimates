# Output: three datasets for calculating coastal effects
# TASK 1. Calc inverse wind raster

inv_wind.r<-calc(final,fun=function(x){return(1/(sqrt(x+1)))})

m.invwind<-1/(sqrt(m.str+1))


# TASK 2. Calc 5kmref-cell difference in coast effect using RASTERS
#get valid 5km gridcells
dir_grids<-"~/Documents/Exeter/Data2015/Templates/"
#dir_grids<-"C:/Data2015/Templates/"
dir_coast<-"~/Documents/Exeter/Data2015/CoastEffect/"

in.file<-paste(dir_grids,"ukhistmask.grd",sep="")
print(in.file)
gridmask<-raster(in.file) #  1=valid cell, NA = sea or not data
gridmask.r<-crop(x=gridmask.r,y=dem) # crop to geographical extent of DEM raster
cell5kmxy<-xyFromCell(gridmask.r,1:ncell(gridmask.r)) # centre xy coordinates for 5km grid cells
radius<-20000 # distance of coastal effect
res<-100 # resolution of raster
sect.width<-10
angle<-180

# CHECK input maps - whether value 1 for all non-coastal etc - NEEDS inversing and checking extent matches DEM !!!!

# For each wind direction and coastal effect map
for (angle in seq(0,360,5))
{
  
  # load coast index maps and extract values for centre of 5km cells
  in.file<-paste(dir_coast,"Coast_Index_",sprintf("%03d",angle,sep=""),"_",(radius/res),"km_",sect.width,"deg.r",sep="")
  print(in.file)
  load(file=in.file) # load final.r
  coast.r<-final.r
  # 100m cell numbers corresponding to 5km centre xy coordinates
  cells<-fourCellsFromXY(coast.r,cell5kmxy)
  vals<-matrix(NA,nrow=dim(cells)[1],ncol=4)
  for(n in 1:4) {
    vals[,n]<-extract(coast.r,c(cells[,n]))
  }
  lref5km<-cbind(cell5kmxy,rowMeans(vals)) # define 5km Lref as mean of 4 adjoining 100m cells
  lref5km.r<-setValues(gridmask.r,lref5km[,3]) # create 5km cell raster
  plot(lref5km.r)
  # convert lref back to 100m cell raster
  lref100.r<-disaggregate(lref5km.r,50)
  plot(lref100.r)
  ldif.r<-lref100.r-coast.r
  plot(ldif.r)
} # end for each wind direction

id surrounding 100m cells
get mean values of surrounding cells

Lref<-resample()
Lref<-mask(Lref,hrtemp) # set cells without historic temp values to NA 

# c. Calc SeaST-5kmref difference 