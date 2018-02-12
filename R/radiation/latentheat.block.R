# This function calculates crop reference evapotranspiration (using the Penman-Monteith equation).
  # Inputs and outputs are as rasters excpet for dn which may be raster or a single value
  # Details of algorithm here: http://www.fao.org/docrep/x0490e/x0490e00.htm
# Input variables:
   # Temp is the temperature at the site in degrees C. You will need to use the anomoly between the
          #5 km grid and the 100m cell in the previous time-step to estimate this (as the CRE function is
          # used to derive the local temperature)
   # Net radiation  - importantly this is in MJ m-2 hour-1 and may require conversion of units from the
                    # satellite derived estimates (typical value ~0.2)
   # RH - relative humidity expressed as a percentage (typical value 80%)
   # P - Atmposheric rressure - in millibars (typical value ~ 1000)
         # currently not assumed to vary by location within a 5km grid cell, but we could do an altitude correction:
         # http://www.fao.org/docrep/x0490e/x0490e07.htm#atmospheric pressure (p)
   # dn  # a binary variable specifying whether it is day (1)  or   night (0) needed as the equation for Soil Heat Flux
         # changes depending on whether it is night or day - can be a Raster or single value
   # u2 is wind speed at at 2 m height [m s-1] (typical value ~5)
# Output variable:
  #Crop reference evapotranspiration [mm m-2 hr-1]
CRE<-function(Temp,Rn,RH,P,dn,u2){
     e0<-0.6108*exp(17.27*Temp/(Temp+237.3)) # saturated vapour pressure
     ea<-e0*(RH/100) # actual vapour pressure
     delta<-4098*(0.6108*exp(17.27*Temp/(Temp+237.3)))/((Temp+237.3)^2)  # slope vapour pressure curve
     if(class(dn)=="RasterLayer") G<-overlay(Rn,dn,fun = function(x, y) ifelse(y>0, 0.1*Rn, 0.5*Rn)) else {
       if (dn>0) G<-0.1*Rn else G<-0.5*Rn } # Soil heat flux
     gamma<-0.000665*(P/10)  # psychrometric constant
     ET0<-(0.408*delta*(Rn-G)+gamma*37/(Temp+273)*u2*(e0-ea))/(delta+gamma*(1+0.34*u2))
     ET0<-calc(ET0,fun=function(x){ifelse(x>0,x,0)})
     ET0
}

# This function calculates the change in relative humidity between two locations as a result of the temperature change
# Input variables:
  # t1 is the reference temperature in degrees C - i.e. the value for each 5 km grid cell
  # t2 is the is the temperature at the site in degrees C. You will need to use the anomoly between the
       #5 km grid and the 100m cell in the previous time-step to estimate this
  # rh is the reference relative humdity, expressed as a percentage - i.e. the value for each 5 km grid cell
# Output variable:
 # the relative humdity at the site, expressed as a percentage - i.e. the value for each 100m grid cell
rh.change<-function(t1,t2,rh){
    # absolute humidity at t1
    e0.1<-0.6108*exp(17.27*t1/(t1+237.3))
    e<-e0.1*(rh/100)
    a.hum.1<-(2165*e)/(t1+273.16) # grams per metre cubed
    # rel humidity at t2
    s.e<-(a.hum.1*(t2+273.16))/2165
    s.e0<-0.6108*exp(17.27*t2/(t2+237.3))
    rhs<-(s.e/s.e0)*100
    rhs
}

# This function calculates the amount of water that can be expected to condense,
# as either a result of a change in temperature from one place to anotehr or through time
# Input variables:all as RASTERS
  # t1 is the reference temperature in degrees C - i.e. the value for each 5 km grid cell
  # t2 is the is the temperature at the site in degrees C. You will need to use the anomoly between the
       #5 km grid and the 100m cell in the previous time-step to estimate this
  # rh is the reference relative humdity, expressed as a percentage - i.e. the value for each 5 km grid cell
# Output variable:
 # the amount of water condensed [mm m-2 hr-1]. Zero if relative humidity is less than 100% - as RASTER
Water.conden<-function(t1,t2,rh){
     # absolute humidity at t2, rh=100
     e0.100<-0.6108*exp(17.27*t2/(t2+237.3))
     e100<-e0.100*(100/100)
     a.hum.100<-(2165*e100)/(t2+273.16) # grams per metre cubed
     # absolute humidity at rh2
     rh2<-rh.change(t1,t2,rh)
     e0.2<-0.6108*exp(17.27*t2/(t2+237.3))
     e2<-e0.2*(rh2/100)
     a.hum.2<-(2165*e2)/(t2+273.16) # grams per metre cubed
     a.hum.2-a.hum.100
     #a.hum.2<-ifelse(rh2>100,a.hum.2,0)
     a.hum.2<-overlay(a.hum.2,rh2,fun=function(x,y) ifelse(y>100,x,0) )
     a.hum.2
}


# Generate example input data for 5 km grid
Temp.5km<-10
Temp.5km.prev<-11 # temperature in previous time step
Rn.5km<-0.2
RH.5km<-95
P<-1003
dn<-1
u2.5km<-5

# Generate example inputs for each 100 grid cell
anom<-matrix(rnorm(2500,0,2),nrow=50,ncol=50) # temperature anomaly from previous time step
Temp.100m<-Temp.5km+anom  # proxy t.100m values using previous anomaly
Rn.100m<-Rn.5km+0.02*anom # net radiation 100m 
RH.100m<-rh.change(Temp.5km,Temp.100m,RH.5km) # RH at 100m using proxy temps in calc
u2.100m<-matrix(rnorm(2500,0,1),nrow=50,ncol=50)+u2.5km # wind speed

# Calcute difference in evapotranspiration
CRE.5km<-CRE(Temp.5km,Rn.5km,RH.5km,P,dn,u2.5km)
CRE.100m<-CRE(Temp.100m,Rn.100m,RH.100m,P,dn,u2.100m)
evapdif<-CRE.100m-CRE.5km # regress this against temperature in model

# Calculate difference in water condensation
wc.5km<-Water.conden(Temp.5km.prev,Temp.5km,RH.5km)
wc.100m<-Water.conden((Temp.5km.prev+anom),Temp.100m,RH.5km)
condendif<-wc.100m-wc.5km   # regress this against temperature in model













