# 1. Coastal effect on temperature

# Set dem scale

# Set time/date scale

# Prior actions
  # Temporal downscale of 5km daily to hrly
  
  # Calculation of windshelter maps = wind_2a_sheltermaps.R
  
  # Downscale wind dir and strength to hr/100m res = wind_1_readdata.R, wind_2b_downscale.R

# Calculation of coastal index
# At each hour:

  # Derive 100res raster from 5km temp data (eventually modified for elevation etc - check existing code)
  
  # Derive 100m coast index raster from central (5km cell value)
  # (disaggregate?? or similar function using )

  


