#!/bin/bash
# YEAR from SGE TASK ID
WKDIR=$HOME"/rscripts" # location of R script
echo "WKDIR= "$WKDIR
cd $WKDIR
. /etc/profile.d/modules.sh
module load shared gdal
module add shared hdf5/1.8.13 udunits2/gcc/2.2.20 netcdf/gcc/64/4.4 R/gcc/3.2.3
which Rscript
Rscript write_app_rasters.R $YEAR
