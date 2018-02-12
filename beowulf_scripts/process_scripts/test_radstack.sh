#!/bin/bash
# Assumes file of list of tar input files before running
# ls $HOME/Data2015/CMSAF-CAL/tar > $HOME/scripts/inputs/CALtarfiles.txt
Syear= 1983
Smonth=1
Sday=1
Eyear=1983
Emonth=3
Emday=30
WKDIR=$HOME"/rscripts" # location of R script
echo "WKDIR= "$WKDIR
cd $WKDIR
. /etc/profile.d/modules.sh
module load shared gdal
module add shared hdf5/1.8.13 udunits2/gcc/2.2.20 netcdf/gcc/64/4.4 R/gcc/3.2.3
which Rscript
Rscript create_radstack_carson.R $Sday $Smonth $Syear $Eday $Emonth $Eyear
