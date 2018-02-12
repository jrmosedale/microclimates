#!/bin/bash
# Assumes file of list of tar input files before running
# ls $HOME/Data2015/CMSAF-CAL/tar > $HOME/scripts/inputs/CALtarfiles.txt
IFS=","
# Read dates from datefile using SGE TASK ID
LINE="1,1,1997,30,3,1997"
echo "Line="$LINE
read -r Sday Smonth Syear Eday Emonth Eyear <<< "$LINE"
echo "Syear= "$Syear
echo "Eyear= "$Eyear
WKDIR=$HOME"/rscripts" # location of R script
echo "WKDIR= "$WKDIR
cd $WKDIR
. /etc/profile.d/modules.sh
module load shared gdal
module add shared hdf5/1.8.13 udunits2/gcc/2.2.20 netcdf/gcc/64/4.4 R/gcc/3.2.3
which Rscript
Rscript create_calstack_v2_carson.R $Sday $Smonth $Syear $Eday $Emonth $Eyear
