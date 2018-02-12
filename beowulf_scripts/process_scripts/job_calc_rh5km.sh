#!/bin/bash
# Assumes file of list of tar input files before running
# ls $HOME/Data2015/CMSAF-CAL/tar > $HOME/scripts/inputs/CALtarfiles.txt
DATEFILE=$HOME/scripts/inputs/datelist1983-2013.txt
echo "DATEFILE= "$DATEFILE
IFS=","
# Read dates from datefile using SGE TASK ID
LINE=`sed -n $SGE_TASK_ID\p $DATEFILE`
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
Rscript calculate_rh5km_carson.R $Sday $Smonth $Syear $Eday $Emonth $Eyear
