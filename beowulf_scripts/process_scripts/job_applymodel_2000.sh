#!/bin/bash
# Assumes file of list of tar input files before running
# ls $HOME/Data2015/CMSAF-CAL/tar > $HOME/scripts/inputs/CALtarfiles.txt
Sday=1
Smonth=1
Syear=2000
Eday=31
Emonth=12
Eyear=2000
# Read cell from SGE TASK ID
CELL=$SGE_TASK_ID
echo "CELL="$CELL
WKDIR=$HOME"/rscripts" # location of R script
echo "WKDIR= "$WKDIR
cd $WKDIR
. /etc/profile.d/modules.sh
module load shared gdal
module add shared hdf5/1.8.13 udunits2/gcc/2.2.20 netcdf/gcc/64/4.4 R/gcc/3.2.3
which Rscript
Rscript apply2_model_carson.R $Sday $Smonth $Syear $Eday $Emonth $Eyear $CELL
