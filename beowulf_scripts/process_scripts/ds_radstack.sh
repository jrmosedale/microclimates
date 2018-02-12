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
echo "Smonth= "$Smonth
echo "Emonth= "$Emonth
WKDIR=$HOME"/rscripts" # location of R script
echo "WKDIR= "$WKDIR
cd $WKDIR
. /etc/profile.d/modules.sh
module load shared gdal
module load shared R
which Rscript
Rscript downscale_radstack_carson.R $Sday $Smonth $Syear $Eday $Emonth $Eyear
