#!/bin/bash
# Assumes file of list of tar input files before running
# ls $HOME/Data2015/CMSAF-CAL/tar > $HOME/scripts/inputs/CALtarfiles.txt
DATEFILE=$HOME/scripts/inputs/datelist1983-2013.txt
echo "DATEFILE= "$DATEFILE
# Read dates from datefile using SGE TASK ID
LINE=`sed -n $SGE_TASK_ID\p $DATEFILE`
echo "Line="$LINE
VARarray=$(echo $IN | tr "," "|n")
Sday=$VARarray[1]
Smonth=$VARarray[2]
Syear=$VARarray[3]
Eday=$VARarray[4]
Emonth=$VARarray[5]
Eyear=$VARarray[6]
WKDIR=$HOME"/rscripts" # location of R script
echo "WKDIR= "$WKDIR
cd $WKDIR
. /etc/profile.d/modules.sh
module load shared gdal
module load shared R
which Rscript
Rscript ds_5kmt_hr_carson.R $Sday $Smonth $Syear $Eday $Emonth $Eyear
