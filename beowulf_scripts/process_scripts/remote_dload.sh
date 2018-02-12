#!/bin/bash
# Assumes file of list of tar input files before running
FILELIST=$HOME/scripts/inputs/download_list.txt
echo "FILELIST= "$FILELIST
# Read dates from datefile using SGE TASK ID
LINE=`sed -n $SGE_TASK_ID\p $FILELIST`
echo "Line="$LINE
WKDIR=$HOME"/rscripts" # location of R script
echo "WKDIR= "$WKDIR
cd $WKDIR
. /etc/profile.d/modules.sh
wget $LINE