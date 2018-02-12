#!/bin/bash

FILELIST=$HOME/scripts/inputs/landsat_images.txt
echo "FILELIST= "$FILELIST
# Read dates from datefile using SGE TASK ID
LINE=`sed -n $SGE_TASK_ID\p $FILELIST`
echo "Line="$LINE
WKDIR=$HOME"/rscripts" # location of R scrip
cd $WKDIR

. /etc/profile.d/modules.sh
module load shared gdal
module load shared R
which Rscript
Rscript albedo_carson.R $LINE
