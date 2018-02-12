#!/bin/bash
ROOT=$HOME
WKDIR=$ROOT"/rscripts" # location of R script
echo "Root = "$ROOT
echo "WKDIR= "$WKDIR
cd $WKDIR
. /etc/profile.d/modules.sh
module load shared gdal
module load shared R
which Rscript
Rscript correct_filenames_carson.R 


