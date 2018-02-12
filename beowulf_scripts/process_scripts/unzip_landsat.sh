#!/bin/bash

WKDIR=$HOME"/rscripts" # location of R script
echo "WKDIR= "$WKDIR
cd $WKDIR

. /etc/profile.d/modules.sh
module load shared gdal
module load shared R
which Rscript
Rscript untar_landsat_carson.R 
