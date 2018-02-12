#!/bin/bash
WKDIR=$HOME"/rscripts" # location of R scrip
cd $WKDIR
. /etc/profile.d/modules.sh
module load shared gdal
module load shared R
which Rscript
Rscript albedo.map.R 
