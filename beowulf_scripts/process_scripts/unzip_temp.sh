#!/bin/bash
ROOT=$HOME
WKDIR=$ROOT"/rscripts" # location of R script
echo "Root = "$ROOT
echo "WKDIR= "$WKDIR
export ZIPDIR=$ROOT"/Data2015/Temp5km/zip/" # location of unzipped files
echo "ZIPIR= "$ZIPDIR
export OUTDIR=$ROOT"/Data2015/Temp5km/unzip/"  # location of extracted files to be written
echo "OUTDIR= "$OUTDIR
cd $WKDIR

. /etc/profile.d/modules.sh
module load shared R
which Rscript
Rscript unzip_temp_carson.R $ZIPDIR $OUTDIR # $HOME = env var 

# Command line: qsub unzip_temp.qsub
