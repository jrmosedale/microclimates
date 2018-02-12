#!/bin/bash
export ROOT=$HOME
WKDIR=$ROOT"/rscripts" # location of R script
echo "Root = "$ROOT
echo "WKDIR= "$WKDIR
cd $WKDIR

. /etc/profile.d/modules.sh
module load shared R
which Rscript
Rscript setup_carson.R $ROOT


