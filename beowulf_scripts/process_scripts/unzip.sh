#!/bin/bash
export ROOT=$HOME"/Data2015/"
WKDIR=$HOME"/rscripts" # location of R script
echo "Root = "$ROOT
echo "WKDIR= "$WKDIR
cd $WKDIR

. /etc/profile.d/modules.sh
module load shared R
which Rscript
Rscript unzip_files_carson.R $ROOT # $HOME = env var 

