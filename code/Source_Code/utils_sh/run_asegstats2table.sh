#!/bin/bash

# Script for pulling together the data in aseg.stats using FS's 'asegstats2table'
# This version does NOT need to operate on an FS-style "SUBJECTS_DIR", because it
# makes use of the '--inputs' argument of 'asegstats2table'

# Author: M.P. Harms; Jan 26, 2024


# Use FS 7.X, for which 'asegstats2table' is Python 3 compatible
module load freesurfer/7.4.1 python/3.10.9

SCRIPTDIR=$PWD
PROJDIR=/ceph/intradb/archive/CinaB/AABC_STG
pushd $PROJDIR

sessions=$(ls -1d HCA*_MR)

todate=`date +%Y.%m.%d`
OUTFILE=$SCRIPTDIR/asegstats_AABC_$todate.csv

inputStr=""
for sess in $sessions; do
    statsfile="${sess}/T1w/${sess}/stats/aseg.stats"
    if [ -e $statsfile ]; then
	inputStr+="$statsfile "
    fi
done

echo "inputStr=$inputStr"

# Gather the stats
asegstats2table --inputs $inputStr --skip --delimiter=comma --tablefile=$SCRIPTDIR/tmp.csv

popd

# Using --inputs in 'asegstats2table' results in an integar in the first column rather than the sess id
# So, replace the first column with the session id
cut -d"," -f2- tmp.csv > tmp2.csv
echo -e "Session\n${sessions}" > sessionList.txt
paste -d"," sessionList.txt tmp2.csv > $OUTFILE
rm sessionList.txt tmp.csv tmp2.csv

# Some variable name replacements, to be more similar to what we used in Connectome DB for HCP-YA
# (still quite dissimilar though...)
cp -p $OUTFILE $OUTFILE.orig
sed -e "s/Left-/FS_L_/g" -e "s/Right-/FS_R_/g" -i $OUTFILE
sed -e "s/lh/FS_L_/g" -e "s/rh/FS_R_/g" -i $OUTFILE
sed -e "s/3rd-Ventricle/FS_3rdVent_Vol/" -e "s/4th-Ventricle/FS_4thVent_Vol/" -e "s/5th-Ventricle/FS_5thVent_Vol/" -i $OUTFILE
