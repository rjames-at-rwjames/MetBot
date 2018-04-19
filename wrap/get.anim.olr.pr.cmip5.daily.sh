#!/bin/bash

#A script to make animations of daily OLR for noaa and cmip5 

basedir=/ouce-home/students/bras2392/CTdata/metbot_multi_dset/plot_ALLdays4anim
dsetdir=/ouce-home/students/bras2392/PyCharm4CT/quicks/dicts4CDO
dsetlist=$dsetdir/dset_info_4CDO.28models.txt
dom=bigtrop

for year in 2001;do
for model in cdr;do
#for model in $(more $dsetlist | gawk '{print $1}');do
dset=$(grep -w $model $dsetlist | gawk '{print $2}' | head -1)
indir=$basedir/$dset/$model

output=$basedir/$dset.$model.olr.pr.daily.anim.$dom.$year.gif

convert -delay 100 -loop 0 $indir/looped_days.$dset.$model.olr.pr.${year}*.png $output

done
done