#!/bin/bash

# Script to get netcdf files of ltmean
#    uses dictionaries in "dicts4CDO" - recently edited to use CDR2 and ERAI


# Find input and output directories
alldir=../../../CTdata
mbdir=$alldir/metbot_multi_dset
dset_dict=./dicts4CDO/dset_info_4CDO.28models.txt
onlynew=True

# Loop datasets - by name because they are original
echo 'Looping datasets and models'

for name in $(more $dset_dict | gawk '{print $1}');do
    echo $name
    dset=$(grep -w $name $dset_dict  | gawk '{print $2}' | head -1)
    echo $dset

    #Loop variables
    echo 'Looping variables'
    for var in olr pr omega gpth u v q T;do

        echo "Running on"
        echo $var
        dict=./dicts4CDO/dset_info_4CDO.28models.$var.txt
        name2=$(grep -w $name $dict | gawk '{print $3}' | head -1)
        dset2=$(grep -w $name $dict | gawk '{print $4}' | head -1)
        ysname=$(grep -w $name $dict | gawk '{print $5}' | head -1)

        indir=$mbdir/$dset2
        outdir=$indir/$name2/

        infile=$indir/$name2.$var.day.mean.$ysname.nc

        if [ "$dset2" == "trmm" ]; then
            year1=1998
            year2=2013
        else
            if [ "$dset2" == "cmip5" ]; then
                year1=1970
                year2=2004
            else
                year1=1979
                year2=2013
            fi
        fi

        outfile=$outdir/$name.$name2.$var.mon.mean.${year1}_${year2}.nc

        if [ "$onlynew" == "True" ] ; then
            echo "Checking if file already exists"

            if [ -f $outfile2 ]; then
                echo "File already exists"
                continue
            else
                echo "No file yet"
            fi
        fi

	    echo "Generating new file"

	    cdo -ymonmean -seldate,${year1}-01-01,${year2}-12-31 $infile $outfile


    done # variable
done # model