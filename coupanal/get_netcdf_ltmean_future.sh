#!/bin/bash

# Script to get netcdf files of ltmean for future
#    for historical I did this from daily data
#    for future I am not sure if I will use daily data for all variables
#    so I am going to get these monthly means from the monthly files on /soge-home/
#    uses dictionaries in "dicts4CDO"


# Find input and output directories
sogdir=/soge-home/data/model/cmip5/
alldir=../../../CTdata
mbdir=$alldir/metbot_multi_dset
dset_dict=./dicts4CDO/dset_info_4CDO.26models.txt
onlynew=False

scen=rcp85
tres=mon
cdom=Amon
run=r1i1p1

# Loop datasets - by name because they are original
echo 'Looping datasets and models'

for name in $(more $dset_dict | gawk '{print $1}');do
    echo $name
    dset=$(grep -w $name $dset_dict  | gawk '{print $2}' | head -1)
    echo $dset

    #Loop variables
    echo 'Looping variables'
    for var in olr pr omega gpth u v q T;do

        # Translate variable names
        vname=$(grep $var var_translate.txt | gawk '{print $1}' | head -1)
        units=$(grep $var var_translate.txt | gawk '{print $3}' | head -1)

        echo "Running on"
        echo $var

        outdir=$mbdir/$dset/$name/

        year1=2065
        year2=2099

        outfile=$outdir/$name.$var.mon.mean.${year1}_${year2}.nc

        if [ "$onlynew" == "True" ] ; then
            echo "Checking if file already exists"

            if [ -f $outfile ]; then
                echo "File already exists"
                continue
            else
                echo "No file yet"
            fi
        fi

        indir=$sogdir/$name/$scen/$tres/$cdom/$run/

        for infile in $(ls $indir/atlas_${vname}_${cdom}_${name}_${scen}_${run}_*.nc);do

	    echo "Generating new file"

            cdo $units -ymonmean -seldate,${year1}-01-01,${year2}-12-31 $infile $outfile

        done # infile listing

    done # variable
done # model
