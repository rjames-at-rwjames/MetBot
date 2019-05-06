#!/bin/bash

# Script to select sampled dates from netcdf files and create new netcdf files
#    based on textfile input of dates calculated in python
#    uses dictionaries in "dicts4CDO" - recently edited to use CDR2 and ERAI

rm foo*.nc
rm foo*.lst

# Find input and output directories
alldir=../../../CTdata
mbdir=$alldir/metbot_multi_dset
dset_dict=./dicts4CDO/dset_info_4CDO.28models.txt
onlynew=True

threshtest=True
sample=blon
lag=False

# Loop thresh
echo "Looping thresholds"
for thname in lower actual upper;do
    echo "Running on thresh"
    echo $thname

    # Loop continental and madagascan samples
    echo 'Looping continental and magascan'
    for wcb in cont mada;do
        echo "Running on"
        echo $wcb

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
                smpdir=$indir/$name2/samples/
                txtdir=$smpdir/txtfiles
                outdir=$smpdir/ncfiles

                mkdir -p $outdir

                infile=$indir/$name2.$var.day.mean.$ysname.nc

                outfile2=$outdir/$name.$name2.$var.sampled_days.$sample.$wcb.$thname.day_0.nc

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

                txtfile=$txtdir/datelist.${dset}_${name}.sample_${sample}.$wcb.$thname.day_0.txt

                if [ "$var" == "pr"  ] ; then

                    # Loop sampled days
                    for date in $(more $txtfile);do

                        echo $date
                		year=$(echo $date | gawk -F- '{print $1}')
		                echo $year
		                if [ "$year" -gt "1997" ]; then

                            cdo -seldate,$date $infile foo_${date}.nc

                            echo foo_${date}.nc > foo.lst
                            cat foo.lst >> foo2.lst

                        fi

                        rm foo.lst

                    done

                else

                    # Loop sampled days
                    for date in $(more $txtfile);do

                        echo $date
                        cdo -seldate,$date $infile foo_${date}.nc

                        echo foo_${date}.nc > foo.lst
                        cat foo.lst >> foo2.lst

                        rm foo.lst

                    done

                fi

                files=$(more foo2.lst)

                rm foo2.lst

                rm $outfile2 # deletes previous file if it exists

                cdo -b 32 -mergetime ${files[@]} $outfile2

                rm foo_*.nc

            done # variable
        done # model
    done # sample dom
done # thresh
