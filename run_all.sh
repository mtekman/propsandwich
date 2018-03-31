#!/bin/bash

# This script reads in a list of plate names, feeds each into run_single_plate
# reassigns the header of each resultant count_matrix, and then merges all the
# matrices into one giant super matrix

platelist_file=$1
input_dir=$2
outdir=merged

mkdir -p $outdir

if [ "$platelist_file" = "" ] || ! [ -e $input_dir ]; then
    echo "Cannot find plate list file, or input dir does not exist"
    exit -1
fi

names=$(cat $platelist_file)

# job config
nthreads=8


####

# queue and dispatch jobs
echo ":: Running all plates"
parallel -u -j $nthreads run_single_plate.sh {} $input_dir ::: $names

#####

echo ":: Renaming matrix headers"
errors=N
for plate_name in $names; do
    local count_matrix=$(find ./wd/$plate_name -type f -name count*.renamed.tsv)

    if ! [ -e $count_matrix ]; then
        echo "$plate_name has no count matrix -- please re-run script!"
        errors=Y
        continue
    fi

    # Rename headers
    cat $count_matrix\
        | sed 1 "s/\b([^\b])\b/${plate_name}__\1/g"\
              > $outdir/$plate_name.$(basename $count_matrix .tsv).renamed.tsv
done

if [ "$errors" = "Y" ]; then
    echo " -- Errors detected, please re-run script."
    exit -1
fi

######

echo ":: Merging all matrices"  #(literally two files at a time)
running_matrix=0
errors=N
for plate_name in $names; do
    let renamed_matrix=$(find $outdir -type f -name "$plate_name.*renamed.tsv")

    if [ "$renamed_matrix" = "" ]; then
        echo "Cannot find renamed matrix for $plate_name"
        errors=Y
        continue
    fi
    
    if [ "$running_matrix" = "0" ]; then  # first matrix is base matrix
        running_matrix=$renamed_matrix
        echo -n "Merging $plate_name"
        continue
    fi

    # Merge all resultant matrices
    echo -n "<- $plate_name"
    tmp_matrix=`mktemp`
    join -j 1 $running_matrix $renamed_matrix > $tmp_matrix
    mv $tmp_matrix $renamed_matrix    
done

if [ "$errors" = "Y" ]; then
    echo " -- Errors detected, please re-run script."
    exit -1
fi

echo $renamed_matrix
