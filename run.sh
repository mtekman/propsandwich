#!/bin/bash

# static - bin
conda=~/miniconda/bin/conda
activate=~/miniconda/bin/activate
deactivate=~/miniconda/bin/deactivate

# static - names
whitelist_out=whitelist.txt
whitelist_raw=whitelist.raw
gzipped=true

# static - dirs
input_files=0_source_files
stage_demulti=1_demultiplex_files

mkdir -p $input_files $stage_demulti

source $activate propsandwich

# Each name represents a single plate (8x12 = 96)
# and the barcodes are unique to a plate only.
plate_names=$1

# @Takes FASTQ from a single plate and generates
#  a count matrix
# $1 - Plate_R1
# $2 - Plate_R2
function umi_whitelist_plate {
    R1_fastq=$1  # We only actually need R1 for the barcodes
    R2_fastq=$2

    umi_tools whitelist\
              --bc-pattern='NNNNNNCCCCCC'\
              --stdin=$R1_fastq\
              --method=reads\
              --plot-prefix=whitelist\
              > $whitelist_out

    cat $whitelist_out\
        | awk '{print NR"\t"$1}'\
              > $whitelist_raw
}


# C (remove barcode sequence from read)
jes_c=true
# ADD (add barcode to read header)
jes_add=true
# BPOS (where the barcode lies)
jes_bpos=READ_1
function jesuite_demultiplex {
    R1_fastq=$1
    R2_fastq=$2
    barcodes=$3
    outp_dir=$4
    
    je demultiplex\
       F1=$R1_fastq F2=$R2_fastq\
       SAME_HEADERS=false\
       BARCODE_FILE=$barcodes\
       BPOS=$jes_bpos C=$jes_c ADD=$jes_add\
       MM=1 MMD=1 Q=10\
       QUALITY_FORMAT=Standard\
       XT=0 ZT=0 RCHAR=: GZ=$gzipped\
       OUTPUT_DIR=$outp_dir\
       KEEP_UNASSIGNED_READ=false\
       STATS_ONLY=false\
       METRICS_FILE_NAME=${stage_demulti}.metrics
}



function stage1_demulti_all {
    names=$1

    while read line; do
        plate_name=$line

        echo "Processing: $plate_name"
        
        R1_fastq=$(find ./ -type f -name "*${plate_name}*.fastq*" | grep "R1" )
        R2_fastq=$(find ./ -type f -name "*${plate_name}*.fastq*" | grep "R2" )

        if [ "$R1_fastq" = "" ] || [ "$R2_fastq" = "" ]; then
            echo "Cannot find R1 or R2 for $plate_name"
            exit -1
        fi

        demultiplex_outdir=$stage_demulti/$plate_name/demultiplex

        mkdir -p $demultiplex_outdir &&
            echo "Extracting barcodes..." &&
            umi_whitelist_plate $R1_fastq $R2_fastq $demultiplex_outdir &&
            echo -e "\n Demultiplexing..." &&
            jesuite_demultiplex $R1_fastq $R2_fastq $whitelist_raw $demultiplex_outdir

    done<$names
}


source $deactivate source propsandwich
