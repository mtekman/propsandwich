#!/bin/bash

function help {
    echo "This script produces a count matrix that is valid for only one plate. To run this for multiple plates, please call this multiple times with different plate_names, or use the run_all.sh wrapper

    `basename $0` <plate_name> <input_dir>

" >&2
    exit -1
}

[ $# -lt 2 ] && help

 
plate_name=$1
input_dir=$2
working_dir=wd/$plate_name

R1_fastq=$(find $input_dir -type f -name "*${plate_name}*.fastq*" | grep "R1" )
R2_fastq=$(find $input_dir -type f -name "*${plate_name}*.fastq*" | grep "R2" )

if [ "$R1_fastq" = "" ] || [ "$R2_fastq" = "" ]; then
    echo "Cannot find R1 or R2 for $plate_name"
    help
fi
mkdir -p $working_dir

echo "Processing: $plate_name with R1 = [$R1_fastq] and R2 = [$R2_fastq]"
sleep 3

# Load bin config and file basenames
source single_plate.config

# Quick check to see whether directory has already been processed
if [ -e $working_dir/$counts_matrix ]; then
    echo "Already processed $plate_name, skipping!"
    exit 0
fi

# Load env
source $activate propsandwich

# Takes FASTQ from a single plate and generates a list of barcodes
function umitools_whitelist_plate {
    local R1_fastq=$1  # We only actually need R1 for the barcodes
    local R2_fastq=$2

    local white_1=$working_dir/$whitelist_out
    local white_2=$working_dir/$whitelist_raw

    umi_tools whitelist\
              --bc-pattern=$bc_pattern\
              --stdin=$R1_fastq\
              --method=reads\
              --plot-prefix=whitelist\
              > $white_1

    cat $white_1
        | awk '{print NR"\t"$1}'\
              > $white_2
}


# Place barcode and umi info into header of sequencing read
#
function umitools_extract_plate {
    local R1_fastq=$1
    local R2_fastq=$2

    local white=$working_dir/$whitelist_raw
    local R1_fix=$working_dir/$R1_fixedfastq
    local R2_fix=$working_dir/$R2_fixedfastq

    umi_tools extract\
              --bc-pattern=$bc_pattern\
              --stdin=$R1_fastq\
              --stdout=$R1_fix\
              --read2-in=$R2_fastq\
              --read2-out=$R2_fix\
              --filter-cell-barcode\
              --whitelist=$white

    # R1_fix and R2_fix generated, but R2_fix is useful to us

}

function rnastar_map {
    local input_sequences=$1

    STAR --runThreadN $star_threads\
         --genomeDir $star_index\
         --readFilesIn $input_sequences\
         --readFilesCommand zcat\
         --outFilterMultimapNmax 1\
         --outSAMtype BAM SortedByCoordinate
}

function featurecounts {
    local star_bam=$1
    featureCounts \
        -a $geneset \
        -o gene_assigned \
        -R BAM $star_bam \
        -T $fc_threads
}

function samtool_sort_index {
    local bam_fc=$1
    local sorted_indexed_bam=$working_dir/$bam_samtool_sortindex
    samtools sort $bam_fc -o $sorted_indexed_bam
    samtools index $sorted_indexed_bam

    echo $sorted_indexed_bam
}

function umitools_count {
    local sorted_bam=$1
    local cm=$working_dir/$counts_matrix
    umi_tools count\
              --per-gene\
              --gene-tag=XT\
              --per-cell\
              -I $sorted_bam\
              -S $cm\
              --wide-format-cell-counts
    echo $cm
}

function processPlate {
    # All functions specify the desired output file as the first argument
    # - Also, all files are tempfiles generated in a dir with plenty of space
    echo "\n\nPlate $plate_name" &&
        echo "Detecting barcodes..." &&
        # make whitelist_raw
        umitools_whitelist_plate $R1_fastq $R2_fastq &&
        echo -e "\nExtracting barcodes..." &&
        # make R1 and R2 fixedfastq
        umitools_extract_plate $R1_fastq $R2_fastq $whitelist_raw &&
        echo -e "\nMapping..." &&
        # make bam_star
        rnastar_map $extracted_reads &&
        echo -e "\nFeature Counts..." &&
        # make bam_featurecounts
        featurecounts $bam_star &&
        echo -e "\nSorting FC output..." &&
        # make bam_sorted
        samtool_sort_index $bam_featurecounts &&
        echo -e "\nCounting molecules..." &&
        # make counts_matrix
        umitools_count $bam_sorted
}

source $deactivate propsandwich




## Depreciated

# # C (remove barcode sequence from read)
# jes_c=true
# # ADD (add barcode to read header)
# jes_add=true
# # BPOS (where the barcode lies)
# jes_bpos=READ_1
# function jesuite_demultiplex {
#     R1_fastq=$1
#     R2_fastq=$2
#     barcodes=$3
#     outp_dir=$4

#     je demultiplex\
#        F1=$R1_fastq F2=$R2_fastq\
#        SAME_HEADERS=false\
#        BARCODE_FILE=$barcodes\
#        BPOS=$jes_bpos C=$jes_c ADD=$jes_add\
#        MM=1 MMD=1 Q=10\
#        QUALITY_FORMAT=Standard\
#        XT=0 ZT=0 RCHAR=: GZ=$gzipped\
#        OUTPUT_DIR=$outp_dir\
#        KEEP_UNASSIGNED_READ=false\
#        STATS_ONLY=false\
#        METRICS_FILE_NAME=${stage_demulti}.metrics
# }
