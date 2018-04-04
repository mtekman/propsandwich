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

# dry run, echo outputs only
dry=""

R1_fastq=$(find -L $input_dir -type f -name "*${plate_name}*.fastq*" | grep "R1" )
R2_fastq=$(find -L $input_dir -type f -name "*${plate_name}*.fastq*" | grep "R2" )

if [ "$R1_fastq" = "" ] || [ "$R2_fastq" = "" ]; then
    echo -e "\nError: Cannot find R1 or R2 for $plate_name"
    exit -1
fi
mkdir -p $working_dir

echo "Processing: $plate_name with:
- R1 = [$R1_fastq] and 
- R2 = [$R2_fastq]"
sleep 1

# Load bin config and file basenames
source single_plate.config


# Quick check to see whether directory has already been processed
if [ -s $working_dir/$counts_matrix ]; then
    echo "Already processed $plate_name, skipping!"
    exit 0
fi


# Takes FASTQ from a single plate and generates a list of barcodes
function umitools_whitelist_plate {
    local R1_fastq=$1  # We only actually need R1 for the barcodes
    local R2_fastq=$2

    local white_1=$working_dir/$whitelist_out
    local white_2=$working_dir/$whitelist_raw

    if [ "$dry" = "" ]; then

            umi_tools whitelist\
                  --bc-pattern=$bc_pattern\
                  --stdin=$R1_fastq\
                  --method=reads\
                  --log2stderr\
                  --plot-prefix=$working_dir/whitelist\
                  2> $working_dir/whitelist_log.txt\
                  > $white_1 &&

            cat $white_1\
                | awk '{print NR"\t"$1}'\
                      > $white_2 &&
            echo $white_2

    else
        echo $white_2
    fi

}


# Place barcode and umi info into header of sequencing read
#
function umitools_extract_plate {
    local R1_fastq=$1
    local R2_fastq=$2
    local white=$3
    local R1_fix=$working_dir/$R1_fixedfastq
    local R2_fix=$working_dir/$R2_fixedfastq

    if [ "$dry" = "" ]; then   
        umi_tools extract\
                  --bc-pattern=$bc_pattern\
                  --stdin=$R1_fastq\
                  --read2-in=$R2_fastq\
                  --stdout=$R1_fix\
                  --read2-out=$R2_fix\
                  --filter-cell-barcode\
                  --whitelist=$white &&
            echo $R2_fix
    fi
    echo $R2_fix
    # R1_fix and R2_fix generated, but R2_fix is useful to us
}


function getSTARIndex {
    star_index=tmp_star/m16_gencode_merged.gtf.gz

    while ! [ -s $star_index ]; do
        echo "Star index not found, generating..."
        sleep 1
        generate_star_index $star_index
    done

    echo $star_index
}

function rnastar_map {
    local input_sequences=$1

    if [ "$dry" = "" ]; then
        STAR --runThreadN $star_threads\
             --genomeDir $(getSTARIndex)\
             --readFilesIn $input_sequences\
             --readFilesCommand zcat\
             --outFilterMultimapNmax 1\
             --outSAMtype BAM SortedByCoordinate &&
            echo TEST_BAM
    else
       echo TEST_BAM 
    fi
}

function featurecounts {
    local star_bam=$1
    local out_counts=$working_dir/$bam_featurecounts

    if [ "$dry" = "" ]; then
        featureCounts \
            -a $geneset \
            -o $out_counts \
            -R BAM $star_bam \
            -T $fc_threads &&
            echo $out_counts
    else
        echo $out_counts
    fi
}

function samtool_sort_index {
    local bam_fc=$1
    local sorted_indexed_bam=$working_dir/$bam_samtool_sortindex

    if [ "$dry" = "" ]; then
        samtools sort $bam_fc -o $sorted_indexed_bam &&
            samtools index $sorted_indexed_bam &&
            echo $sorted_indexed_bam
    else
        echo $sorted_indexed_bam
    fi

}

function umitools_count {
    local sorted_bam=$1
    local cm=$working_dir/$counts_matrix

    if [ "$dry" = "" ]; then
        umi_tools count\
                  --per-gene\
                  --gene-tag=XT\
                  --per-cell\
                  -I $sorted_bam\
                  -S $cm\
                  --wide-format-cell-counts &&
            echo $cm
    else
        echo $cm
    fi
}

function processPlate {
    # All functions specify the desired output file as the first argument
    # - Also, all files are tempfiles generated in a dir with plenty of space
    echo "Plate $plate_name" &&
        echo "Detecting barcodes..." &&
        # make whitelist_raw
        dry="T"
        whitelist_raw=$(umitools_whitelist_plate $R1_fastq $R2_fastq) &&
        echo -e "\nExtracting barcodes..." &&
        # make R1 and R2 fixedfastq
        extracted_reads=$(umitools_extract_plate $R1_fastq $R2_fastq $whitelist_raw) &&
        echo -e "\nMapping..." &&
        # make bam_star
        dry=""
        bam_star=$(rnastar_map $extracted_reads) &&
        echo -e "\nFeature Counts..." &&
        # make bam_featurecounts
        bam_featurecounts=$(featurecounts $bam_star) &&
        echo -e "\nSorting FC output..." &&
        # make bam_sorted
        bam_sorted=$(samtool_sort_index $bam_featurecounts) &&
        echo -e "\nCounting molecules..." &&
        # make counts_matrix
        counts=$(umitools_count $bam_sorted) &&

        echo $counts
}


function generate_star_index {
    local tmp=tmp_star
    mkdir -p $tmp
    # Links
    local fastam16_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M16/GRCm38.primary_assembly.genome.fa.gz
    local annotm16_url=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz

    #Names
    local genome_fasta=$tmp/$(basename $fastam16_url)
    local annot_gtf=$tmp/$(basename $annotm16_url)
    
    echo "1. Generate GTF file of merged transcripts"
    outname=$1
    tmp1=tmp_star/m16_gencode.gtf.gz
    tmp2=tmp_star/star_m16

    wget -c $annotm16_url -O $tmp1 &&
        # Tool doesn't work unless we look only for lines with transcript_id
        zcat $tmp1 | head | grep "^#" > $tmp1.edited &&
        zcat $tmp1 | grep "transcript_id" >> $tmp1.edited &&
        cgat gtf2gtf --method=merge-exons -I $tmp1.edited | cgat gtf2gtf --method=set-transcript-to-gene | gzip > $tmp2 &&

        mv $star_gtf_tmp $annot_gtf &&

        echo "2. Download Fasta and perform Index generation" &&
        wget -c $fastam16_url -O $genome_fasta &&
        gunzip $genome_fasta &&
        genome_fasta=$tmp/$(basename $genome_fasta .fa.gz).fa &&
        mkdir $tmp/star_tmp &&
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $tmp/star_tmp \
--genomeFastaFiles $genome_fasta --sjdbGTFfile $annot_gtf --sjdbOverhang 99

    

   
#    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir genomedir --genomeFastaFiles GRCm38.primary_assembly.genome.fa  ../tmp_star/m16_gencode.gtf.gz --sjdbOverhang 99
    
#   STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $(getSTARIndex)\
--genomeFastaFiles $fasta

}




function generate_star_flattened_gtf {

    
}


# Load env + execute
source $activate propsandwich
#processPlate
generate_star_index tmp_star/m16_gencode_merged.gtf.gz
source $deactivate propsandwich
