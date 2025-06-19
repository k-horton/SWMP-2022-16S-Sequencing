#!/bin/bash
# Set the working directory 
cd [Working Directory]
conda init

# Load the cutadapt module
conda activate cutadapt

# After looking at the distribution, the primer sequences can be trimmed before assessing read quality.
# The forward primer (341f) is a combination of two primers: CCTACGGGDGGCWGCAG and CCTAYGGGGYGCWGCAG
#       This can be combined into one primer for the purposes of trimming the sequences: CCTAYGGGDBGCWGCAG
# The reverse primer (806r) is: GACTACNVGGGTMTCTAATCC

# Aside from primers, there are also Illumina adapters, random hexamers, and barcodes on each amplicon.
# The combined oligos are:
# F: 5ʹ – AATGATACGGCGACCACCGAGATCTACAXXXXXXXXACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNCCTAYGGGDBGCWGCAG – 3ʹ 
# R: 5′ – CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNGACTACNVGGGTMTCTAATCC – 3′ 

# Note that the barcode (XXXXXXXX) is unique to each sample, so it will be specified when we run the multi-job shell script

FBARCODE=$1 # 5’ forward barcode sequence
RBARCODE=$2 # 5’ reverse barcode sequence. Do not reverse complement
FASTQR1=$3 # name of sample R1 fastq file (everything before .fastq.gz)
FASTQR2=$4 # name of sample R2 fastq file (everything before .fastq.gz)


# Trim primer sequences from the 5’ end of R1 and R2 fastq files using -g and -G respectively.
# -m specifies minimum read length to keep.
# The following code will only trim the primer sequences
cutadapt\
 -g AATGATACGGCGACCACCGAGATCTACA${FBARCODE}ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNCCTAYGGGDBGCWGCAG\
 -G CAAGCAGAAGACGGCATACGAGAT${RBARCODE}GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNGACTACNVGGGTMTCTAATCC\
 -m 1\
 -o $PWD/initial_trim/$FASTQR1.noprimer.fastq\
 -p $PWD/initial_trim/$FASTQR2.noprimer.fastq\
 $PWD/extracted_sequences/$FASTQR1.fastq $PWD/extracted_sequences/$FASTQR2.fastq #input R1 and R2 files
