#!/bin/bash
cd [Working Directory]

# Load the cutadapt module
conda activate cutadapt

# Define variables 
FBARCODE=$1 # 5’ forward barcode sequence
RBARCODE=$2 # 5’ reverse barcode sequence. Do not reverse complement
FASTQR1=$3 # name of sample R1 fastq file (everything before .fastq.gz)
FASTQR2=$4 # name of sample R2 fastq file (everything before .fastq.gz)
RLENGTHCUT=$5 #specifies length to trim reads

# Trim primer sequences from the 5’ end of R1 and R2 fastq files using -g and -G respectively.
# -m specifies minimum read length to keep.
cutadapt\
 -g AATGATACGGCGACCACCGAGATCTACA${FBARCODE}ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNCCTAYGGGDBGCWGCAG\
 -G CAAGCAGAAGACGGCATACGAGAT${RBARCODE}GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNGACTACNVGGGTMTCTAATCC\
 -m 1\
 -l $RLENGTHCUT\
 -o $PWD/quality_trim/$FASTQR1.noprimer.ltrim.fastq\
 -p $PWD/quality_trim/$FASTQR2.noprimer.ltrim.fastq\
 $PWD/extracted_sequences/$FASTQR1.fastq	$PWD/extracted_sequences/$FASTQR2.fastq #input R1 and R2 files