#!/bin/bash
cd [Working Directory]

#Script variables
FASTQ=$1 #name of sample fastq file 

# Load the pandaseq module
conda activate pandaseq

# Pandaseq allows merging forward and reverse sequences from paired-end and concatenating those that don't merge. Minimum merging overlap (-o) is 20 because that is DADA2 default.

# Merge/concatenate R1 and R2 primer-removed untrimmed sequences
pandaseq\
  -f $PWD/quality_trim/${FASTQ}_L001_R1_001.noprimer.ltrim.fastq\
  -r $PWD/quality_trim/${FASTQ}_L001_R2_001.noprimer.ltrim.fastq\
  -o 20\
  -g $PWD/merge_concat/${FASTQ}.merged.concat.log.txt\
  -F\
  -w $PWD/merge_concat/${FASTQ}.merged.fastq\
  -U $PWD/merge_concat/${FASTQ}.concat.fastq
  
# Combine these merged and concat reads into one file
cat merge_concat/${FASTQ}.merged.fastq merge_concat/${FASTQ}.concat.fastq > merge_concat/${FASTQ}.merged.concat.fastq

# Compress
gzip merge_concat/${FASTQ}.merged.fastq merge_concat/${FASTQ}.concat.fastq merge_concat/${FASTQ}.merged.concat.fastq
		
