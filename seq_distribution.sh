#!/bin/bash
# Set the working directory 
cd [Working Directory]

source ./conda.sh

conda activate qiime2-amplicon-2024.5
# Sequences are in Casava 1.8 paired-end demultiplexed fastq format
# There are forward (R1) and reverse (R2) reads for each sample

# First, import sequences to QIIME 2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path All_Sequences \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path seq_dist.input_sequences.qza

# After importing, generate a summary of the sequences so distribution and quality can be assessed
qiime demux summarize \
  --i-data seq_dist.input_sequences.qza \
  --o-visualization seq_dist.raw_seqs_quality.qzv

# The output file can be uploaded to https://view.qiime2.org/ to view the results
