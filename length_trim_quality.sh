#!/bin/bash

cd [Working Directory]

conda activate qiime2-amplicon-2024.5

# Create QIIME 2 artifact file (.qza) from the manifest file we created
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path ./fastq_manifest_noprimer.ltrim.txt\
  --output-path ./qiime_output/length_trim.input_sequences.qza\
  --input-format PairedEndFastqManifestPhred33V2
 
# Visualize read quality
qiime demux summarize\
  --i-data ./qiime_output/length_trim.input_sequences.qza\
  --o-visualization ./qiime_output/length_trim.quality.qzv
