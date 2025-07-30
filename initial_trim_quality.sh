#!/bin/bash          	# specifies interpreter for the script, which is the bash shell
 
# load QIIME 2 environment
conda activate qiime2-amplicon-2024.5

# Create QIIME 2 artifact file (.qza) from the manifest file we created
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path fastq_manifest_noprimer.txt\
  --output-path initial_trim.input_sequences.qza\
  --input-format PairedEndFastqManifestPhred33V2
 
# Visualize read quality
qiime demux summarize\
  --i-data initial_trim.input_sequences.qza\
  --o-visualization initial_trim.raw_seqs_quality.qzv
