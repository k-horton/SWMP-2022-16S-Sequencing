#!/bin/bash   

cd [Working Directory]

# load QIIME 2 environment
conda activate qiime2-amplicon-2024.5

# Create QIIME 2 artifact file (.qza) from the manifest file we created
qiime tools import\
  --type 'SampleData[JoinedSequencesWithQuality]'\
  --input-path ./merge_concat/fastq_manifest_noprimer.merge.concat.txt\
  --output-path ./merge_concat/qiime_output/merge_concat.input_sequences.qza\
  --input-format SingleEndFastqManifestPhred33V2
 
# Visualize read quality
qiime demux summarize\
  --i-data ./merge_concat/qiime_output/merge_concat.input_sequences.qza\
  --o-visualization ./merge_concat/qiime_output/merge_concat.quality.qzv
