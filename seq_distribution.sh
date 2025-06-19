# Sequences are in Casava 1.8 paired-end demultiplexed fastq format
# There are forward (R1) and reverse (R2) reads for each sample

# First, import sequences to QIIME 2
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path zipped_sequences \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path sequences.qza

# After importing, generate a summary of the sequences so distribution and quality can be assessed
qiime demux summarize \
  --i-data sequences.qza \
  --o-visualization sequences.qzv

# The output file can be uploaded to https://view.qiime2.org/ to view the results
