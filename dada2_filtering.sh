#DADA2 processing for denoising and chimera removal
qiime dada2 denoise-single\
  --i-demultiplexed-seqs $PWD/merge_concat.input_sequences.qza\
  --p-trim-left 0\
  --p-trunc-len 0\
  --o-table $PWD/merge_concat.seqs_count_table.qza\
  --o-representative-sequences $PWD/merge_concat.representative_seqs.qza\
  --o-denoising-stats $PWD/merge_concat.dada2_stats.qza
  
		#Visualize denoising stats to see how many reads lost at each step.
qiime metadata tabulate\
  --m-input-file $PWD/merge_concat.dada2_stats.qza\
  --o-visualization $PWD/merge_concat.dada2_stats.qzv
  
qiime feature-table summarize \
  --i-table merge_concat.seqs_count_table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data merge_concat.representative_seqs.qza \
  --o-visualization rep-seqs.qzv
