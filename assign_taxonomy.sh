qiime feature-classifier classify-sklearn\
  --i-classifier $PWD/silva-138-99-nb-classifier.qza \
  --i-reads $PWD/merge_concat.representative_seqs.qza\
  --p-read-orientation auto\
  --o-classification $PWD/merge_concat.silva_taxonomy.qza
		  
			#Visualise taxonomy
qiime metadata tabulate\
  --m-input-file $PWD/merge_concat.silva_taxonomy.qza\
  --o-visualization $PWD/merge_concat.silva_taxonomy.qzv
  
#Create .biom files and tree.nwk file and export them to be used in R for analyses
	
	#OTU count table- output file named by default is feature-table.biom
qiime tools export\
  --input-path $PWD/merge_concat.seqs_count_table.qza\
  --output-path $PWD/Biom_files
		#Fix OTU table header 
			#First rename then export the .biom file as a txt, edit, then convert back to .biom
			
mv $PWD/Biom_files/feature-table.biom $PWD/Biom_files/merge_concat.seqs_count_table.biom

biom convert\
  --i $PWD/Biom_files/merge_concat.seqs_count_table.biom\
  --o $PWD/Biom_files/merge_concat.seqs_count_table.txt\
  --to-tsv
  
		#Delete first line of text from this table as it is not needed. Then fix the new first line of the table -
			#need to make sure there are tabs between #OTUID and each sample name. May have to perform 'sed' manually
			#instead of copying and pasting from here. On command line, press [ctrl] and 'v' together and then press [tab] -
			#this will result in a tab.
sed -i '1d' $PWD/Biom_files/merge_concat.seqs_count_table.txt
sed -i "1 s/.*/#OTUID	$SAMPLESTABBED/" $PWD/Biom_files/merge_concat.seqs_count_table.txt
		#Convert this edited feature-table.txt back to a .biom file.
biom convert\
  --i $PWD/Biom_files/merge_concat.seqs_count_table.txt\
  --o $PWD/Biom_files/merge_concat.seqs_count_table.biom\
  --to-hdf5\
  --table-type="OTU table"
  
	#Tree- output file named by default as tree.nwk
qiime tools export\
  --input-path $PWD/merge_concat.rooted_tree.qza\
  --output-path $PWD/Biom_files
		#Rename tree.nwk
mv $PWD/Biom_files/tree.nwk $PWD/Biom_files/merge_concat.tree.nwk


	#Taxonomy GG- output file named by default as taxonomy.tsv
qiime tools export\
  --input-path $PWD/merge_concat.silva_taxonomy.qza\
  --output-path $PWD/Biom_files
		#Rename taxonomy.tsv
mv $PWD/Biom_files/taxonomy.tsv $PWD/Biom_files/merge_concat.silva_taxonomy.tsv
		#Fix gg_taxonomy.tsv table header.
sed -i '1 s/.*/#OTUID	taxonomy	confidence/' $PWD/Biom_files/merge_concat.silva_taxonomy.tsv

