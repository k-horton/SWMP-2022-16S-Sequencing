## QIIME 2 Pipeline Overview
Sequencing data was processed and analyzed using a QIIME 2 pipeline adapted from the previously outlined pipeline (Dacey and Chain, 2021) to trim and combine sequences using a combination of merging and concatenating. Assigned taxonomy and sequences were further analysed using the microeco package in R. 

#### Additional Notes
1. All shell script files were generated on a Windows PC, so before running any scripts in Linux or a Windows Subsystem for Linux (WSL) such as Ubuntu, the command **'dos2unix  [file name]'*** must be run to convert the script to Linux format.

   *replace **[file name]** with the name of the script file

3.  Where files say '[Working Directory]', the location of all files being used should be inserted in place of '[Working Directory]'. For example, 'cd /mnt/c/Users/username/Documents/folder1/folder2'

4. There are several points in the pipeline where manifest files were created because the sequences were not in the typical .qza QIIME 2 format.

### Steps:
1. Sequence distribution across samples was evaluated.

     *(File: seq_distribution.sh)*

2. Primer and barcode sequences were removed using cutadapt.

     *(Files: initial_primer_trim.sh and multi_job_initial_primer_trim.sh)*
   
  - Manifest file with locations of sequences created before next step.
  
       *(Files: manifest_filecreate_initial_trim.sh, manifest_file_initial_trim.sh,  multi_job_manifest_initial_trim.sh)*
       
3. Phred quality score (Q-score) distribution plots were used to determine trimming length.

     *(File: initial_quality.sh)*
   
4. The sequence length was trimmed when the median Q-score falls below 20. All sequences were trimmed to the same length to avoid any bias and allow for all ASVs that are the same to be appropriately grouped.

   *(Files: length_trim.sh, multi_job_length_trim.sh)*
   
- Quality of the length-trimmed sequences was verified by creating a new manifest file with locations of length-trimmed sequences and creating Q-score distribution plots.

     *(Files: manifest_filecreate_length_trim.sh, manifest_file_length_trim.sh,  multi_job_manifest_length_trim.sh)*
    
5. Sequences that could be merged (overlap with R1 and R2 reads) were merged with PANDAseq, and any unmerged reads (no overlap) were subsequently concatenated with PANDAseq.

   *(Files: merge_concat.sh, multi_job_merge_coincat.sh)*

- Create a manifest file with the locations of the merged and concatenated reads for each sample before next step.

     *(Files: manifest_filecreate_merge_concat.sh, manifest_file_merge_concat.sh,  multi_job_manifest_merge_concat.sh)*
 
6. After generating merged and concatenated sequences, they were imported to QIIME as single-end reads. DADA2 was used to perform quality filtering, denoising, and chimera removal.

   *(File: dada2_filtering.sh )*
   
7. Taxonomy was assigned to ASVs using the QIIME 2 plugin feature-classifier based on the Silva 138 feature classifer.

    *(File: assign_taxonomy.sh)*

##### Citations
Dacey, D. P., and F. J. J. Chain. 2021. Concatenation of paired-end reads improves taxonomic classification of amplicons for profiling microbial communities. BMC Bioinformatics 22.
