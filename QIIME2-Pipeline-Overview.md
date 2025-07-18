## QIIME 2 Pipeline Overview
Sequencing data was processed and analyzed using a QIIME 2 (Boylen, E. et al, 2019) pipeline adapted from the previously outlined pipeline (Dacey and Chain, 2021) to trim and combine sequences using a combination of merging and concatenating. Assigned taxonomy and sequences were further analysed using the microeco package in R (see *R-microeco-overview.md*). 

#### Additional Notes
1. All shell script files were generated on a Windows PC, so before running any scripts in Linux or a Windows Subsystem for Linux (WSL) such as Ubuntu, the command **'dos2unix  [file name]'*** must be run to convert the script to Linux format.

   *replace **[file name]** with the name of the script file

3.  Where files say '[Working Directory]', the location of all files being used should be inserted in place of '[Working Directory]'. For example, 'cd /mnt/c/Users/username/Documents/folder1/folder2'

4. There are several points in the pipeline where manifest files were created because the sequences were not in the typical .qza QIIME 2 format.

### Steps:
1. Sequence distribution across samples was evaluated.

     *(File: seq_distribution.sh)*

2. Primer and barcode sequences were removed using cutadapt (Martin, M. 2010).

     *(Files: initial_primer_trim.sh and multi_job_initial_primer_trim.sh)*
   
      - Manifest file with locations of sequences created before next step.
  
      *(Files: manifest_filecreate_initial_trim.sh, manifest_file_initial_trim.sh,  multi_job_manifest_initial_trim.sh)*
       
3. Phred quality score (Q-score) distribution plots were used to determine trimming length.

     *(File: initial_quality.sh)*
   
4. The sequence length was trimmed when the median Q-score falls below 20. All sequences were trimmed to the same length to avoid any bias and allow for all ASVs that are the same to be appropriately grouped.

   *(Files: length_trim.sh, multi_job_length_trim.sh)*
   
    - Quality of the length-trimmed sequences was verified by creating a new manifest file with locations of length-trimmed sequences and creating Q-score distribution plots.

     *(Files: manifest_filecreate_length_trim.sh, manifest_file_length_trim.sh,  multi_job_manifest_length_trim.sh)*
    
5. Sequences that could be merged (overlap with R1 and R2 reads > 20 bases) were merged with PANDAseq (Masella, A. P. et al, 2012), and unmerged reads (overlap < 20) were subsequently concatenated with PANDAseq where possible.

   *(Files: merge_concat.sh, multi_job_merge_coincat.sh)*

    - Create a manifest file with the locations of the merged and concatenated reads for each sample before next step.

     *(Files: manifest_filecreate_merge_concat.sh, manifest_file_merge_concat.sh,  multi_job_manifest_merge_concat.sh)*
 
6. After generating merged and concatenated sequences, they were imported to QIIME as single-end reads. DADA2 was used to perform quality filtering, denoising, and chimera removal (Callahan, B. J., 2016).

   *(File: dada2_filtering.sh )*
   
7. Taxonomy was assigned to ASVs using the QIIME 2 plugin feature-classifier based on the Silva 138 feature classifer (Pruesse et al. 2007, Quast et al. 2013, Robeson et al. 2020).

    *(File: assign_taxonomy.sh)*

### Citations
Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumuham M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodriguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz  C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen lB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, UI-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Votgmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37:852–857.

Callahan, B. J., P. J. McMurdie, M. J. Rosen, A. W. Han, A. J. A. Johnson, and S. P. Holmes. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods 13:581–583.

Dacey, D. P., and F. J. J. Chain. 2021. Concatenation of paired-end reads improves taxonomic classification of amplicons for profiling microbial communities. BMC Bioinformatics 22.

Guerrini, C. J., J. R. Botkin, and A. L. McGuire. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37:852–857.

Martin, M. 2010. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17.

Masella, A. P., A. K. Bartram, J. M. Truszkowski, D. G. Brown, and J. D. Neufeld. 2012. PANDAseq: Paired-end assembler for illumina sequences. BMC Bioinformatics 13.

Pruesse, E., C. Quast, K. Knittel, B. M. Fuchs, W. Ludwig, J. Peplies, and F. O. Glöckner. 2007. SILVA: A comprehensive online resource for quality checked and aligned ribosomal RNA sequence data compatible with ARB. Nucleic Acids Research 35:7188–7196.

Quast, C., E. Pruesse, P. Yilmaz, J. Gerken, T. Schweer, P. Yarza, J. Peplies, and F. O. Glöckner. 2013. The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research 41.

Robeson, M. S., D. R. O’Rourke, B. D. Kaehler, M. Ziemski, M. R. Dillon, J. T. Foster, and N. A. Bokulich. 2020, October 5. RESCRIPt: Reproducible sequence taxonomy reference database management for the masses.
