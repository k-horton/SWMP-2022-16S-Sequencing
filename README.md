# Read Me
The code outlined in this project was used to analyze DNA sequences obtained from water samples collected from stormwater management ponds in 2022 as part of academic research at Ontario Tech University. Some elements were adapted from previously published code provided by Dacey and Chain, 2021. 

A DNA library was prepared of the V3-V4 region of the 16S gene using the Quick-16S Plus NGS Library Prep Kit (V3-V4, UDI) (Zymo Research), and then sequenced using paired-end Illumina MiSeq sequencing using the MiSeq Reagent Kit v3 (600-cycle). 

Sequencing data was processed and analyzed using a QIIME 2 pipeline adapted from the previously outlined pipeline (Dacey and Chain, 2021) to trim and combine sequences using a combination of merging and concatenating. Assigned taxonomy and sequences were further analysed using the microeco package in R (Liu et al. 2021).

**Citations:**

Dacey, D. P., and F. J. J. Chain. 2021. Concatenation of paired-end reads improves taxonomic classification of amplicons for profiling microbial communities. BMC Bioinformatics 22.

Liu, C., Cui, Y., Li, X., Yao, M. microeco: an R package for data mining in microbial community ecology. FEMS Microbiology Ecology, 2021, Volume 97, Issue 2, fiaa255
