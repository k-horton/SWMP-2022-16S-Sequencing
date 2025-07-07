# Read Me
The code outlined in this project was used to analyze DNA sequences obtained from water samples collected from stormwater management ponds in 2022 as part of academic research at Ontario Tech University. Some elements were adapted from previously published code provided by Dacey and Chain, 2021. Additionally, water quality data was also analyzed, see the repository k-horton/SWMP-2022-Water-Quality (https://github.com/k-horton/SWMP-2022-Water-Quality) for more information.

A DNA library was prepared of the V3-V4 region of the 16S gene using the Quick-16S Plus NGS Library Prep Kit (V3-V4, UDI) (Zymo Research 2024), and then sequenced using paired-end Illumina MiSeq sequencing using the MiSeq Reagent Kit v3 (600-cycle). 

Sequencing data was processed and analyzed using a QIIME 2 (Bolyen, E. et al, 2019) pipeline adapted from the previously outlined pipeline (Dacey and Chain, 2021) to trim and combine sequences using a combination of merging and concatenating. Assigned taxonomy and sequences were further analysed using the microeco and phyloseq packages (Liu et al. 2021, McMurdie and Holmes 2013) in R (R Core Team 2024). Please see the individual documents "R-microeco-overview.md" and "QIIME-2-overview.md" for a full list of references and packages used.

**Citations:**
Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumuham M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodriguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz  C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen lB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, UI-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Votgmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37:852–857.

Dacey, D. P., and F. J. J. Chain. 2021. Concatenation of paired-end reads improves taxonomic classification of amplicons for profiling microbial communities. BMC Bioinformatics 22.

Guerrini, C. J., J. R. Botkin, and A. L. McGuire. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37:852–857.

Liu, C., Cui, Y., Li, X., Yao, M. microeco: an R package for data mining in microbial community ecology. FEMS Microbiology Ecology, 2021, Volume 97, Issue 2, fiaa255

McMurdie, P. J., and S. Holmes. 2013. Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE 8.

R Core Team. 2024. R: A Language and Environment for Statistical Computing. Vienna, Austria.

Zymo Research. 2024. ZymoBIOMICS Microbial Community DNA Standard
