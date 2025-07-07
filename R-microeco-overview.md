## Analysis of 16S Sequences in R using microeco

After processing sequences with QIIME 2, data analysis was completed with the R package microeco (Liu, C. et al. 2021). Note that water quality data has already been analyzed (see repository SWMP-2022-Water-Quality), and will be used to understand the context of the microbial community observed.

The following R (R Core Team, 2024) packages were used for this data analysis and producing figures:
- decontam (Lazarevic et al. 2016)
- dplyr (Wickham et al. 2023a)
- EnvStats (Millard 2013)
- file2meco (Liu et al. 2022a)
- ggplot2 (Wickham 2016a)
- ggpubr (Kassambara 2016)
- gridExtra (Auguie 2010)
- microeco (Liu et al. 2021a)
- phyloseq (McMurdie and Holmes 2013a)
- writexl (Ooms 2017)
- readxl (Wickham and Bryan 2015)
- rstatix (Kassambara 2019)
- stringr (Wickham 2009)
- tidyr (Wickham et al. 2024a)
- vegan (Oksanen et al. 2025)

### Steps:
1. QIIME 2 processed sequences were loaded into R and combined with sample metadata to create a 'microtable' object using functions from the file2meco R package. Mitochondrial and chloroplast sequences were filtered out.

     *(File: create_microtable.R)*

2. Assess the sequencing error rate with the microeco package using the sequenced positive control sample, ZymoBIOMICS Microbial Community DNA Standard (Zymo Research, 2024). The sequencing error rate was determined to be +/- 1.20%, which is within the manufacturers stated error rate (< 15%).

   *(File: sequencing_error_rate.R)*

3. The negative control was assessed for potential false positive hits. Note that the negative control sample is not a 'true' negative control, as it was not subjected to DNA extraction. The negative control was supplied by the sequencing laboratory. Several taxa were detected in the negative control sample. Further evaluation of the potential contaminants is required. 

   *(File: negative_control_sample.R)*

4. The decontam R package was used to identify contaminating sequences through both sequence frequency and sequence prevalance approaches as outlined by Davis et al. 2018. Contaminating sequences were removed using the prevalence-based method provided in the decontam R package because a negative-control was used in the study. Additionally, several of the samples had relatively low biomass and prevalence, making the frequency-based method less reliable. The ggplot2 package was used to plot the prevalence of taxa in true samples against prevalence in negative control (Wickham 2016). In total, 3 taxa were identified as contaminants and removed.

    *(File: decontam_freq.R, decontam_prev.R)*

   Now that the dataframe has been cleaned up, use the file "prepare_dataframe.R" to quickly prepare microtable and WQ data for analyses. The abundance data was normalized by dividing the observed number of reads for the taxon of interest by the total number of observed reads in the library (sample). 

6.  The abundance of cyanobacteria detected with microscopy, qPCR, and 16S sequencing were assessed across field days. 

    *(File: cyano_abund.R)*

Figure 1: Relative abundance of cyanobacteria compared to total bacteria over time based on sequenced V3-V4 region of DNA extracted from stormwater management pond (SWMP) water samples.
![cyano_abund_box](https://github.com/user-attachments/assets/8b1f21a3-841c-4005-a857-d32bb04c9640)

Table 1: p-values for paired-sample t-test of cyanobacterial abundance for each field day conducted following ANOVA.
|  Date    |June 23rd  | July 20th |  Aug 3rd |  Aug 23rd | Aug 31st|
| ---      | ---       | ---       | ---      | ---       | ---     | 
|July 20th |0.0534     |-          |-         |-          |-        |
|Aug 3rd   |0.5958     |1.0000     |-         |-          |-        |
|Aug 23rd  |1.0000     |0.4621     |1.0000    |-          |-        |  
|Aug 31st  |0.8811     |8.6e-05    |0.0027    |0.1584     |-        |  
|Sept 27th |0.4616     |3.8e-05    |0.0012    |0.0769     |1.0000 |

6. Abundance of each cyanobacteria morphology observed with microscopy were compared across field days.

    *(File: community_composition.R)*

7. Potential nitrogen fixers were identified based on 16S sequencing data and the presence of heterocysts observed in the sample. Abundance and presence of potential nitrogen fixers was compared with WQ variables (T-test, Pearson correlation) and across dates (ANOVA) to identify trends. 

    *(File: nitrogen_limitation.R)*

9. Potential microcystin (MC) producers were identified based on 16S sequencing data and data from a qPCR targeting the microcystin synthetase E coding region (*mcyE*). Abundance and presence of potential MC-producers was compared with WQ variables (T-test, Pearson correlation) and across dates (ANOVA) to identify trends. Other toxin producers were also identified using the 16S sequencing data.
 
    *(File: mcyE_presence.R)*

10. The abundance of potential toxin producers and potential n-fixers were compared to determine if there is a link between nitrogen limitation and the proliferation of toxic cyanobacteria in stormwater management ponds. 

     *(File: n_lim_mcyE_relationships.R)*

##### Citations
Auguie, B. 2010, June 20. gridExtra: Miscellaneous Functions for “Grid” Graphics.

Davis, N. M., Di. M. Proctor, S. P. Holmes, D. A. Relman, and B. J. Callahan. 2018. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome 6.

Kassambara, A. 2016, July 20. ggpubr: “ggplot2” Based Publication Ready Plots.

Kassambara, A. 2019, May 27. rstatix: Pipe-Friendly Framework for Basic Statistical Tests.

Lazarevic, V., N. Gaïa, M. Girard, and J. Schrenzel. 2016. Decontamination of 16S rRNA gene amplicon sequence datasets based on bacterial load assessment by qPCR. BMC Microbiology 16.

Liu, C., Cui, Y., Li, X., Yao, M. microeco: an R package for data mining
  in microbial community ecology. FEMS Microbiology Ecology, 2021, Volume 97, Issue 2,
  fiaa255

Liu, C., Li, X., Mansoldo, F.R.P., An, J., Kou, Y., Zhang, X., Wang, J., Zeng, J.,
  Vermelho, A.B., Yao, M. Microbial habitat specificity largely affects microbial
  co-occurrence patterns and functional profiles in wetland soils. Geoderma, 2022. 418, 115866.

McMurdie, P. J., and S. Holmes. 2013. Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE 8.

Millard, S. P. 2013. An R Package for Environmental Statistics. Page R Documentation. Springer, Seattle, Washington.

Oksanen, J., G. L. Simpson, F. G. Blanchet, R. Kindt, P. Legendre, P. R. Minchin, R. B. O’Hara, P. Solymos, M. H. H. Stevens, E. Szoecs, H. Wagner, M. Barbour, M. Bedward, B. Bolker, D. Borcard, T. Borman, G. Carvalho, M. Chirico, M. De Caceres, S. Durand, H. B. A. Evangelista, R. FitzJohn, M. Friendly, B. Furneaux, G. Hannigan, M. O. Hill, L. Lahti, C. Martino, D. McGlinn, M.-H. Ouellette, E. Ribeiro Cunha, T. Smith, A. Stier, C. J. F. Ter Braak, and J. Weedon. 2025. vegan: Community Ecology Package.
Ooms, J. 2017, August 30. writexl: Export Data Frames to Excel “xlsx” Format.

R Core Team. 2024. R: A Language and Environment for Statistical Computing. Vienna, Austria.

Wickham, H. 2009, November 9. stringr: Simple, Consistent Wrappers for Common String Operations.

Wickham, H. 2016. ggplot2: Elegant Graphics for Data Analysis Second Edition. Page (R. Gentleman, K. Hornik, and G. Parmigiani, Eds.). 2nd edition. Springer, Houston, TX

Wickham, H., R. François, L. Henry, K. Müller, and D. Vaughan. 2023, November 17. dplyr: A Grammar of Data Manipulation.

Wickham, H., D. Vaughan, and M. Girlich. 2024, January 24. tidyr: Tidy Messy Data.

Zymo Research, ZymoBIOMICS Microbial Community DNA Standard CAT# D6305, 2024. https://zymoresearch.eu/products/zymobiomics-microbial-community-dna-standard. 

