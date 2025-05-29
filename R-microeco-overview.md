## Analysis of Sequences in R using microeco

After processing sequences with QIIME 2, data analysis was completed with the R package microeco (Liu, C. et al. 2021). Note that water quality data has already been analyzed (see repository SWMP-2022-Water-Quality), and will be used to understand the context of the microbial community observed here.

### Steps:
1. QIIME 2 processed sequences were loaded into R and combined with sample metadata to create a 'microtable' object using functions from the file2meco R package (Liu, C. et al. 2022). Mitochondrial and chloroplast sequences were filtered out.

     *(File: create_microtable.R)*

2. Assess the sequencing error rate with the microeco package (Liu, C. et al. 2022) using the sequenced positive control sample, ZymoBIOMICS Microbial Community DNA Standard (Zymo Research, 2024). The sequencing error rate was determined to be +/- 1.20%, which is within the manufacturers stated error rate (< 15%).

   *(File: sequencing_error_rate.R)*

3. Assess the negative control for potential false positive hits using the packages microeco, dplyr, and tidyr (Liu et al. 2021, Wickham et al. 2023, 2024). Note that the negative control sample is not a 'true' negative control, as it was not subjected to DNA extraction. The negative control was supplied by the sequencing laboratory. Several taxa were detected in the negative control sample. Further evaluation of the potential contaminants is required. 

   *(File: negative_control_sample.R)*

4. The decontam R package was used to identify contaminating sequences through both sequence frequency and sequence prevalance approaches as outlined by Davis et al. 2018. Contaminating sequences were removed using the prevalence-based method provided in the decontam R package because a negative-control was used in the study. Additionally, several of the samples had relatively low biomass and prevalence, making the frequency-based method less reliable. The ggplot2 package was used to plot the prevalence of taxa in true samples against prevalence in negative control (Wickham 2016). In total, 3 taxa were identified as contaminants and removed.

    *(File: decontam_freq.R, decontam_prev.R)*

5. Abundance of cyanobacteria compared to all bacteria detected was assessed across samples. The file2meco, microeco, ggplot2, dplyr, and phyloseq packages were used (Liu et al. 2022, 2021, Wickham et al. 2016, 2023, McMurdie and Holmes 2013). 

    *(File: cyano_abund.R)*

Figure 1: Relative abundance of cyanobacteria compared to total bacteria over time based on sequenced V3-V4 region of DNA extracted from stormwater management pond (SWMP) water samples.
![cyano_abund_box](https://github.com/user-attachments/assets/e0fc3af7-f71c-468d-868c-5ae9e8d74bdb)

Table 1: p-values for paired-sample t-test of cyanobacterial abundance for each field day conducted following ANOVA.
|  Date    |June 23rd  | July 20th |  Aug 3rd |  Aug 23rd | Aug 31st|
| ---      | ---       | ---       | ---      | ---       | ---     | 
|July 20th |0.0534     |-          |-         |-          |-        |
|Aug 3rd   |0.5958     |1.0000     |-         |-          |-        |
|Aug 23rd  |1.0000     |0.4621     |1.0000    |-          |-        |  
|Aug 31st  |0.8811     |8.6e-05    |0.0027    |0.1584     |-        |  
|Sept 27th |0.4616     |3.8e-05    |0.0012    |0.0769     |1.0000 |

Figure 2: Total cyanobacterial cells observed over time in stormwater management pond (SWMP) water samples with microscopy.
![t_cyano_abund_box](https://github.com/user-attachments/assets/ce02a93b-49a2-40aa-9221-41a2f24fd46f)


6. Assess community composition of all cyanobacteria detected.

    *(File: community_composition.R)*

7. Identify nitrogen fixers and relationship with nitrogen conditions.

    *(File: nitrogen_fixation.R)*

9. Identify potential microcystin producers and relationship with environmental variables.
 
    *(File: mcyE_producers.R)*

10. Compare molecular methods (qPCR and sequencing) to traditional method (microscopy) used in the study.

     *(File: method_comparison.R)*

##### Citations
Davis, N. M., Di. M. Proctor, S. P. Holmes, D. A. Relman, and B. J. Callahan. 2018. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome 6.

Liu, C., Cui, Y., Li, X., Yao, M. microeco: an R package for data mining
  in microbial community ecology. FEMS Microbiology Ecology, 2021, Volume 97, Issue 2,
  fiaa255

Liu, C., Li, X., Mansoldo, F.R.P., An, J., Kou, Y., Zhang, X., Wang, J., Zeng, J.,
  Vermelho, A.B., Yao, M. Microbial habitat specificity largely affects microbial
  co-occurrence patterns and functional profiles in wetland soils. Geoderma, 2022. 418, 115866.

McMurdie, P. J., and S. Holmes. 2013. Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE 8.

Wickham, H. 2016. ggplot2: Elegant Graphics for Data Analysis Second Edition. Page (R. Gentleman, K. Hornik, and G. Parmigiani, Eds.). 2nd edition. Springer, Houston, TX

Wickham, H., R. François, L. Henry, K. Müller, and D. Vaughan. 2023, November 17. dplyr: A Grammar of Data Manipulation.

Wickham, H., D. Vaughan, and M. Girlich. 2024, January 24. tidyr: Tidy Messy Data.

Zymo Research, ZymoBIOMICS Microbial Community DNA Standard CAT# D6305, 2024. https://zymoresearch.eu/products/zymobiomics-microbial-community-dna-standard. 

