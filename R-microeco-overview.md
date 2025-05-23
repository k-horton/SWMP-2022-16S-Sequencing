## Analysis of Sequences in R using microeco

After processing sequences with QIIME 2, data analysis was completed with the R package microeco (Liu, C. et al. 2021). Note that water quality data has already been analyzed (see repository SWMP-2022-Water-Quality), and will be used to understand the context of the microbial community observed here.

### Steps:
1. QIIME 2 processed sequences were loaded into R and combined with sample metadata to create a 'microtable' object using functions from the file2meco R package (Liu, C. et al. 2022). Mitochondrial and chloroplast sequences were filtered out.

     *(File: create_microtable.R)*

2. Assess the sequencing error rate using the sequenced positive control sample, ZymoBIOMICS Microbial Community DNA Standard (Zymo Research, 2024). The sequencing error rate was determined to be +/- 1.20%, which is within the manufacturers stated error rate (< 15%).

   *(File: sequencing_error_rate.R)*

3. Assess the negative control for potential false positive hits. Note that the negative control sample is not a 'true' negative control, as it was not subjected to DNA extraction. The negative control was supplied by the sequencing laboratory. Several taxa were detected in the negative control sample. Further evaluation of the potential contaminants is required. 

   *(File: negative_control_sample.R)*

4. The decontam R package was used to identify contaminating sequences through both sequence frequency and sequence prevalance approaches as outlined by Davis et al. 2018. Contaminating sequences were removed using the prevalence-based method provided in the decontam R package because a negative-control was used in the study. Additionally, several of the samples had relatively low biomass and prevalence, making the frequency-based method less reliable. In total, 3 taxa were identified as contaminants and removed.

    *(File: decontam_freq.R, decontam_prev.R)*

5. Abundance of cyanobacteria compared to all bacteria detected was assessed across samples.

    *(File: cyano_abund.R)*

![cyano_abund](https://github.com/user-attachments/assets/faaa3add-9ce5-486d-b2c4-bffd1c4f45d6)

6. Assess community composition of all cyanobacteria detected.

    *(File: community_composition.R)*

7. Identify nitrogen fixers and relationship with nitrogen conditions.

    *(File: nitrogen_fixation.R)*

9. Identify potential microcystin producers and relationship with environmental variables.
 
    *(File: mcyE_producers.R)*

10. Compare molecular methods (qPCR and sequencing) to traditional method (microscopy) used in the study.

     *(File: method_comparison.R)*

##### Citations
Davis, N. M., Di. M. Proctor, S. P. Holmes, D. A. Relman, and B. J. Callahan. 2018. Simple statistical identification and removal of         contaminant sequences in marker-gene and metagenomics data. Microbiome 6.

Liu, C., Cui, Y., Li, X., Yao, M. microeco: an R package for data mining
  in microbial community ecology. FEMS Microbiology Ecology, 2021, Volume 97, Issue 2,
  fiaa255

Liu, C., Li, X., Mansoldo, F.R.P., An, J., Kou, Y., Zhang, X., Wang, J., Zeng, J.,
  Vermelho, A.B., Yao, M. Microbial habitat specificity largely affects microbial
  co-occurrence patterns and functional profiles in wetland soils. Geoderma, 2022. 418, 115866.

Zymo Research, ZymoBIOMICS Microbial Community DNA Standard CAT# D6305, 2024. https://zymoresearch.eu/products/zymobiomics-microbial-community-dna-standard. 

