## Analysis of Sequences in R using microeco

After processing sequences with QIIME 2, data analysis was completed with the R package microeco (Liu, C. et al. 2021). 

### Steps:
1. QIIME 2 processed sequences were loaded into R and combined with sample metadata to create a 'microtable' object using functions from the file2meco R package (Liu, C. et al. 2022). Mitochondrial and chloroplast sequences were filtered out.

     *(File: create_microtable.R)*

2. Assess the sequencing error rate using the sequenced positive control sample, ZymoBIOMICS Microbial Community DNA Standard (Zymo Research, 2024. The sequencing error rate was determined to be +/- 1.20%, which is within the manufacturers stated error rate (< 15%).

   *(File: sequencing_error_rate.R)*

3. Assess the negative control for potential false positive hits. Note that the negative control sample is not a 'true' negative control, as it was not subjected to DNA extraction. The negative control was supplied by the sequencing laboratory. Several taxa were detected in the negative control sample. Further evaluation of the potential contaminants is required. 

   *(File: negative_control_sample.R)*

4. The decontam R package was used to identify contaminating sequences through both sequence frequency and sequence prevalance approaches as outlined by Davis et al. 2018. Contaminating sequences were removed using the prevalence-based method provided in the decontam R package because a negative-control was used in the study. Additionally, several of the samples had relatively low biomass and prevalence, making the frequency-based method less reliable. In total, 3 taxa were identified as contaminants and removed.

    *(File: decontam_freq.R, decontam_prev.R)*

5. Abundance of cyanobacteria compared to all bacteria detected was assessed across samples.

    *(File: cyano_abund.R)*

   ![cyano_abund_zoom](https://github.com/user-attachments/assets/44efb70f-5878-4663-b520-b2cc02739e8b)
