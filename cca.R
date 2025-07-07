# Constrained ordination (CCA and RDA) and variance inflation factor
library(decontam)
library(dplyr)
library(file2meco)
library(ggplot2)
library(microeco)
library(phyloseq)
library(readxl)
library(stringr)
library(tidyr)
library(vegan)
#Functions goodness and inertcomp can be used to assess the goodness of fit for 
#individual sites or species. 
#Function vif.cca and alias.cca can be used to analyse linear dependencies among 
#constraints and conditions. 

#The two most commonly used constrained ordination techniques are:
#    1. Redundancy Analysis (RDA) and 
#    2. Canonical Correspondence Analysis (CCA). 
#RDA is the constrained form of PCA, and is inappropriate under the unimodal model. 
#CCA is the constrained form of CA, and therefore is preferred for most ecological data sets 
#(since unimodality is common). CCA also is appropriate under a linear model, as 
#long as one is interested in species composition rather than absolute abundance
#### Prepare microtable as before ####
#create microtable 
{# Assign current working directory to 'dir'
  dir<-getwd()
  # Define the path to each of the files needed to create the microtable object. Replace [Path to file] with the actual path
  abund_file_path <- paste0(dir, "/merge_concat.seqs_count_table.qza")
  sample_file_path <- paste0(dir, "/metadata_decontam.xlsx")
  taxonomy_file_path <- paste0(dir, "/merge_concat.silva_taxonomy.qza")
  tree_data <- paste0(dir, "/merge_concat.rooted_tree.qza")
  rep_data <- paste0(dir, "/merge_concat.representative_seqs.qza")
  tcontam <- qiime2meco(abund_file_path, sample_table = sample_file_path, 
                        taxonomy_table = taxonomy_file_path, phylo_tree = tree_data, 
                        rep_fasta = rep_data, auto_tidy = TRUE)
  # drop the positive control sample from the microtable
  table2<-clone(tcontam)
  table2$sample_table<-subset(table2$sample_table, table2$sample_table$Control!="Pos")
  # Remove pollution (mitochondrial and chloroplast sequences) from the microtable
  table2$filter_pollution(taxa = c("mitochondria", "chloroplast"))
  # convert microtable to phyloseq object
  physeq <- meco2phyloseq(table2)
  physeq
}
#and run decontamination 
{sample_data(physeq)$is.neg <- sample_data(physeq)$Control == "Neg"
  contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold = 0.10)
  # Subset only the sequences designated as contaminants
  prev.contam<-subset(contamdf.prev, contamdf.prev$contaminant == TRUE)
  # Make a list of contaminant sequence IDs
  prev.contam.list<-row.names(prev.contam)
  # Check what the matching taxa are for the contaminant sequence IDs
  prev.contam.taxa<-clone(table2)
  prev.contam.taxa$otu_table<-prev.contam.taxa$otu_table[row.names(prev.contam.taxa$otu_table) %in% c(prev.contam.list),]
  prev.contam.taxa$tidy_dataset()
  prev.taxa.name<-trans_abund$new(dataset = prev.contam.taxa, taxrank="Genus")
  # Remove contaminants from dataframe
  ps.prev.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)
  ps.prev.noncontam
}
#convert phyloseq object back to a microtable and remove control samples
{df <- phyloseq2meco(ps.prev.noncontam)
  no_ctrl <- clone(df)
  no_ctrl$sample_table <- subset(no_ctrl$sample_table, Control == "FALSE")
  no_ctrl$tidy_dataset()}

#### read in WQ data ####
cyano_data<-read_excel(paste0(dir,"/wq_overall_cyano_ASV_abundance.xlsx"))

#### Genus-level ####
# Clone microtable
ord<-clone(no_ctrl)
# Subset to cyanobacteria only
ord$tax_table <- subset(ord$tax_table, (Kingdom == "k__Bacteria" &
                                          Class == "c__Cyanobacteriia"))

ord$tax_table$rows<-rownames(ord$tax_table)
# use updated (2025) taxa names determined through literature review to clean up remaining names
name_change<-read_excel(paste0(dir,"/Taxa_name_changes2.xlsx"))
name_change$Genus<-name_change$Orig_Genus
name_change$Family<-name_change$Orig_Family
name_change_sub<-name_change[c(3:4, 10:11)]

ord$tax_table<- merge(ord$tax_table, name_change_sub, by=c("Genus","Family"), all=TRUE) 
rownames(ord$tax_table)<-ord$tax_table$rows
ord$tax_table<-ord$tax_table[c(3:6, 9,10)]
ord$tax_table<-rename(ord$tax_table, "Genus"="New_Genus",
                      "Family" ="New_Family")
table1<-ord$tax_table

## Genus - filtered
# Filter out taxa that are not identified to the genus level
ord$tax_table<- subset(ord$tax_table, (Genus!="Unidentified Phormidiaceae"&
                                         Genus!="Unidentified Nostocaceae"&
                                         Genus!="Unidentified Nodosilineaceae"&
                                         Genus!="Unidentified Microcystaceae"&
                                         Genus!="Unidentified Leptolyngbyaceae"&
                                         Genus!="Unidentified Cyanobiaceae"&
                                         Genus!="Unidentified Cyanobacteria"))
ord$tidy_dataset()
metadata_ord_con<- cyano_data[,c(1:3,25,26,27,32,35:38)]
# Rename columns 
{colnames(metadata_ord_con)[which(names(metadata_ord_con) == "Cl_mgL")] <- "Cl-"
  colnames(metadata_ord_con)[which(names(metadata_ord_con) == "TP_ugL")] <- "TP"
  colnames(metadata_ord_con)[which(names(metadata_ord_con) == "TOSS_gL")] <- "TOSS"
  colnames(metadata_ord_con)[which(names(metadata_ord_con) == "Nitrate_mgL")] <- "NO3"
  colnames(metadata_ord_con)[which(names(metadata_ord_con) == "DO_mgL")] <- "DO"
  colnames(metadata_ord_con)[which(names(metadata_ord_con) == "Temp_C")] <- "Temp."
  colnames(metadata_ord_con)[which(names(metadata_ord_con) == "Cond_uScm")] <- "Cond."
}

metadata<-as.data.frame(metadata_ord_con)
rownames(metadata) <- metadata[,1]

drop<-drop_na(metadata[, 4:11])

cca_genus <- trans_env$new(dataset = ord, add_data = drop)
cca_genus$cal_ordination(method = "CCA", taxa_level = "Genus")

cca_genus$res_ordination_R2
#r.squared adj.r.squared 
#0.2758682     0.1341800  
vif.cca(cca_genus$res_ordination)
#TOSS    `Cl-`       TP      NO3       pH    Cond.    Temp.       DO    
#5.164814 1.371711 3.318317 2.180924 1.691225 2.290605 1.523556 1.907926 

cca_genus$cal_ordination_anova(taxa_level = "Genus")
cca_genus$res_ordination_terms
# Model: cca(formula = use_data ~ TOSS + `Cl-` + TP + NO3 + pH + Cond. + Temp. + DO, data = env_data)
#         Df ChiSquare      F Pr(>F)  
#TOSS      1   0.04863 0.4615  0.503  
#`Cl-`     1   0.02858 0.2713  0.702  
#TP        1   0.05052 0.4795  0.606  
#NO3       1   0.25655 2.4350  0.040 *
#pH        1   0.15372 1.4591  0.095 .
#Cond.     1   0.35150 3.3363  0.027 *
#Temp.     1   0.10555 1.0018  0.187  
#DO        1   0.16894 1.6035  0.039 *
#Residual 29   3.05538   

aov<-anova(cca_genus$res_ordination)
aov
#Permutation test for cca under reduced model
#Permutation: free
#Number of permutations: 999
# Model: cca(formula = use_data ~ TOSS + `Cl-` + TP + NO3 + pH + Cond. + Temp. + DO, data = env_data)
# Df ChiSquare     F Pr(>F)   
# Model     8    1.1640 1.381  0.002 **
# Residual 29    3.0554         

cca_genus$cal_ordination_envfit(taxa_level = "Genus")
cca_genus$res_ordination_envfit
#          CCA1     CCA2     r2 Pr(>r)   
#TOSS  -0.98643  0.16421 0.0302  0.643   
#Cl-    0.99176 -0.12810 0.0071  0.888   
#TP    -0.93677 -0.34995 0.0036  0.936   
#NO3    0.67570  0.73718 0.3459  0.042 * 
#pH     0.95765 -0.28793 0.0302  0.649   
#Cond.  0.99984 -0.01770 0.8584  0.003 **
#Temp.  0.98536  0.17049 0.0671  0.393   
#DO     0.98623 -0.16537 0.2886  0.104 

summary(cca_genus$res_ordination)
#Call:
#  cca(formula = use_data ~ TOSS + `Cl-` + TP + NO3 + pH + Cond. +      Temp. + DO, data = env_data) 

#Partitioning of scaled Chi-square:
#               Inertia Proportion
#Total           4.219     1.0000
#Constrained     1.164     0.2759
#Unconstrained   3.055     0.7241

#Eigenvalues, and their contribution to the scaled Chi-square 

#Importance of components:
#                       CCA1    CCA2    CCA3    CCA4    CCA5     CCA6     CCA7      CCA8    CA1    CA2    CA3
#Eigenvalue            0.5937 0.20870 0.15803 0.11788 0.05808 0.018952 0.008026 0.0005800 0.9995 0.5351 0.4494
#Proportion Explained  0.1407 0.04946 0.03745 0.02794 0.01376 0.004492 0.001902 0.0001375 0.2369 0.1268 0.1065
#Cumulative Proportion 0.1407 0.19018 0.22764 0.25557 0.26934 0.273829 0.275731 0.2758682 0.5128 0.6396 0.7461
#                     CA4     CA5     CA6     CA7     CA8     CA9     CA10    CA11     CA12      CA13
#Eigenvalue            0.27282 0.22114 0.17096 0.15848 0.08898 0.06450 0.034370 0.03329 0.021215 0.0041997
#Proportion Explained  0.06466 0.05241 0.04052 0.03756 0.02109 0.01529 0.008146 0.00789 0.005028 0.0009953
#Cumulative Proportion 0.81075 0.86316 0.90368 0.94124 0.96233 0.97762 0.985762 0.99365 0.998680 0.9996756
#                     CA14      CA15
#Eigenvalue            0.0008609 0.0005079
#Proportion Explained  0.0002040 0.0001204
#Cumulative Proportion 0.9998796 1.0000000

#Accumulated constrained eigenvalues
#                     Importance of components:
#                       CCA1   CCA2   CCA3   CCA4    CCA5    CCA6     CCA7      CCA8
#Eigenvalue            0.5937 0.2087 0.1580 0.1179 0.05808 0.01895 0.008026 0.0005800
#Proportion Explained  0.5101 0.1793 0.1358 0.1013 0.04990 0.01628 0.006895 0.0004983
#Cumulative Proportion 0.5101 0.6894 0.8252 0.9264 0.97633 0.99261 0.999502 1.0000000

cca_genus$trans_ordination(show_taxa = 40, adjust_arrow_length = TRUE, max_perc_env = 1,
                           max_perc_tax = 1, min_perc_env = 0.2, min_perc_tax = 0.2)

cca_genus$res_ordination_trans$df_sites$Date<-
  factor(cca_genus$res_ordination_trans$df_sites$Date,
         levels = c("June 23rd",  "July 20th", "Aug 3rd",
                    "Aug 23rd",   "Aug 31st",  "Sept 27th"),
         labels = c("June 23rd",    "July 20th",   "Aug 3rd",
                    "Aug 23rd",   "Aug 31st","Sept 27th"))

plot_genus3<-cca_genus$plot_ordination(plot_color = "Date", taxa_text_color="steelblue4",
                                       taxa_text_size=3.5,
                                       taxa_arrow_color = "steelblue4",
                                       env_text_color = "black", env_text_size = 4.5,
                                       env_arrow_color = "black") +
  geom_jitter() + theme(axis.title=element_text(size=15), 
                        legend.text=element_text(size=12),
                        legend.title=element_text(size=14))+
  scale_color_viridis_d()+
  annotate(geom="text", x=max(cca_genus$res_ordination_trans$df_arrows$x/1.025), 
           y=max(cca_genus$res_ordination_trans$df_arrows_spe$y/1.025), 
           label =paste("p < 0.05", "\n(ANOVA)"), 
           colour ="black", size=6)

#### Family level ####
cca_family <- trans_env$new(dataset = ord, add_data = drop)
cca_family$cal_ordination(method = "CCA", taxa_level = "Family")
cca_family$res_ordination_R2
# r.squared adj.r.squared 
#0.4420846     0.2578751 
vif.cca(cca_family$res_ordination)
#     TOSS    `Cl-`       TP      NO3       pH    Cond.    Temp.       DO 
# 5.164814 1.371711 3.318317 2.180924 1.691225 2.290605 1.523556 1.907926  

aov<-anova(cca_family$res_ordination)
aov
#Permutation test for cca under reduced model
#Permutation: free
#Number of permutations: 999
#   Model: cca(formula = use_data ~ TSS + `Cl-` + TP + NO3 + pH + Cond. + Temp. + DO, data = env_data)
#           Df  ChiSquare    F      Pr(>F)   
# Model     8    1.0988 2.8724  0.002 **
# Residual 29    1.3866  

cca_family$trans_ordination(show_taxa = 40, adjust_arrow_length = TRUE, max_perc_env = 1,
                            max_perc_tax = 1, min_perc_env = 0.2, min_perc_tax = 0.2)
cca_family$res_ordination_trans$eigval
cca_family$res_ordination$CCA$eig

cca_family$res_ordination_trans$df_sites$Date<-
  factor(cca_family$res_ordination_trans$df_sites$Date,
         levels = c("June 23rd",  "July 20th", "Aug 3rd",
                    "Aug 23rd",   "Aug 31st",  "Sept 27th"),
         labels = c("June 23rd",    "July 20th",   "Aug 3rd",
                    "Aug 23rd",   "Aug 31st","Sept 27th"))

plot_family<-cca_family$plot_ordination(plot_color = "Date", taxa_text_color="steelblue4",
                                        taxa_text_size=3.5,
                                        taxa_arrow_color = "steelblue4",
                                        env_text_color = "black", env_text_size = 4.5,
                                        env_arrow_color = "black") +
  geom_jitter() + theme(axis.title=element_text(size=15), 
                        legend.text=element_text(size=12),
                        legend.title=element_text(size=14))+
  scale_color_viridis_d()+
  annotate(geom="text", x=max(cca_family$res_ordination_trans$df_arrows$x/1.025), 
           y=max(cca_family$res_ordination_trans$df_arrows_spe$y/1.025), 
           label =paste("p < 0.05", "\n(ANOVA)"), 
           colour ="black", size=6)
