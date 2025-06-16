# Constrained ordination (CCA and RDA) and variance inflation factor
library(vegan)
library(ggplot2)
library(file2meco)
library(decontam)
library(microeco)
library(phyloseq)
library(dplyr)
library(readxl)
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

#### Prepare data ####
# Clone microtable
ord<-clone(no_ctrl)

# Subset to cyanobacteria only
ord$tax_table <- subset(ord$tax_table, (Kingdom == "k__Bacteria" &
                                          Class == "c__Cyanobacteriia"))

# If wanted, can rename taxa to simplify the names
{ord$tax_table$Genus <- ifelse(ord$tax_table$Genus == "g__Cyanothece_PCC-8801", "Cyanothece spp.", 
                              ifelse(ord$tax_table$Genus == "g__Microcystaceae", "Microcystaceae",
                                     ifelse(ord$tax_table$Genus == "g__Microcystis_PCC-7914", "Microcystis spp.",
                                            ifelse(ord$tax_table$Genus == "g__Snowella_0TU37S04",
                                                   "Snowella spp.",
                                                   ifelse(ord$tax_table$Genus == "g__Synechocystis_PCC-6803",
                                                          "Synechocystis spp.",
                                                          ifelse(ord$tax_table$Genus =="g__Aphanizomenon_NIES81",
                                                                 "Aphanizomenon spp.",
                                                                 ifelse(ord$tax_table$Genus =="g__Cuspidothrix_LMECYA_163",
                                                                        "Cuspidothrix spp.",
                                                                        ifelse(ord$tax_table$Genus =="g__Gloeotrichia_SAG_32.84",
                                                                               "Gloeotrichia spp.",
                                                                               ifelse(ord$tax_table$Genus == "g__Richelia_HH01",
                                                                                      "Richelia spp.",
                                                                                      ifelse(ord$tax_table$Genus == "g__Rivularia_PCC-7116",
                                                                                             "Rivularia spp.",
                                                                                             ifelse(ord$tax_table$Genus == "g__Planktothricoides_SR001",
                                                                                                    "Planktothricoides spp.",
                                                                                                    ifelse(ord$tax_table$Genus =="g__Planktothrix_NIVA-CYA_15",
                                                                                                           "Plankothrix spp.",
                                                                                                           ifelse(ord$tax_table$Genus == "g__Tychonema_CCAP_1459-11B",
                                                                                                                  "Tychonema spp.",
                                                                                                                  ifelse(ord$tax_table$Genus == "g__CENA359" | 
                                                                                                                           ord$tax_table$Genus == "g__JSC-12" |
                                                                                                                           ord$tax_table$Genus == "g__LB3-76" |
                                                                                                                           ord$tax_table$Genus == "g__MIZ36" |
                                                                                                                           ord$tax_table$Genus == "g__RD011",
                                                                                                                         "Leptolyngbyaceae",
                                                                                                                         ifelse(ord$tax_table$Genus == "g__Calothrix_KVSF5",
                                                                                                                                "Calothrix spp.",
                                                                                                                                ifelse(ord$tax_table$Genus =="g__Leptolyngbya_ANT.L52.2" | 
                                                                                                                                         ord$tax_table$Genus== "g__Leptolyngbya_PCC-6306" |
                                                                                                                                         ord$tax_table$Genus == "g__Leptolyngbya_SAG_2411",
                                                                                                                                       "Leptolyngbya spp.",
                                                                                                                                       ifelse(ord$tax_table$Genus == "g__Leptolyngbyaceae",
                                                                                                                                              "Leptolyngbyaceae",
                                                                                                                                              ifelse(ord$tax_table$Genus == "g__Phormidesmis_ANT.L52.6",
                                                                                                                                                     "Phormidesmis spp.",
                                                                                                                                                     ifelse(ord$tax_table$Genus =="g__Phormidium_SAG_37.90",
                                                                                                                                                            "Phormidium spp.",
                                                                                                                                                            ifelse(ord$tax_table$Genus == "g__Limnothrix",
                                                                                                                                                                   "Limnothrix spp.",
                                                                                                                                                                   ifelse(ord$tax_table$Genus == "g__Nodosilinea_PCC-7104",
                                                                                                                                                                          "Nodosilinea spp.",
                                                                                                                                                                          ifelse(ord$tax_table$Genus == "g__Nodosilineaceae",
                                                                                                                                                                                 "Nodosilineaceae",
                                                                                                                                                                                 ifelse(ord$tax_table$Genus == "g__Pseudanabaena_PCC-7429",
                                                                                                                                                                                        "Pseudanabaena spp.",
                                                                                                                                                                                        ifelse(ord$tax_table$Genus =="g__SepB-3",
                                                                                                                                                                                               "Other Cyanobacteria",
                                                                                                                                                                                               ifelse(ord$tax_table$Genus ==  "g__Cyanobium_PCC-6307",
                                                                                                                                                                                                      "Cyanobium spp.",
                                                                                                                                                                                                      ifelse(ord$tax_table$Genus =="g__Schizothrix_LEGE_07164",
                                                                                                                                                                                                             "Schizothrix spp.",
                                                                                                                                                                                                             ifelse(ord$tax_table$Genus == "g__unidentified" & ord$tax_table$Family == "f__Nostocaceae",
                                                                                                                                                                                                                    "Nostocaceae",
                                                                                                                                                                                                                    ifelse(ord$tax_table$Genus =="g__unidentified" & ord$tax_table$Family == "f__Phormidiaceae",
                                                                                                                                                                                                                           "Phormidiaceae",
                                                                                                                                                                                                                           ifelse(ord$tax_table$Genus == "g__unidentified" & ord$tax_table$Family == "f__Leptolyngbyaceae",
                                                                                                                                                                                                                                  "Leptolyngbyaceae",
                                                                                                                                                                                                                                  ifelse(ord$tax_table$Genus == "g__unidentified" & ord$tax_table$Family == "f__Cyanobiaceae",
                                                                                                                                                                                                                                         "Cyanobiaceae",
                                                                                                                                                                                                                                         ifelse(ord$tax_table$Genus == "g__" & ord$tax_table$Class == "c__Cyanobacteriia",
                                                                                                                                                                                                                                                "Other Cyanobacteria",
                                                                                                                                                                                                                                                ifelse(ord$tax_table$Genus == "g__unidentified" & ord$tax_table$Class == "c__Cyanobacteriia",
                                                                                                                                                                                                                                                       "Other Cyanobacteria",
                                                                                                                                                                                                                                                      "ERROR"))))))))))))))))))))))))))))))))} 
# Make sure there are no errors:
table(ord$tax_table$Genus)

# Filter out taxa that are not identified to the genus level
ord$tax_table<- subset(ord$tax_table, (Genus!="Other Cyanobacteria"&
                                         Genus!="Leptolyngbyaceae"&
                                         Genus!="Microcystaceae"&
                                         Genus!="Nodosilineaceae"))
ord$tidy_dataset()

# load in file with metadata
dir<-getwd()
metadata_ord<-read_excel(paste0(dir,"/total_cyano_df.xlsx"))

# subset dataframe
# subset just the columns that are needed
metadata_ord_con<- metadata_ord[,c(2,3:4,17,19,20,25,28:31)]

library(reporter)
colnames(metadata_ord_con)[which(names(metadata_ord_con) == "TSS_gL")] <- "TSS"
colnames(metadata_ord_con)[which(names(metadata_ord_con) == "Cl_mgL")] <- "Cl-"
colnames(metadata_ord_con)[which(names(metadata_ord_con) == "TP_ugL")] <- "TP"
colnames(metadata_ord_con)[which(names(metadata_ord_con) == "Nitrate_mgL")] <- "NO3"
colnames(metadata_ord_con)[which(names(metadata_ord_con) == "DO_mgL")] <- "DO"
colnames(metadata_ord_con)[which(names(metadata_ord_con) == "Temp_C")] <- "Temp."
colnames(metadata_ord_con)[which(names(metadata_ord_con) == "Cond_uScm")] <- "Cond."

# fix sample names
library(stringr)
metadata_ord_fix<-metadata_ord_con %>%
  mutate(
    parts = str_split(Sample, "_", n = 4), 
    Sample = purrr::map_chr(parts, ~ {
      paste0(.x[1], "_", .x[2], .x[3], "_", .x[4])
    })
  ) %>%
  select(-parts) 

metadata<-as.data.frame(metadata_ord_fix)
rownames(metadata) <- metadata[,1]
#### CCA - Genus level ####
# TSS has missing values, so will have to drop some observations to include it. First omit TSS, then include
# No TSS
library(tidyr)
no_tss<-drop_na(metadata[, 5:11])
cca1 <- trans_env$new(dataset = ord, add_data = no_tss)
cca1$cal_ordination(method = "CCA", taxa_level = "Genus")
vif.cca(cca1$res_ordination)
# `Cl-`       TP      NO3       pH    Cond.    Temp.       DO 
# 1.091348 1.690619 1.835965 1.509794 1.586257 1.406352 1.685307 

anova(cca1$res_ordination)
#Permutation test for cca under reduced model
#Permutation: free
#Number of permutations: 999
#   Model: cca(formula = use_data ~ `Cl-` + TP + NO3 + pH + Cond. + Temp. + DO, data = env_data)
#               Df      ChiSquare      F    Pr(>F)   
#   Model       7       1.2278       1.3584  0.005 **
#   Residual   34       4.3903   

# WITH TSS
w_tss<-drop_na(metadata[, 4:11])
cca2 <- trans_env$new(dataset = ord, add_data = w_tss)
cca2$cal_ordination(method = "CCA", taxa_level = "Genus")
vif.cca(cca2$res_ordination)
# TSS    `Cl-`       TP      NO3       pH    Cond.    Temp.       DO 
# 4.681571 1.149135 2.367628 2.130413 1.938407 1.744396 1.855217 2.203724 

aov<-anova(cca2$res_ordination)
#Permutation test for cca under reduced model
#Permutation: free
#Number of permutations: 999
#   Model: cca(formula = use_data ~ TSS + `Cl-` + TP + NO3 + pH + Cond. + Temp. + DO, data = env_data)
#           Df  ChiSquare    F      Pr(>F)   
#  Model     8    1.1791    1.4059  0.007 **
#  Residual 29    3.0402      
cca2$res_ordination_trans$df_sites<-factor(cca2$res_ordination_trans$df_sites,
                                           levels =c("June 23rd",  "July 20th","Aug 3rd",
                                                     "Aug 23rd", "Aug 31st",  "Sept 27th"),
                                           labels = c("June 23rd",  "July 20th","Aug 3rd",
                                                      "Aug 23rd", "Aug 31st",  "Sept 27th"))
cca2$trans_ordination(show_taxa = 40, adjust_arrow_length = TRUE, max_perc_env = 1,
                      max_perc_tax = 1, min_perc_env = 0.2, min_perc_tax = 0.2)

cca2$plot_ordination(plot_color = "Date")

cca2$plot_ordination(plot_color = "Seq_MC_pres")

cca2$plot_ordination(plot_color = "mcyE_Pres")

cca2$plot_ordination(plot_color = "Seq_NF_pres")

cca2$plot_ordination(plot_color = "Het_Pres")

cca2$plot_ordination(plot_color = "Storm_Base")


summary(cca2$res_ordination)



#### CCA - Family level ####
cca3 <- trans_env$new(dataset = ord, add_data = w_tss)
cca3$cal_ordination(method = "CCA", taxa_level = "Family")

vif.cca(cca3$res_ordination)
#   TSS    `Cl-`       TP      NO3       pH    Cond.    Temp.       DO 
#   4.681571 1.149135 2.367628 2.130413 1.938407 1.744396 1.855217 2.203724 

cca3$res_ordination_trans$df_sites<-factor(cca3$res_ordination_trans$df_sites,
                                           levels =c("June 23rd",  "July 20th","Aug 3rd",
                                                     "Aug 23rd", "Aug 31st",  "Sept 27th"),
                                           labels = c("June 23rd",  "July 20th","Aug 3rd",
                                                      "Aug 23rd", "Aug 31st",  "Sept 27th"))

cca3$trans_ordination(show_taxa = 40, adjust_arrow_length = TRUE, max_perc_env = 1,
                      max_perc_tax = 1, min_perc_env = 0.1, min_perc_tax = 0.1)

#### Plot CCA ####
cca2$res_ordination_trans$df_sites$Date<-
  factor(cca2$res_ordination_trans$df_sites$Date,
         levels = c("June 23rd",
                    "July 20th",
                    "Aug 3rd",
                    "Aug 23rd",
                    "Aug 31st",
                    "Sept 27th"),
         labels = c("June 23rd",
                    "July 20th",
                    "Aug 3rd",
                    "Aug 23rd",
                    "Aug 31st",
                    "Sept 27th"))


plot1<-cca2$plot_ordination(plot_color = "Date", taxa_text_color="red3", taxa_text_size=3.5,
                            taxa_arrow_color = "red3",
                            env_text_color = "black", env_text_size = 4.5,
                            env_arrow_color = "black") +
  geom_jitter()

plot1

plot1 + theme(axis.title=element_text(size=15), 
              legend.text=element_text(size=12),
              legend.title=element_text(size=14))+
  scale_color_viridis_d()+
  annotate(geom="text", x=max(cca2$res_ordination_trans$df_arrows$x/1.025), 
           y=max(cca2$res_ordination_trans$df_arrows_spe$y/1.025), 
           label =paste("p = ", aov[1,]$`Pr(>F)`, "\n(ANOVA)"), 
           colour ="black", size=6)

library(svglite)
ggsave("CCA.svg",  plot = last_plot(), 
       path = paste0(dir,"/Ordination"),
       width = 12, height = 8, units = "in")