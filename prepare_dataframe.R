### Prepare dataframes
library(decontam)
library(file2meco)
library(microeco)
library(phyloseq)
library(writexl)
library(readxl)
library(tidyr)
library(dplyr)
library(EnvStats)

#### Prepare Microtable ####
# First, load sequence data and remove contaminant samples as determined previously 
# (see files decontam_prev.R, decontam_freq.R)
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
  # 2 samples with 0 abundance are removed from the otu_table ...
  
  # Drop the positive control sample from the microtable
  table2<-clone(tcontam)
  table2$sample_table<-subset(table2$sample_table, table2$sample_table$Control!="Pos")
  
  # Remove pollution (mitochondrial and chloroplast sequences) from the microtable
  table2$filter_pollution(taxa = c("mitochondria", "chloroplast"))
  # Total 366 features are removed from tax_table ...
  
  # convert microtable to phyloseq object
  physeq <- meco2phyloseq(table2)
  physeq
  
  # Specify negative control sample
  sample_data(physeq)$is.neg <- sample_data(physeq)$Control == "Neg"
  contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold = 0.10)
  table(contamdf.prev$contaminant)
  
  # Identify contaminants
  contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
  # Remove contaminants from dataframe
  ps.prev.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)
  
  # convert phyloseq object back to a microtable:
  df <- phyloseq2meco(ps.prev.noncontam)
  # 14 taxa with 0 abundance are removed from the otu_table ...
  
  # remove the negative control sample
  no_ctrl <- clone(df)
  no_ctrl$sample_table <- subset(no_ctrl$sample_table, Control == "FALSE")
  no_ctrl$tidy_dataset()
}

#### Read in WQ data  ####
wq_data<-read_excel(paste0(dir,"/swmp.xlsx"))
# Fix sample column of WQ data
library(stringr)
wq_data$Sample<-str_replace(wq_data$Sample, "_(\\d)", "\\1")

#### Calculate additional variables ####
# Dissolved inorganic nitrogen
wq_data$DIN<-(wq_data$Nitrate_mgL+wq_data$Nitrite_mgL+wq_data$Ammonia_Ammonium_mgL)
# Separate NH4 + NH3
wq_data$pka <- 0.09018 + (2729.92/(wq_data$Temp_C +273.15))
wq_data$f_NH3 <- 1/(10^(wq_data$pka - wq_data$pH)+1)
wq_data$NH3_mgL <- wq_data$Ammonia_Ammonium_mgL*wq_data$f_NH3
wq_data$NH4_mgL<-wq_data$Ammonia_Ammonium_mgL - wq_data$NH3_mgL
# Nitrogen species : TN ration
wq_data$NO3_TN<-wq_data$Nitrate_mgL / wq_data$Total_N_mgL
wq_data$NO2_TN<-wq_data$Nitrite_mgL / wq_data$Total_N_mgL
wq_data$ON_TN<-wq_data$Organic_N_mgL / wq_data$Total_N_mgL
wq_data$NH3NH4_TN<-wq_data$Ammonia_Ammonium_mgL / wq_data$Total_N_mgL
wq_data$NH3_TN <- wq_data$NH3_mgL / wq_data$Total_N_mgL
wq_data$NH4_TN <- wq_data$NH4_mgL / wq_data$Total_N_mgL
# TOSS : TSS ratio
wq_data$TOSS_TSS <- wq_data$TOSS_gL / wq_data$TSS_gL
# TN : TP ratio
wq_data$TN_TP <- wq_data$Total_N_mgL / (wq_data$TP_ugL*0.001) # <- need to convert to mg/L
# change any infinite values generated to NA
is.na(wq_data) <- sapply(wq_data, is.infinite)

#### Check normality of the data #####
# Chlorophyll-a #
res_aov <- aov(Chla_gL ~ Date, data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Chla_gL+0.01) ~ Date,   data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# This is better - use log transformation 

# Unicellular cyanobacteria cells #
res_aov <- aov(Unicellularcells_L ~ Date,    data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Unicellularcells_L+0.1) ~ Date,  data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
wq_data_sub<-wq_data
wq_data_sub$Unicellularcells_L[wq_data_sub$Unicellularcells_L==0] <- NA
wq_data_sub<-wq_data_sub[complete.cases(wq_data_sub),]

res_aov <- aov(Unicellularcells_L ~ Date,  data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Unicellularcells_L+0.01) ~ Date,
               data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

#this looks good! Use log transformation

# Colonial cyanobacteria cells #
res_aov <- aov(Colonialcells_L ~ Date,  data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Colonialcells_L+0.01) ~ Date, data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# This is better - use log transformation 


# Filamentous cyanobacteria cells #
res_aov <- aov(Filamentouscells_L ~ Date, data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Filamentouscells_L+0.01) ~ Date,  data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# This is questionable, try some other transformations:

# root transformation -> x^(1/2)
res_aov <- aov((Filamentouscells_L+0.01)^(1/2) ~ Date, data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# reciprocal transformation -> 1/x
res_aov <- aov(1/(Filamentouscells_L+0.01) ~ Date,  data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# squared transformation -> x^2
res_aov <- aov((Filamentouscells_L+0.01)^(2) ~ Date,
               data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
wq_data_sub<-wq_data
wq_data_sub$Filamentouscells_L[wq_data_sub$Filamentouscells_L==0] <- NA
wq_data_sub<-wq_data_sub[complete.cases(wq_data_sub),]

res_aov <- aov(Filamentouscells_L ~ Date,   data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Filamentouscells_L+0.01) ~ Date, data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

#this looks good! Use log transformation

# Total cyanobacteria cells #
res_aov <- aov(Totalcells_L ~ Date, data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Totalcells_L+0.01) ~ Date, data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# This is better - use log transformation 

# 16S copies / L #
res_aov <- aov(copies_16S_L ~ Date, data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(copies_16S_L+0.01) ~ Date,  data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
wq_data_sub<-wq_data
wq_data_sub$copies_16S_L[wq_data_sub$copies_16S_L==0] <- NA
wq_data_sub<-wq_data_sub[complete.cases(wq_data_sub),]

res_aov <- aov(copies_16S_L ~ Date,data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

res_aov <- aov(log10(copies_16S_L+0.01) ~ Date, data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# this is good, use log transformation

# mcyE copies / L #
res_aov <- aov(copies_mcyE_L ~ Date,   data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(copies_mcyE_L+0.01) ~ Date,    data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
wq_data_sub<-wq_data
wq_data_sub$copies_mcyE_L[wq_data_sub$copies_mcyE_L==0] <- NA
wq_data_sub<-wq_data_sub[complete.cases(wq_data_sub),]

res_aov <- aov(copies_mcyE_L ~ Date, data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

res_aov <- aov(log10(copies_mcyE_L+0.01) ~ Date,   data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# this is better, use log transformation

# mcyE copies : 16S copies #
res_aov <- aov(mcyE_16S_ratio ~ Date,    data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(mcyE_16S_ratio+0.01) ~ Date, data = wq_data)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
wq_data_sub<-wq_data
wq_data_sub$mcyE_16S_ratio[wq_data_sub$mcyE_16S_ratio==0] <- NA
wq_data_sub<-wq_data_sub[complete.cases(wq_data_sub),]

res_aov <- aov(mcyE_16S_ratio ~ Date,  data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

res_aov <- aov(log10(mcyE_16S_ratio+0.01) ~ Date,data = wq_data_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# this is better, use log transformation
# So, to summarize, the following transformations will be applied when conducting
#   parametric statistical tests:
#   - Chlorophyll-a     ->  log10(Chla_gL+0.01) 
#   - Unicellular cells ->  log10(Unicellularcells_L+0.01)
#   - Colonial cells    ->  log10(Colonialcells_L+0.01) 
#   - Filamentous cells ->  log10(Filamentouscells_L+0.01)
#   - Total cells       ->  log10(Totalcells_L+0.01)
#   - 16S copies        ->  log10(copies_16S_L+0.01) 
#   - mcyE copies       ->  log10(copies_mcyE_L+0.01) 
#   - mcyE : 16S copies ->  log10(mcyE_16S_ratio+0.01) 


#### Transform data ####
wq_data$LogAbs_440<-log10(wq_data$Abs_440 + 0.001)
wq_data$LogAbs_750<-log10(wq_data$Abs_750 + 0.001)
wq_data$LogTSS<-log10(wq_data$TSS_gL + 0.001)
wq_data$LogTOSS<-log10(wq_data$TOSS_gL + 0.001)
wq_data$LogCl<-log10(wq_data$Cl_mgL + 0.1)
wq_data$LogTP<-log10(wq_data$TP_ugL + 0.1)
wq_data$LogTDP<-log10(wq_data$TDP_ugL + 0.1)
wq_data$LogTKN<-log10(wq_data$Total_Kjeldahl_N_mgL + 0.1)
wq_data$LogNH3NH4<-log10(wq_data$Ammonia_Ammonium_mgL + 0.01)
wq_data$LogNO3<-log10(wq_data$Nitrate_mgL + 0.001)
wq_data$LogNO2<-log10(wq_data$Nitrite_mgL + 0.001)
wq_data$LogON<-log10(wq_data$Organic_N_mgL + 0.1)
wq_data$LogTN<-log10(wq_data$Total_N_mgL + 0.01)
wq_data$LogDIN<-log10(wq_data$DIN + 0.001)
wq_data$LogChla<-log10(wq_data$Chla_gL + 0.001)
wq_data$LogNH3<-log10(wq_data$NH3_mgL + 0.0001)
wq_data$LogNH4<-log10(wq_data$NH4_mgL + 0.01)
wq_data$LogNO3_TN<-log10(wq_data$NO3_TN + 0.001)
wq_data$LogNO2_TN<-log10(wq_data$NO2_TN + 0.01)
wq_data$LogNH3NH4_TN<-log10(wq_data$NH3NH4_TN + 0.01)
wq_data$LogNH3_TN<-log10(wq_data$NH3_TN + 0.001)
wq_data$LogNH4_TN<-log10(wq_data$NH4_TN + 0.01)
wq_data$LogTN_TP<-log10(wq_data$TN_TP + 0.1)

wq_data$logChla<-log10(wq_data$Chla_gL+0.01)

wq_data$log16S<-log10(wq_data$copies_16S_L+0.01)
wq_data$logmcyE<-log10(wq_data$copies_mcyE_L+0.01)
wq_data$logU_cyano<-log10(wq_data$Unicellularcells_L+0.01)
wq_data$logC_cyano<-log10(wq_data$Colonialcells_L+0.01)
wq_data$logF_cyano<-log10(wq_data$Filamentouscells_L+0.01)
wq_data$logT_cyano<-log10(wq_data$Totalcells_L+0.01)
wq_data$logmcyE_16S<-log10(wq_data$mcyE_16S_ratio+0.01)
wq_data$logHet<-log10(wq_data$Total_Het_L+0.01)

#### Assign labels to categorical variables ####
wq_data$Storm_Base<-factor(wq_data$Storm_Base ,
                               levels = c("Storm","Base"),
                               labels = c("Stormflow","Low Flow"))

wq_data$Site<-factor(wq_data$Site ,
                         levels = c("in","pond","out"),
                         labels = c("Inflow","Pond","Outflow"))

wq_data$Dredge_Cat<-factor(wq_data$Dredge_Cat,
                               levels = c("Recent", "Needs Dredging"),
                               labels = c("Recently Dredged", "Needs Dredging"))

wq_data$Pre_2003<-factor(wq_data$Pre_2003,
                             levels = c("Yes", "No"),
                             labels = c("Built Pre-2003", "Built Post-2003"))


#### Calculate abundance data ####

# There are a few things we need to extract here:
#   1. Overall bacteria library size and cyanobacteria library size for each sample
#   2. Relative abundance of cyanobacteria as a whole out of entire library for each sample.
#   3. Relative abundance of individual cyano taxa out of entire cyanobacteria library for each sample

# 1. Overall library size for each sample
# Bacteria library size
all_bact<-clone(no_ctrl)
all_bact$tax_table<-subset(all_bact$tax_table, Kingdom=="k__Bacteria")
all_bact$tidy_dataset()

bact_otu_table<-all_bact$otu_table
bact_library_size<-colSums(bact_otu_table[c(1:69)])
print(bact_library_size)

bact_library_size<-as.data.frame(bact_library_size)
bact_library_size<- data.frame(Sample = row.names(bact_library_size), 
                               bact_library_size, row.names=NULL)

# Cyanobacteria library size
all_cyano<-clone(no_ctrl)
all_cyano$tax_table<-subset(all_cyano$tax_table, Class=="c__Cyanobacteriia")
all_cyano$tidy_dataset()

cyano_otu_table<-all_cyano$otu_table
cyano_library_size<-colSums(cyano_otu_table[c(1:46)])
print(cyano_library_size)

cyano_library_size<-as.data.frame(cyano_library_size)
cyano_library_size<- data.frame(Sample = row.names(cyano_library_size), 
                                cyano_library_size, row.names=NULL)

library_size<-left_join(bact_library_size, cyano_library_size, by="Sample")

dir<-getwd()
write_xlsx(library_size, path = paste0(dir,"/library_size.xlsx"))

#   2. Relative abundance of cyanobacteria as a whole out of entire library for each sample.

# reads
head(bact_library_size)
# associated taxonomy
bact_tax_table<-all_bact$tax_table

# Ensure the row names of otu_table match the row names of tax_table and re-order if necessary
if (!all(rownames(bact_otu_table) == rownames(bact_tax_table))) {
  warning("Rownames of otu_table and tax_table do not match exactly.
          Reordering otu_table based on tax_table rownames.")
  bact_otu_table <- bact_otu_table[rownames(bact_tax_table), ]
}

# merge taxa and otu tables
merged_table <- merge(bact_tax_table, bact_otu_table, by=0, all=TRUE) 
write_xlsx(merged_table, path = paste0(dir,"/all_Bacteria_ASV_abundance.xlsx"))

# Define function to clean up the taxonomic prefixes
clean_taxa_name <- function(x) {
  ifelse(grepl("__", x), sub("^[^_]*__", "", x), x)
}

# drop species and rownames columns
table_rename<-merged_table[c(2:7,9:77)]

table_rename2 <- table_rename %>%
  mutate( 
    # Apply cleaning to all relevant columns first if they contain prefixes like "g__"
    across(c(Genus, Family, Order, Class, Phylum, Kingdom), clean_taxa_name),
    # redefine the Genus column based on conditions
    Genus = case_when(
      # Class = Cyanobacteriia
      (!is.na(Class) & (Class == "Cyanobacteriia" )) ~ paste0("Cyanobacteria"),
      
      # Class =/= Cyanobacteria
      (!is.na(Class) & (Class != "Cyanobacteriia" )) ~ paste0("Other Bacteria"),
      TRUE ~ "Unassigned_Taxa")
  )

cyano_vs_bact_sample_reads<-table_rename2[c(6:75)] %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), sum))


pivot<-pivot_longer(cyano_vs_bact_sample_reads, cols=c(2:70),
                    names_to = "Sample", values_to = "Reads")

pivot <- pivot %>% select(Sample, everything()) %>%pivot_wider(names_from="Genus", 
                                                               values_from="Reads")
# merge table w/ the total library size for each sample
overall_cyano_abund <- merge(pivot, library_size, by=1, all=TRUE) 

#   3. Relative abundance of individual cyano taxa out of entire library for each sample

cyano_otu_table<-all_cyano$otu_table
cyano_tax_table<-all_cyano$tax_table

# merge taxa and otu tables
merged_table <- merge(cyano_tax_table, cyano_otu_table, by=0, all=TRUE) 

# drop species, class, phylum, kingdom, and rownames columns
table_rename<-merged_table[c(5:7,9:54)]
print(unique(table_rename$Genus))
#[1] "g__Leptolyngbya_ANT.L52.2"   "g__"                         "g__Cyanobium_PCC-6307"      
#[4] "g__Phormidesmis_ANT.L52.6"   "g__RD011"                    "g__Synechocystis_PCC-6803"  
#[7] "g__Calothrix_KVSF5"          "g__Phormidium_SAG_37.90"     "g__Aphanizomenon_NIES81"    
#[10] "g__Cuspidothrix_LMECYA_163"  "g__Snowella_0TU37S04"        "g__Pseudanabaena_PCC-7429"  
#[13] "g__Gloeotrichia_SAG_32.84"   "g__Leptolyngbya_PCC-6306"    "g__LB3-76"                  
#[16] "g__Leptolyngbyaceae"         "g__Leptolyngbya_SAG_2411"    "g__Microcystis_PCC-7914"    
#[19] "g__Planktothrix_NIVA-CYA_15" "g__MIZ36"                    "g__Nodosilinea_PCC-7104"    
#[22] "g__JSC-12"                   "g__Schizothrix_LEGE_07164"   "g__Rivularia_PCC-7116"      
#[25] "g__Microcystaceae"           "g__Limnothrix"               "g__SepB-3"                  
#[28] "g__CENA359"                  "g__Tychonema_CCAP_1459-11B"  "g__Planktothricoides_SR001" 
#[31] "g__Richelia_HH01"            "g__Nodosilineaceae"          "g__Cyanothece_PCC-8801"   

# use updated (2025) taxa names determined through literature review to clean up remaining names
name_change<-read_excel(paste0(dir,"/Taxa_name_changes2.xlsx"))
name_change$Genus<-name_change$Orig_Genus
name_change$Family<-name_change$Orig_Family

seq_data<-name_change[c(4:9)]
seq_data<-rename(seq_data, "Genus"="New_Genus")
seq_data<-unique(seq_data)

merged_table <- merge(table_rename, name_change, by=c("Genus","Family"), all=TRUE) 
#clean up table
merged_table$Genus<-merged_table$New_Genus
merged_table$Family<-merged_table$New_Family
merged_table<-merged_table[c(1:2,4:49)]

genus_sample_reads<-merged_table[c(1, 3:48)] %>% 
  group_by(Genus) %>% 
  summarise(across(c(1:46), sum))

family_sample_reads<-merged_table[c(2:48)] %>% 
  group_by(Family) %>% 
  summarise(across(c(1:46), sum))

genus_seq<- left_join(genus_sample_reads, seq_data, by="Genus")

family_pivot<-family_sample_reads %>%
  pivot_longer(cols=c(2:47), names_to="Sample", values_to="Abundance") 
family_pivot<-family_pivot[family_pivot$Abundance!=0,]

genus_pivot<-genus_seq %>%
  pivot_longer(cols=c(2:47), names_to="Sample", values_to="Abundance")
#drop zeroes
genus_pivot<-genus_pivot[genus_pivot$Abundance!=0,]

# Merge the abundance data with the total sample read count
genus_merge <- left_join(genus_pivot, library_size, by="Sample") 
family_merge <- merge(family_pivot, library_size, by="Sample", all=TRUE) 
#### Proportion normalization ####
# Equation:
#       Relative Abundance(ASV1) = N(ASV1) / (N(ASV1)+N(ASV2)+...N(ASVn))

# General cyanobacteria
overall_cyano_abund$Prop_abund_cyano<-overall_cyano_abund$Cyanobacteria /overall_cyano_abund$bact_library_size
overall_cyano_abund$Prop_abund_bact<-overall_cyano_abund$'Other Bacteria' / overall_cyano_abund$bact_library_size

# Genus
genus_merge$genus_prop_abund<-genus_merge$Abundance / genus_merge$cyano_library_size

# Family
family_merge$family_prop_abund<-family_merge$Abundance / family_merge$cyano_library_size

#### Check normality of abundance data ####
# check normality
wq_cyano_abund <- merge(overall_cyano_abund, wq_data, by="Sample", all=TRUE) 
res_aov <- aov(Prop_abund_cyano ~ Date, data = wq_cyano_abund)
par(mfrow = c(1, 2)) 
# histogram
hist(res_aov$residuals)
# QQ-plot
qqPlot(res_aov$residuals, add.line = TRUE)
# not completely normal, more normal than data that was not summed

# try log transformation
res_aov <- aov(log10(Prop_abund_cyano+0.01) ~ Date,  data = wq_cyano_abund)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# This is much more normal - use log transformation moving forward
#      for parametric tests

#double check genus and family data
wq_genus_merge <- merge(genus_merge, wq_data, by="Sample", all=TRUE) 
wq_genus_drop<-wq_genus_merge[wq_genus_merge$Abundance!=0,]
res_aov <- aov(genus_prop_abund ~ Date, data = wq_genus_drop)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals, add.line = TRUE)

res_aov <- aov(log10(genus_prop_abund+0.01) ~ Date,  data = wq_genus_drop)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

wq_family_merge <- merge(family_merge, wq_data, by="Sample", all=TRUE) 
wq_family_drop<-wq_family_merge[wq_family_merge$Abundance!=0,]
res_aov <- aov(family_prop_abund ~ Date, data = wq_family_drop)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals, add.line = TRUE)

res_aov <- aov(log10(family_prop_abund+0.01) ~ Date,  data = wq_family_drop)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# Apply transformations
family_merge$logP_Abund_family<-log10(family_merge$family_prop_abund+0.01)
genus_merge$logP_Abund_genus<-log10(genus_merge$genus_prop_abund+0.01)
overall_cyano_abund$logP_Abund_cyano<-log10(overall_cyano_abund$Prop_abund_cyano+0.01)
overall_cyano_abund$logP_Abund_bact<-log10(overall_cyano_abund$Prop_abund_bact+0.01)

#### Output dataframes ####
# Merge WQ data with abundance data
wq_genus_merge <- merge(genus_merge, wq_data, by="Sample", all=TRUE) 
wq_family_merge <- merge(family_merge, wq_data, by="Sample", all=TRUE) 
wq_cyano_abund <- merge(overall_cyano_abund, wq_data, by="Sample", all=TRUE) 

# re-order data frames
wq_genus_merge <- wq_genus_merge %>%
  select(Sample, Date, IDL, everything())

wq_family_merge <- wq_family_merge %>%
  select(Sample, Date, IDL, everything())

wq_cyano_abund<- wq_cyano_abund %>%
  select(Sample, Date, IDL, everything())

# Abundance data only
write_xlsx(overall_cyano_abund, path = paste0(dir,"/overall_cyano_ASV_abundance.xlsx"))
write_xlsx(genus_merge, path = paste0(dir,"/genus_cyano_ASV_abundance.xlsx"))
write_xlsx(family_merge, path = paste0(dir,"/family_cyano_ASV_abundance.xlsx"))

# Abundance + WQ data
write_xlsx(wq_cyano_abund, path = paste0(dir,"/wq_overall_cyano_ASV_abundance.xlsx"))
write_xlsx(wq_genus_merge, path = paste0(dir,"/wq_genus_cyano_ASV_abundance.xlsx"))
write_xlsx(wq_family_merge, path = paste0(dir,"/wq_family_cyano_ASV_abundance.xlsx"))
