### Prepare dataframes
library(decontam)
library(file2meco)
library(microeco)
library(writexl)
library(readxl)
library(tidyr)
library(dplyr)
library(EnvStats)

#### Prepare Taxonomy Data ####
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
# Calculate abundance for all bacteria in the samples
{bact<-clone(no_ctrl)

bact$tax_table <- subset(bact$tax_table, Kingdom == "k__Bacteria")
all_bact_abund <- trans_abund$new(dataset = bact, taxrank = "Genus",
                                  delete_taxonomy_lineage = FALSE,
                                  delete_taxonomy_prefix = TRUE) 
# separate taxonomy
s_all_bact_abund<-separate_wider_delim(all_bact_abund$data_abund, cols = Taxonomy, delim = "|", 
                                       names = c("Kingdom", "Phylum", "Class", "Order",
                                                 "Family", "Genus"))
}
# Subset to just cyanobacteria
{cyano_abund <- s_all_bact_abund[which(s_all_bact_abund$Class=="Cyanobacteriia"), ]
# create dataframe simplified to class
gen_cyano<-cyano_abund[c(3,7:55)]

con_data1<- cyano_abund[c(3,7:55)]%>%
  group_by(Sample) %>%
  summarize(across(Abundance, sum), .groups="drop")

con_data2<-unique(cyano_abund[c(7, 12:53)])

gen_cyano<-merge(con_data1, con_data2, by = "Sample")}
dir<-getwd()
excel_output_path <- paste0(dir,"/class_lev_cyano.xlsx")
write_xlsx(gen_cyano, path = excel_output_path)

# create a dataframe with all class, order, family, genus info
spec_cyano<-cyano_abund[c(3:6,8,10,12,13)]

#drop zero values
spec_cyano<- spec_cyano[which(spec_cyano$Abundance!=0), ]

excel_output_path <- paste0(dir,"/all_lev_cyano.xlsx")
write_xlsx(spec_cyano, path = excel_output_path)

#### Read in WQ data and merge w/ other microbial data  ####
dir<-getwd()
seq<-read_excel(paste0(dir,"/class_lev_cyano.xlsx"))
df1<-read_excel(paste0(dir,"/swmp.xlsx"))

seq$key<-paste(seq$Date,"-" ,seq$IDL) 
df1$key<-paste(df1$Date,"-" ,df1$IDL) 

# merge seq data w/ swmp data
seq_sub<-seq[c(2,41,42,44,45)]
total_cyano<-merge(x = df1, y = seq_sub, by = "key", all.x = TRUE)

#### Create a 2nd dataframe with Genus-specific abundance data ####
df2<-total_cyano
seq2<-read_excel(paste0(dir,"/all_lev_cyano.xlsx"))

# merge the WQ + micro data with the sequencing data
seq2$key<-paste(seq2$Date,"-" ,seq2$IDL) 
df2$key<-paste(df2$Date,"-",df2$IDL) 
#rename abundance columns in both dataframes
df2<-rename(df2,t_cyano_abundance = Abundance)
seq2<-rename(seq2,genus_abundance = Abundance)

seq_sub<-seq2[c(1:6, 9)]
genus_cyano<-merge(x = df2, y = seq_sub, by = "key", all.x = TRUE)

#### Assign labels to categorical variables ####
total_cyano$Storm_Base<-factor(total_cyano$Storm_Base ,
                        levels = c("Storm","Base"),
                        labels = c("Stormflow","Low Flow"))

total_cyano$Site<-factor(total_cyano$Site ,
                  levels = c("in","pond","out"),
                  labels = c("Inflow","Pond","Outflow"))

total_cyano$Dredge_Cat<-factor(total_cyano$Dredge_Cat,
                        levels = c("Recent", "Needs Dredging"),
                        labels = c("Recently Dredged", "Needs Dredging"))

total_cyano$Pre_2003<-factor(total_cyano$Pre_2003,
                      levels = c("Yes", "No"),
                      labels = c("Built Pre-2003", "Built Post-2003"))

genus_cyano$Storm_Base<-factor(genus_cyano$Storm_Base ,
                       levels = c("Storm","Base"),
                       labels = c("Stormflow","Low Flow"))

genus_cyano$Site<-factor(genus_cyano$Site ,
                 levels = c("in","pond","out"),
                 labels = c("Inflow","Pond","Outflow"))

genus_cyano$Dredge_Cat<-factor(genus_cyano$Dredge_Cat,
                       levels = c("Recent", "Needs Dredging"),
                       labels = c("Recently Dredged", "Needs Dredging"))

genus_cyano$Pre_2003<-factor(genus_cyano$Pre_2003,
                     levels = c("Yes", "No"),
                     labels = c("Built Pre-2003", "Built Post-2003"))


#### Calculate additional variables for both datasets ####
# Dissolved inorganic nitrogen
total_cyano$DIN<-(total_cyano$Nitrate_mgL+total_cyano$Nitrite_mgL+total_cyano$Ammonia_Ammonium_mgL)
# Separate NH4 + NH3
total_cyano$pka <- 0.09018 + (2729.92/(total_cyano$Temp_C +273.15))
total_cyano$f_NH3 <- 1/(10^(total_cyano$pka - total_cyano$pH)+1)
total_cyano$NH3_mgL <- total_cyano$Ammonia_Ammonium_mgL*total_cyano$f_NH3
total_cyano$NH4_mgL<-total_cyano$Ammonia_Ammonium_mgL - total_cyano$NH3_mgL
# Nitrogen species : TN ration
total_cyano$NO3_TN<-total_cyano$Nitrate_mgL / total_cyano$Total_N_mgL
total_cyano$NO2_TN<-total_cyano$Nitrite_mgL / total_cyano$Total_N_mgL
total_cyano$ON_TN<-total_cyano$Organic_N_mgL / total_cyano$Total_N_mgL
total_cyano$NH3NH4_TN<-total_cyano$Ammonia_Ammonium_mgL / total_cyano$Total_N_mgL
total_cyano$NH3_TN <- total_cyano$NH3_mgL / total_cyano$Total_N_mgL
total_cyano$NH4_TN <- total_cyano$NH4_mgL / total_cyano$Total_N_mgL
# TOSS : TSS ratio
total_cyano$TOSS_TSS <- total_cyano$TOSS_gL / total_cyano$TSS_gL
# TN : TP ratio
total_cyano$TN_TP <- total_cyano$Total_N_mgL / (total_cyano$TP_ugL*0.001) # <- need to convert to mg/L
# change any infinite values generated to NA
is.na(total_cyano) <- sapply(total_cyano, is.infinite)

# Dissolved inorganic nitrogen
genus_cyano$DIN<-(genus_cyano$Nitrate_mgL+genus_cyano$Nitrite_mgL+genus_cyano$Ammonia_Ammonium_mgL)
# Separate NH4 + NH3
genus_cyano$pka <- 0.09018 + (2729.92/(genus_cyano$Temp_C +273.15))
genus_cyano$f_NH3 <- 1/(10^(genus_cyano$pka - genus_cyano$pH)+1)
genus_cyano$NH3_mgL <- genus_cyano$Ammonia_Ammonium_mgL*genus_cyano$f_NH3
genus_cyano$NH4_mgL<-genus_cyano$Ammonia_Ammonium_mgL - genus_cyano$NH3_mgL
# Nitrogen species : TN ration
genus_cyano$NO3_TN<-genus_cyano$Nitrate_mgL / genus_cyano$Total_N_mgL
genus_cyano$NO2_TN<-genus_cyano$Nitrite_mgL / genus_cyano$Total_N_mgL
genus_cyano$ON_TN<-genus_cyano$Organic_N_mgL / genus_cyano$Total_N_mgL
genus_cyano$NH3NH4_TN<-genus_cyano$Ammonia_Ammonium_mgL / genus_cyano$Total_N_mgL
genus_cyano$NH3_TN <- genus_cyano$NH3_mgL / genus_cyano$Total_N_mgL
genus_cyano$NH4_TN <- genus_cyano$NH4_mgL / genus_cyano$Total_N_mgL
# TOSS : TSS ratio
genus_cyano$TOSS_TSS <- genus_cyano$TOSS_gL / genus_cyano$TSS_gL
# TN : TP ratio
genus_cyano$TN_TP <- genus_cyano$Total_N_mgL / (genus_cyano$TP_ugL*0.001) # <- need to convert to mg/L
# change any infinite values generated to NA
is.na(genus_cyano) <- sapply(genus_cyano, is.infinite)

#### Check normality of the data #####
# check normality
res_aov <- aov(Abundance ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
# histogram
hist(res_aov$residuals)
# QQ-plot
qqPlot(res_aov$residuals, add.line = TRUE)
# not completely normal, more normal than data that was not summed

# try log transformation

res_aov <- aov(log10(Abundance+0.1) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# This is much more normal - use log transformation moving forward
#      for parametric tests

# Chlorophyll-a #
res_aov <- aov(Chla_gL ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Chla_gL+0.01) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# This is better - use log transformation 

# Unicellular cyanobacteria cells #
res_aov <- aov(Unicellularcells_L ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Unicellularcells_L+0.1) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
total_cyano_sub<-total_cyano
total_cyano_sub$Unicellularcells_L[total_cyano_sub$Unicellularcells_L==0] <- NA
total_cyano_sub<-total_cyano_sub[complete.cases(total_cyano_sub),]

res_aov <- aov(Unicellularcells_L ~ Date,
               data = total_cyano_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Unicellularcells_L+0.01) ~ Date,
               data = total_cyano_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

#this looks good! Use log transformation


# Colonial cyanobacteria cells #
res_aov <- aov(Colonialcells_L ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Colonialcells_L+0.01) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# This is better - use log transformation 


# Filamentous cyanobacteria cells #
res_aov <- aov(Filamentouscells_L ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Filamentouscells_L+0.01) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# This is questionable, try some other transformations:

# root transformation -> x^(1/2)
res_aov <- aov((Filamentouscells_L+0.01)^(1/2) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# reciprocal transformation -> 1/x
res_aov <- aov(1/(Filamentouscells_L+0.01) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# squared transformation -> x^2
res_aov <- aov((Filamentouscells_L+0.01)^(2) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
total_cyano_sub<-total_cyano
total_cyano_sub$Filamentouscells_L[total_cyano_sub$Filamentouscells_L==0] <- NA
total_cyano_sub<-total_cyano_sub[complete.cases(total_cyano_sub),]

res_aov <- aov(Filamentouscells_L ~ Date,
               data = total_cyano_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Filamentouscells_L+0.01) ~ Date,
               data = total_cyano_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

#this looks good! Use log transformation

# Total cyanobacteria cells #
res_aov <- aov(Totalcells_L ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Totalcells_L+0.01) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# This is better - use log transformation 

# 16S copies / L #
res_aov <- aov(copies_16S_L ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(copies_16S_L+0.01) ~ Date,
               data = total_cyano)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
total_cyano_sub<-total_cyano
total_cyano_sub$copies_16S_L[total_cyano_sub$copies_16S_L==0] <- NA
total_cyano_sub<-total_cyano_sub[complete.cases(total_cyano_sub),]

res_aov <- aov(copies_16S_L ~ Date,
               data = total_cyano_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

res_aov <- aov(log10(copies_16S_L+0.01) ~ Date,
               data = total_cyano_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# this is good, use log transformation


# So, to summarize, the following transformations will be applied when conducting
#   parametric statistical tests:
#   - Abundance         ->  log10(Abundance+0.1)
#   - Chlorophyll-a     ->  log10(Chla_gL+0.01) 
#   - Unicellular cells ->  log10(Unicellularcells_L+0.01)
#   - Colonial cells    ->  log10(Colonialcells_L+0.01) 
#   - Filamentous cells ->  log10(Filamentouscells_L+0.01)
#   - Total cells       ->  log10(Totalcells_L+0.01)
#   - 16S copies        ->  log10(copies_16S_L+0.01) 

# normality will be checked for additional variables as they come up





##### Apply Transformations #####
total_cyano$LogAbs_440<-log10(total_cyano$Abs_440 + 0.001)
total_cyano$LogAbs_750<-log10(total_cyano$Abs_750 + 0.001)
total_cyano$LogTSS<-log10(total_cyano$TSS_gL + 0.001)
total_cyano$LogTOSS<-log10(total_cyano$TOSS_gL + 0.001)
total_cyano$LogCl<-log10(total_cyano$Cl_mgL + 0.1)
total_cyano$LogTP<-log10(total_cyano$TP_ugL + 0.1)
total_cyano$LogTDP<-log10(total_cyano$TDP_ugL + 0.1)
total_cyano$LogTKN<-log10(total_cyano$Total_Kjeldahl_N_mgL + 0.1)
total_cyano$LogNH3NH4<-log10(total_cyano$Ammonia_Ammonium_mgL + 0.01)
total_cyano$LogNO3<-log10(total_cyano$Nitrate_mgL + 0.001)
total_cyano$LogNO2<-log10(total_cyano$Nitrite_mgL + 0.001)
total_cyano$LogON<-log10(total_cyano$Organic_N_mgL + 0.1)
total_cyano$LogTN<-log10(total_cyano$Total_N_mgL + 0.01)
total_cyano$LogDIN<-log10(total_cyano$DIN + 0.001)
total_cyano$LogChla<-log10(total_cyano$Chla_gL + 0.001)
total_cyano$LogNH3<-log10(total_cyano$NH3_mgL + 0.0001)
total_cyano$LogNH4<-log10(total_cyano$NH4_mgL + 0.01)
total_cyano$LogNO3_TN<-log10(total_cyano$NO3_TN + 0.001)
total_cyano$LogNO2_TN<-log10(total_cyano$NO2_TN + 0.01)
total_cyano$LogNH3NH4_TN<-log10(total_cyano$NH3NH4_TN + 0.01)
total_cyano$LogNH3_TN<-log10(total_cyano$NH3_TN + 0.001)
total_cyano$LogNH4_TN<-log10(total_cyano$NH4_TN + 0.01)
total_cyano$LogTN_TP<-log10(total_cyano$TN_TP + 0.1)

total_cyano$logChla<-log10(total_cyano$Chla_gL+0.01)
total_cyano$logAbund<-log10(total_cyano$Abundance + 0.1)
total_cyano$log16S<-log10(total_cyano$copies_16S_L+0.01)
total_cyano$logmcyE<-log10(total_cyano$copies_mcyE_L+0.01)
total_cyano$logU_cyano<-log10(total_cyano$Unicellularcells_L+0.01)
total_cyano$logC_cyano<-log10(total_cyano$Colonialcells_L+0.01)
total_cyano$logF_cyano<-log10(total_cyano$Filamentouscells_L+0.01)
total_cyano$logT_cyano<-log10(total_cyano$Totalcells_L+0.01)


write_xlsx(total_cyano,paste0(dir,"/total_cyano_df.xlsx"))



genus_cyano$LogAbs_440<-log10(genus_cyano$Abs_440 + 0.001)
genus_cyano$LogAbs_750<-log10(genus_cyano$Abs_750 + 0.001)
genus_cyano$LogTSS<-log10(genus_cyano$TSS_gL + 0.001)
genus_cyano$LogTOSS<-log10(genus_cyano$TOSS_gL + 0.001)
genus_cyano$LogCl<-log10(genus_cyano$Cl_mgL + 0.1)
genus_cyano$LogTP<-log10(genus_cyano$TP_ugL + 0.1)
genus_cyano$LogTDP<-log10(genus_cyano$TDP_ugL + 0.1)
genus_cyano$LogTKN<-log10(genus_cyano$Total_Kjeldahl_N_mgL + 0.1)
genus_cyano$LogNH3NH4<-log10(genus_cyano$Ammonia_Ammonium_mgL + 0.01)
genus_cyano$LogNO3<-log10(genus_cyano$Nitrate_mgL + 0.001)
genus_cyano$LogNO2<-log10(genus_cyano$Nitrite_mgL + 0.001)
genus_cyano$LogON<-log10(genus_cyano$Organic_N_mgL + 0.1)
genus_cyano$LogTN<-log10(genus_cyano$Total_N_mgL + 0.01)
genus_cyano$LogDIN<-log10(genus_cyano$DIN + 0.001)
genus_cyano$LogChla<-log10(genus_cyano$Chla_gL + 0.001)
genus_cyano$LogNH3<-log10(genus_cyano$NH3_mgL + 0.0001)
genus_cyano$LogNH4<-log10(genus_cyano$NH4_mgL + 0.01)
genus_cyano$LogNO3_TN<-log10(genus_cyano$NO3_TN + 0.001)
genus_cyano$LogNO2_TN<-log10(genus_cyano$NO2_TN + 0.01)
genus_cyano$LogNH3NH4_TN<-log10(genus_cyano$NH3NH4_TN + 0.01)
genus_cyano$LogNH3_TN<-log10(genus_cyano$NH3_TN + 0.001)
genus_cyano$LogNH4_TN<-log10(genus_cyano$NH4_TN + 0.01)
genus_cyano$LogTN_TP<-log10(genus_cyano$TN_TP + 0.1)

genus_cyano$logChla<-log10(genus_cyano$Chla_gL+0.01)
genus_cyano$logTAbund<-log10(genus_cyano$t_cyano_abundance + 0.1)
genus_cyano$logGAbund<-log10(genus_cyano$genus_abundance + 0.1)
genus_cyano$log16S<-log10(genus_cyano$copies_16S_L+0.01)
genus_cyano$logmcyE<-log10(genus_cyano$copies_mcyE_L+0.01)
genus_cyano$logU_cyano<-log10(genus_cyano$Unicellularcells_L+0.01)
genus_cyano$logC_cyano<-log10(genus_cyano$Colonialcells_L+0.01)
genus_cyano$logF_cyano<-log10(genus_cyano$Filamentouscells_L+0.01)
genus_cyano$logT_cyano<-log10(genus_cyano$Totalcells_L+0.01)


write_xlsx(genus_cyano,paste0(dir,"/genus_cyano_df.xlsx"))
