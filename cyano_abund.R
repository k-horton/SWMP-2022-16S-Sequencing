library(decontam)
library(file2meco)
library(microeco)
library(dplyr)
library(ggplot2)
library(phyloseq)

# First, load sequence data and remove contaminant samples as determined previously 
    # (see files decontam_prev.R, decontam_freq.R)
    # If you already have the ps.prev.noncontam object from previous scripts, skip next 
    # section and continue to "If you already have ps.prev.noncontam object"

#### If you don't have ps.prev.noncontam object: ####

#create microtable as before
# Assign current working directory to 'dir'
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

# Skip the next section "If you already have ps.prev.noncontam object" and
 # continue on at section "Now continue with "no_ctrl" object"

#### If you already have ps.prev.noncontam object: ####

# convert phyloseq object back to a microtable:
df <- phyloseq2meco(ps.prev.noncontam)

no_ctrl<- df
no_ctrl$sample_table<-subset(no_ctrl$sample_table, Control == "FALSE")
no_ctrl$tidy_dataset()


#### Now continue with "no_ctrl" object: ####

# Calculate abundance for all bacteria in the samples
bact<-clone(no_ctrl)

bact$tax_table <- subset(bact$tax_table, Kingdom == "k__Bacteria")
all_bact_abund <- trans_abund$new(dataset = bact, taxrank = "Genus") 


# make a list of taxa that are part of the cyanobacteriia class
cyano<-clone(no_ctrl)
cyano$tax_table <- subset(cyano$tax_table, Class == "c__Cyanobacteriia")
unique(cyano$tax_table$Genus)
# [1] "g__"                         "g__Cyanobium_PCC-6307"       "g__Nodosilineaceae"         
# [4] "g__Nodosilinea_PCC-7104"     "g__Schizothrix_LEGE_07164"   "g__CENA359"                 
# [7] "g__Leptolyngbya_ANT.L52.2"   "g__Calothrix_KVSF5"          "g__Planktothricoides_SR001" 
# [10] "g__Planktothrix_NIVA-CYA_15" "g__Tychonema_CCAP_1459-11B"  "g__Rivularia_PCC-7116"      
# [13] "g__Richelia_HH01"            "g__Aphanizomenon_NIES81"     "g__Cuspidothrix_LMECYA_163" 
# [16] "g__Gloeotrichia_SAG_32.84"   "g__LB3-76"                   "g__Phormidesmis_ANT.L52.6"  
# [19] "g__JSC-12"                   "g__MIZ36"                    "g__Phormidium_SAG_37.90"    
# [22] "g__RD011"                    "g__Snowella_0TU37S04"        "g__Synechocystis_PCC-6803"  
# [25] "g__Microcystis_PCC-7914"     "g__Leptolyngbya_PCC-6306"    "g__Leptolyngbyaceae"        
# [28] "g__Limnothrix"               "g__Leptolyngbya_SAG_2411"    "g__Cyanothece_PCC-8801"     
# [31] "g__Microcystaceae"           "g__SepB-3"                   "g__Pseudanabaena_PCC-7429"  

# generate simplified taxa names for cyanobacteria (allows filtering of bacterial abundance data)
cyano_abund <- trans_abund$new(dataset = cyano, taxrank = "Genus") 
unique(cyano_abund$data_abund$Taxonomy)
# [1] "Aphanizomenon_NIES81"     "CENA359"                  "Calothrix_KVSF5"         
# [4] "Cuspidothrix_LMECYA_163"  "Cyanobium_PCC-6307"       "Cyanothece_PCC-8801"     
# [7] "Gloeotrichia_SAG_32.84"   "JSC-12"                   "LB3-76"                  
# [10] "Leptolyngbya_ANT.L52.2"   "Leptolyngbya_PCC-6306"    "Leptolyngbya_SAG_2411"   
# [13] "Leptolyngbyaceae"         "Limnothrix"               "MIZ36"                   
# [16] "Microcystaceae"           "Microcystis_PCC-7914"     "Nodosilinea_PCC-7104"    
# [19] "Nodosilineaceae"          "Phormidesmis_ANT.L52.6"   "Phormidium_SAG_37.90"    
# [22] "Planktothricoides_SR001"  "Planktothrix_NIVA-CYA_15" "Pseudanabaena_PCC-7429"  
# [25] "RD011"                    "Richelia_HH01"            "Rivularia_PCC-7116"      
# [28] "Schizothrix_LEGE_07164"   "SepB-3"                   "Snowella_0TU37S04"       
# [31] "Synechocystis_PCC-6803"   "Tychonema_CCAP_1459-11B"  "unidentified"  
cyano_list<-cyano_abund$data_abund$Taxonomy #1518 rows
un_cyano_list<-unique(cyano_list) #33 unique names

#make a list of all bacteria w/o cyanobacteria included
bact_list<-clone(no_ctrl)
bact_list$tax_table <- subset(bact_list$tax_table, Kingdom == "k__Bacteria")
bact_list$tax_table <- subset(bact_list$tax_table, Class != "c__Cyanobacteriia")
bact_list$tidy_dataset()

bact_abund <- trans_abund$new(dataset = bact_list, taxrank = "Genus") 
bact_list<-bact_abund$data_abund$Taxonomy #56097 rows
un_bact_list<-unique(bact_list) #813 unique names


# Use lists to rename taxa in original dataset
cyano_label<-rep(list("Cyanobacteria"), 33) #33 unique names
bact_label<-rep(list("Other Bacteria"), 813) #813 unique names

all_bact_abund$data_abund$Taxonomy<-factor(all_bact_abund$data_abund$Taxonomy,
                                       levels=c(un_bact_list, un_cyano_list),
                                       labels = c(bact_label, cyano_label))
table(all_bact_abund$data_abund$Taxonomy)
# Other Bacteria  Cyanobacteria 
# 56097           2208

# set colours for each group
my_colours_cyano<-c("Cyanobacteria" = "turquoise",
                   "Other Bacteria"="grey40")

# set order of dates on plot
date<-c("June 23rd",
     "July 20th",
     "Aug 3rd",
     "Aug 23rd",
     "Aug 31st",
     "Sept 27th")

#shows each observation independently - probably best representation of the data
ggplot(all_bact_abund$data_abund, aes(x= IDL,y= Abundance, 
                          fill = Taxonomy)) + 
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=my_colours_cyano, name = "Group") + 
  xlab("SWMP") +
  ylab("Relative Abundance") +
  facet_wrap(factor(all_bact_abund$data_abund$Date, level=date))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(strip.background = element_rect(fill="gray85"),
        panel.margin = unit(0.3, "lines"),
        axis.text.x=element_text(angle=0,vjust=0.5),
        text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        plot.margin=margin(0.5,0.5,0.5,0.5,unit="cm"))


# Trim the plot to 50% abundance so it's easier to see abundance of cyanobacteria
ggplot(all_bact_abund$data_abund, aes(x= IDL,y= Abundance, 
                                      fill = Taxonomy)) + 
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=my_colours_cyano, name = "Group") + 
  xlab("SWMP") +
  ylab("Relative Abundance") +
  facet_wrap(factor(all_bact_abund$data_abund$Date, level=date))+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(strip.background = element_rect(fill="gray85"),
        panel.margin = unit(0.3, "lines"),
        axis.text.x=element_text(angle=0,vjust=0.5),
        text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        plot.margin=margin(0.5,0.5,0.5,0.5,unit="cm"))+
  coord_cartesian(ylim=c(0,0.50))
