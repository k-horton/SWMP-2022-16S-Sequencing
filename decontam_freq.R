#### Decontam Frequency ####
library(decontam)
library(phyloseq)
library(file2meco)
library(microeco)
library(ggplot2)

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

# For frequency based analysis of contamaninants, no controls are used.
# Drop the positive and negative control samples from the microtable. 
table2<-clone(tcontam)
table2$sample_table<-subset(table2$sample_table, table2$sample_table$Control=="FALSE")

# Remove pollution (mitochondrial and chloroplast sequences) from the microtable
table2$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# Total 366 features are removed from tax_table ...

# convert microtable to phyloseq object
physeq <- meco2phyloseq(table2)
physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8653 taxa and 69 samples ]
# sample_data() Sample Data:       [ 69 samples by 45 sample variables ]
# tax_table()   Taxonomy Table:    [ 8653 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8653 tips and 8594 internal nodes ]
# refseq()      DNAStringSet:      [ 8653 reference sequences ]

#In our phyloseq object, "q_conc" is the sample variable that holds the DNA concentration information (ng/L)
contamdf.freq <- isContaminant(physeq, method="frequency", conc="q_conc", threshold=0.10)
# contamdf.freq is a data.frame containing p-values for classifying each sequence as a contaminant or not.
# The contaminant column contains TRUE/FALSE classification values where TRUE indicates the associated sequence feature is a contaminant 
table(contamdf.freq$contaminant)
# FALSE  TRUE 
# 8574    79
# There are 79 sequences determined to be contaminants based on a threshold of 0.10 (if p < 0.10, sequence is classified as a contaminant)

# Subset only the sequences designated as contaminants
freq.contam<-subset(contamdf.freq, contamdf.freq$contaminant == TRUE)
head(freq.contam)
#                                      freq     prev   p.freq    p.prev   p contaminant
# 12f424753516c34ec248982406945990 2.596054e-03   22 0.073364072     NA 0.073364072        TRUE
# b551e403a6113ffd5c8c8ebdc15935d9 1.136352e-04    4 0.076021579     NA 0.076021579        TRUE
# d9efa00159b63446d8d517f4d3e79d93 1.201215e-04    2 0.051410178     NA 0.051410178        TRUE
# 266b16ce9637ef826a0c400922fc1a61 1.045005e-04    2 0.006263649     NA 0.006263649        TRUE
# 1718155e71be8ffea9a18b1d0dddc4f4 1.944074e-05    2 0.049323691     NA 0.049323691        TRUE
# d74ae06c0bcb72d35316529fa0a672e5 2.286096e-03   11 0.023979927     NA 0.023979927        TRUE

# Create a list of contaminant sequence IDs
freq.contam.list<-row.names(freq.contam)
head(freq.contam.list)
# [1] "12f424753516c34ec248982406945990" "b551e403a6113ffd5c8c8ebdc15935d9"
# [3] "d9efa00159b63446d8d517f4d3e79d93" "266b16ce9637ef826a0c400922fc1a61"
# [5] "1718155e71be8ffea9a18b1d0dddc4f4" "d74ae06c0bcb72d35316529fa0a672e5"

# Match taxa to the contaminant sequence IDs
freq.contam.taxa<-clone(table2)
freq.contam.taxa$otu_table<-freq.contam.taxa$otu_table[row.names(freq.contam.taxa$otu_table) %in% c(freq.contam.list),]
freq.contam.taxa$tidy_dataset()
# 4 samples with 0 abundance are removed from the otu_table ...

freq.taxa.name<-trans_abund$new(dataset = freq.contam.taxa, taxrank="Phylum")
unique(freq.taxa.name$data_abund$Taxonomy)
# [1] "Acidobacteriota"   "Actinobacteriota"  "Bacteroidota"      "Bdellovibrionota" 
# [5] "Campylobacterota"  "Chloroflexi"       "Cyanobacteria"     "Dependentiae"     
# [9] "Firmicutes"        "Myxococcota"       "Patescibacteria"   "Planctomycetota"  
# [13] "Proteobacteria"    "Verrucomicrobiota"

freq.taxa.name<-trans_abund$new(dataset = freq.contam.taxa, taxrank="Genus")
unique(freq.taxa.name$data_abund$Taxonomy)
# [1] "Aquicella"               "CL500-3"                 "CPR2"                   
# [4] "Candidatus_Megaira"      "Candidatus_Ovatusbacter" "Chthoniobacter"         
# [7] "Corynebacterium"         "Cyanobium_PCC-6307"      "Dinghuibacter"          
# [10] "Enterococcus"            "Ferribacterium"          "FukuN57"                
# [13] "Gastranaerophilales"     "Haliangium"              "Ilumatobacter"          
# [16] "Legionella"              "Longivirga"              "Micrococcus"            
# [19] "Mycobacterium"           "Paracoccus"              "Phreatobacter"          
# [22] "Pir4_lineage"            "Planktothricoides_SR001" "Prevotella_7"           
# [25] "Pseudarcobacter"         "Pseudolabrys"            "Psychrobacter"          
# [28] "Rhodoferax"              "Saccharimonadales"       "Solitalea"              
# [31] "Staphylococcus"          "Subdoligranulum"         "Subgroup_17"            
# [34] "Sulfurovum"              "Terrimonas"              "Veillonella"            
# [37] "Vermiphilaceae"          "env.OPS_17"              "hgcI_clade"             
# [40] "unidentified"    

# Look at the distributions of a clear non-contaminant and a clear contaminant look like:
# First, output observation numbers for clear non-contaminants:
head(which((contamdf.freq$contaminant == FALSE) & 
             (contamdf.freq$freq > 0.001) & 
             (contamdf.freq$prev > 10)))
# [1]  59  83  84 649 691 699

# Then, output observations numbers for contaminants
which((contamdf.freq$contaminant == TRUE) & (contamdf.freq$prev > 4))
#[1]   71  663 1545 3865 3891 4503 5167 6265

# Using the 61st ASV as the representative non-contaminant and the 8839th ASV 
# as the representative contaminant, plot the frequency of the ASV against the DNA concentration
plot_frequency(physeq, taxa_names(physeq)[c(83,663)], conc="q_conc") + 
  xlab("DNA Concentration (ng/L)")
# The contaminant has a negative linear trend, while the non-contaminant has a random distribution. For non-contaminants, frequency is 
# independent of DNA concentration, while contaminant frequency has an inverse relationship with DNA concentration. 
# This concept is explained in more detail in the paper by Davis et al. 2018.

# For curiosity, determine the taxa names of the clear contaminant and clear true observation we just plotted
# true observation
true<-clone(table2)
true.otu<-taxa_names(physeq)[(83)]
true$otu_table<-subset(true$otu_table, row.names(true$otu_table) == true.otu)
true$tidy_dataset()
#42 samples with 0 abundance are removed from the otu_table ...
true.name<-trans_abund$new(dataset = true, taxrank="Genus")
unique(true.name$data_abund$Taxonomy)
#[1] "Arcobacter"

# contaminant observation
contam<-clone(table2)
contam.otu<-taxa_names(physeq)[(663)]
contam$otu_table<-subset(contam$otu_table, row.names(contam$otu_table) == contam.otu)
contam$tidy_dataset()
#58 samples with 0 abundance are removed from the otu_table ...
contam.name<-trans_abund$new(dataset = contam, taxrank="Genus")
unique(contam.name$data_abund$Taxonomy)
#[1] "Cyanobium_PCC-6307"

# Inspect some more of the ASVs that were classified as contaminants to assess distribution:
set.seed(100)
plot_frequency(physeq, taxa_names(physeq)[sample(which(contamdf.freq$contaminant),30)], conc="q_conc") +
  xlab("DNA Concentration (ng/L)")

# many of the contaminants have a prevalence of only 2 or 3, this can result in falsely 
# designating a sequence as a contaminant

# If we want to remove these contaminants, use the following code:
ps.freq.noncontam <- prune_taxa(!contamdf.freq$contaminant, physeq)
ps.freq.noncontam
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8574 taxa and 69 samples ]
# sample_data() Sample Data:       [ 69 samples by 45 sample variables ]
# tax_table()   Taxonomy Table:    [ 8574 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8574 tips and 8517 internal nodes ]
# refseq()      DNAStringSet:      [ 8574 reference sequences ]

# However, prevalence-based contaminant removal may be better to use in this situation due 
# to low prevalence. 
# Assess sequences based on prevalence and then decide which approach to apply.
