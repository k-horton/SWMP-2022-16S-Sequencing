##### Decontam Prevalence ####
# Prevalence-based contaminant identification is preferred for low biomass environments, 
# and relies on the use of a negative control sample.

library(decontam)
library(file2meco)
library(ggplot2)
library(microeco)
library(phyloseq)

#create microtable as before
{dir<-getwd()

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
}

# First drop the positive control sample from the microtable
table2<-clone(tcontam)
table2$sample_table<-subset(table2$sample_table, table2$sample_table$Control!="Pos")

# Remove pollution (mitochondrial and chloroplast sequences) from the microtable
table2$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# Total 366 features are removed from tax_table ...

# convert microtable to phyloseq object
physeq <- meco2phyloseq(table2)
physeq

# In the phyloseq object, "Sample_or_Control" is the sample variable that holds the 
# negative control information. 

sample_data(physeq)$is.neg <- sample_data(physeq)$Control == "Neg"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold = 0.10)
table(contamdf.prev$contaminant)
# FALSE  TRUE 
#  8657   3  
# There are 3 contaminants detected
head(which(contamdf.prev$contaminant))
## [1]  1979  6275  6490

# Subset only the sequences designated as contaminants
prev.contam<-subset(contamdf.prev, contamdf.prev$contaminant == TRUE)
prev.contam
#                                       freq prev p.freq     p.prev          p contaminant
# 101637aba474a6614786486ad42acb11 0.008030306   12     NA 0.08571429 0.08571429        TRUE
# 0a1f5ed70dcb37a86924dbd5fde981af 0.005845732   11     NA 0.07857143 0.07857143        TRUE
# d1d2369e5d24e487e8c865eff9161dcb 0.002636087    6     NA 0.04285714 0.04285714        TRUE

# Make a list of contaminant sequence IDs
prev.contam.list<-row.names(prev.contam)
prev.contam.list
# [1] "101637aba474a6614786486ad42acb11" "0a1f5ed70dcb37a86924dbd5fde981af"
# [3] "d1d2369e5d24e487e8c865eff9161dcb"

# Check what the matching taxa are for the contaminant sequence IDs
prev.contam.taxa<-clone(table2)
prev.contam.taxa$otu_table<-prev.contam.taxa$otu_table[row.names(prev.contam.taxa$otu_table) %in% c(prev.contam.list),]
prev.contam.taxa$tidy_dataset()
# 45 samples with 0 abundance are removed from the otu_table ...
prev.taxa.name<-trans_abund$new(dataset = prev.contam.taxa, taxrank="Genus")
unique(prev.taxa.name$data_abund$Taxonomy)
# [1] "Enhydrobacter" "Pelomonas"     "Streptococcus"

# Since there are only 3 contaminants identified with a 0.1 threshold, try using a more aggressive threshold to see if 
# more contaminants are identified.
# When threshold = 0.5, this will identify contaminants that are more prevalent in negative controls than in true samples. 

contamdf.prev05 <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
# FALSE  TRUE 
#  8657   3     
# Even with a more aggressive threshold, the same number of taxa are designated as contaminants. 

# Look at the number of times several of these taxa were observed in negative controls and positive samples.
# Make phyloseq object of presence-absence in negative controls and true samples
physeq.pa <- transform_sample_counts(physeq, function(abund) 1*(abund>0))
physeq.pa.neg <- prune_samples(sample_data(physeq.pa)$Control == "Neg", physeq.pa)
physeq.pa.pos <- prune_samples(sample_data(physeq.pa)$Control == "FALSE", physeq.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(physeq.pa.pos), pa.neg=taxa_sums(physeq.pa.neg), contaminant=contamdf.prev$contaminant)

# Plot the prevalence of taxa in true samples against prevalence in negative control
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


# Remove contaminants from dataframe
ps.prev.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)
ps.prev.noncontam
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8657 taxa and 70 samples ]
# sample_data() Sample Data:       [ 70 samples by 46 sample variables ]
# tax_table()   Taxonomy Table:    [ 8657 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8657 tips and 8598 internal nodes ]
# refseq()      DNAStringSet:      [ 8657 reference sequences ]

# Moving forward with analysis, the ps.prev.noncontam object will be used for statistical analyses.