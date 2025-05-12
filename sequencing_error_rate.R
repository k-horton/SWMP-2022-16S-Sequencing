library(microeco)

# Use the same microtable generated previously, but clone and subset only the positive control
pos_ctrl <- clone(table1)
pos_ctrl$sample_table <- subset(pos_ctrl$sample_table, Control == "Pos")

# use tidy_dataset to trim all the basic files
pos_ctrl$tidy_dataset()
pos_ctrl
#microtable-class object:
  # sample_table have 1 rows and 45 columns
  # otu_table have 17 rows and 1 columns
  # tax_table have 17 rows and 7 columns
  # phylo_tree have 17 tips
  # rep_fasta have 17 sequences

# Calculate abundance of the positive control at the Genus level 
pos_abund <- trans_abund$new(dataset = pos_ctrl, taxrank = "Genus")
pos_abund

#data_taxanames:  Bacillus, Listeria, Limosilactobacillus, Staphylococcus, Enterococcus, 
#                 Escherichia-Shigella, Pseudomonas, Paenibacillus 

# Note that the manufacturer states that the positive control contains the following taxa:
# Pseudomonas aeruginosa, Bacillus subtilis, Escherichia coli, Salmnonella enterica
# Lactobacillus fermentum, Enterococcus faecalis, Staphylococcus aureus,
# Saccharomyces cerevisiae, Cryptococcus neoformans

# However, with the V3-V4 16S sequencing kit used, 
#    Saccharomyces and Cryptococcus (Eukaryotes) will not be detected

# Generate a pie chart to visualize the results
pos_abund$plot_pie(facet_nrow = 1, add_label = TRUE, legend_text_italic = TRUE)

# Now, create a dataframe with only the taxonomy and experimental % abundance data
Taxonomy<-pos_abund$data_abund$Taxonomy
e_Abundance<-pos_abund$data_abund$Abundance
experimental_abundance<-data.frame(Taxonomy, e_Abundance)
experimental_abundance
#             Taxonomy e_Abundance
# 1             Bacillus 19.35023106
# 2         Enterococcus 10.71278532
# 3 Escherichia-Shigella  9.55608458
# 4  Limosilactobacillus 16.03136816
# 5             Listeria 16.11398964
# 6        Paenibacillus  0.00700182
# 7          Pseudomonas  3.72356813
# 8       Staphylococcus 14.86206414
# 9         unidentified  9.64290716

# Some notes on the experimental taxonomy:
# 1. Limosilactobacillus was previously a part of the genus Lactobacillus, so it will be designated
# as Lactobacillus for determining error rate

# 2. There is Escherichia coli but not Shigella in the positive control. 
# Escherichia-Shigella will be assumed to be composed of Escherichia coli only. 

# 3. There is an unidentified taxa. The only taxa missing that we know should be present is Salmonella.
# The expected abundance of Salmonella closely matches the experimental abundance of the unidentified taxa, 
# so it will be designated as Salmonella for determining the error rate.

# 4. Paenibacillus was not one of the expected taxa. It would likely be classified as Bacillus, so add its 
# abundance to the Bacillus abundance (it has such low abundance it could be probably be removed altogether)

# Re-organize the experimental dataframe based on above notes:
# Re-name Escherichia-Shigella and Unidentified
experimental_abundance$Taxonomy<-factor(experimental_abundance$Taxonomy,
                            levels = c("Bacillus", "Enterococcus", "Escherichia-Shigella",
                                       "Limosilactobacillus", "Listeria", "Paenibacillus",
                                       "Pseudomonas", "Staphylococcus", "unidentified"),
                            labels = c("Bacillus", "Enterococcus", "Escherichia",
                                       "Lactobacillus", "Listeria", "Paenibacillus",
                                       "Pseudomonas", "Staphylococcus", "Salmonella"))

# Add Paenibacillus abundance to Bacillus
experimental_abundance$e_Abundance[1] <- experimental_abundance$e_Abundance[1] + experimental_abundance$e_Abundance[6]

# Remove Paenibacillus
experimental_abundance <- experimental_abundance[-6, ]
experimental_abundance
#         Taxonomy e_Abundance
# 1       Bacillus   19.357233
# 2   Enterococcus   10.712785
# 3    Escherichia    9.556085
# 4  Lactobacillus   16.031368
# 5       Listeria   16.113990
# 7    Pseudomonas    3.723568
# 8 Staphylococcus   14.862064
# 9     Salmonella    9.642907

# create dataframe with theoretical % abundance for the positive control 
#     (from manufactor's specifications)
Taxonomy<-c("Bacillus", "Enterococcus","Escherichia","Lactobacillus", "Listeria","Pseudomonas", "Staphylococcus", "Salmonella")
t_Abundance<-c(17.4, 9.9, 10.1, 18.4, 14.1, 4.2, 15.5, 10.4)

theo_abundance<-data.frame(Taxonomy, t_Abundance)
theo_abundance
#        Taxonomy t_Abundance
# 1       Bacillus        17.4
# 2   Enterococcus         9.9
# 3    Escherichia        10.1
# 4  Lactobacillus        18.4
# 5       Listeria        14.1
# 6    Pseudomonas         4.2
# 7 Staphylococcus        15.5
# 8     Salmonella        10.4

# Join both datasets
abundance<- merge(experimental_abundance,theo_abundance,by="Taxonomy")
abundance
#        Taxonomy e_Abundance t_Abundance
# 1       Bacillus   19.357233        17.4
# 2   Enterococcus   10.712785         9.9
# 3    Escherichia    9.556085        10.1
# 4  Lactobacillus   16.031368        18.4
# 5       Listeria   16.113990        14.1
# 6    Pseudomonas    3.723568         4.2
# 7     Salmonella    9.642907        10.4
# 8 Staphylococcus   14.862064        15.5

# calculate deviation between experimental and theoretical abundance
abundance$Difference<- abs(abundance$e_Abundance - abundance$t_Abundance)
abundance

#        Taxonomy e_Abundance t_Abundance Difference
# 1       Bacillus   19.357233        17.4  1.9572329
# 2   Enterococcus   10.712785         9.9  0.8127853
# 3    Escherichia    9.556085        10.1  0.5439154
# 4  Lactobacillus   16.031368        18.4  2.3686318
# 5       Listeria   16.113990        14.1  2.0139896
# 6    Pseudomonas    3.723568         4.2  0.4764319
# 7     Salmonella    9.642907        10.4  0.7570928
# 8 Staphylococcus   14.862064        15.5  0.6379359

# calculate mean of the difference in abundance
mean(abundance$Difference)
#[1] 1.196002

# So, the error rate is only +- 1.20%, this is within the manufacturers error (< 15%)
