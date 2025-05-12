library(microeco)
library(dplyr)
library(tidyr)

# Use the same microtable generated previously, but clone and subset only the negative control
neg_ctrl <- clone(table1)
neg_ctrl$sample_table <- subset(neg_ctrl$sample_table, Control == "Neg")
# use tidy_dataset to trim all the basic files
neg_ctrl$tidy_dataset()
neg_ctrl 
#microtable-class object:
  #sample_table have 1 rows and 45 columns
  #otu_table have 10 rows and 1 columns
  #tax_table have 10 rows and 7 columns
  #phylo_tree have 10 tips
  #rep_fasta have 10 sequences

#only 10 sequences are actually in the negative control

# Assess the abundance of the negative control at the Genus level check for foreign contamination
neg_abund <- trans_abund$new(dataset = neg_ctrl, taxrank = "Genus")
neg_abund
# data_taxanames:  Phyllobacterium, Enhydrobacter, Pelomonas, Prevotella_9, Streptococcus

# Generate a pie chart to visualize the results
neg_abund$plot_pie(facet_nrow = 1, add_label = TRUE, legend_text_italic = TRUE)

#Let's check all other samples for these taxa:
#genera:  Phyllobacterium, Enhydrobacter, Pelomonas, Prevotella_9, Streptococcus

# Create a dataframe with only the taxonomy and negative control abundance data
Taxonomy<-neg_abund$data_abund$Taxonomy
Abundance<-neg_abund$data_abund$Abundance
neg_abundance<-data.frame(Taxonomy, Abundance)
neg_abundance
#         Taxonomy   Abundance
# 1   Enhydrobacter 26.70916044
# 2       Pelomonas 16.21926789
# 3 Phyllobacterium 28.71243854
# 4    Prevotella_9 14.87525041
# 5   Streptococcus 13.47295575
# 6    unidentified  0.01092697

# Clone the original microtable again and remove the control samples
test_neg_ctrl <- clone(table1)
test_neg_ctrl$sample_table <- subset(test_neg_ctrl$sample_table, Control == "FALSE")
test_neg_ctrl$tidy_dataset()

# Calculate the abundance for all samples (except controls)
sample_abund <- trans_abund$new(dataset = test_neg_ctrl, taxrank = "Genus") #59271 observations

# Drop all observations with abundance of 0
sample_abund$data_abund$Abundance[sample_abund$data_abund$Abundance==0] <- NA
sample_abund$data_abund<-drop_na(sample_abund$data_abund, Abundance) #5036 observations

# Now, check if the genera from our negative control are in our actual samples
# It is not possible to determine if the unidentified taxa are the same between the negative control and the 
# actual samples, and they will be removed during analysis anyways, so don't include 'unidentified'
# taxa here
sample_abund$data_abund<- subset(sample_abund$data_abund, (Taxonomy == "Enhydrobacter" |
                                                             Taxonomy == "Pelomonas"|
                                                             Taxonomy == "Phyllobacterium" |
                                                             Taxonomy == "Prevotella_9" |
                                                             Taxonomy == "Streptococcus"))
sample_abund$data_abund #45 observations
# generate a list of samples that contain the taxa
sample_abund_list<-unique(sample_abund$data_abund$Sample)
sample_abund_list
#[1] "F1_SWMP1_S1"   "F1_SWMP7_S7"   "F1_SWMP10_S10" "F2_SWMP1_S13"  "F2_SWMP8_S20"  "F2_SWMP12_S24"
#[7] "F4_SWMP2_S38"  "F4_SWMP6_S42"  "F5_SWMP1_S49"  "F5_SWMP12_S60" "F1_SWMP8_S8"   "F3_SWMP2_S26" 
#[13] "F5_SWMP3_S51"  "F6_SWMP3_S63"  "F6_SWMP9_S68"  "F1_SWMP4_S4"   "F2_SWMP2_S14"  "F5_SWMP10_S58"
#[19] "F6_SWMP4_S64"  "F6_SWMP6_S65"  "F6_SWMP11_S70" "F1_SWMP2_S2"   "F1_SWMP5_S5"   "F1_SWMP6_S6"  
#[25] "F1_SWMP11_S11" "F1_SWMP12_S12" "F2_SWMP9_S21"  "F3_SWMP1_S25"  "F3_SWMP6_S30"  "F4_SWMP7_S43" 
#[31] "F4_SWMP9_S45"  "F5_SWMP4_S52"  "F5_SWMP5_S53" 

# If abundance of negative control taxa is very low, it should not impact results too much
# check for samples that have abundance of negative control taxa > 1%
above_1<-clone(sample_abund)
above_1$data_abund<- subset(above_1$data_abund, Abundance > 1)
above_1$data_abund # 10 observations
#        Sample     Taxonomy Abundance
# 1687 F2_SWMP12_S24 Enhydrobacter  1.184939
# 1688  F4_SWMP2_S38 Enhydrobacter 13.109005
# 1690  F5_SWMP1_S49 Enhydrobacter  3.424357
# 3263  F5_SWMP3_S51     Pelomonas  1.469981
# 3456  F6_SWMP6_S65  Prevotella_9  1.433107
# 4367   F1_SWMP2_S2 Streptococcus  1.623551
# 4368   F1_SWMP4_S4 Streptococcus  1.538529
# 4373 F1_SWMP11_S11 Streptococcus  1.306731
# 4377  F3_SWMP1_S25 Streptococcus  1.382801
# 4381  F5_SWMP1_S49 Streptococcus 44.084275

# Now, subset orignal dataset (table1) with only the samples listed above, then 
# Visualize results on a pie chart
above_1_plot<-clone(table1)
above_1_plot$tax_table <- subset(above_1_plot$tax_table, Kingdom == "k__Bacteria")
# Calculate abundance
above_1_plot<-trans_abund$new(dataset = above_1_plot, taxrank = "Genus", ntaxa=50)

# Drop all observations with abundance of 0
above_1_plot$data_abund$Abundance[above_1_plot$data_abund$Abundance==0] <- NA
above_1_plot$data_abund<-drop_na(above_1_plot$data_abund, Abundance)


# Subset samples
above_1_plot$data_abund<-subset(above_1_plot$data_abund, (Sample == "F2_SWMP12_S24" |
                                                            Sample == "F4_SWMP2_S38"|
                                                            Sample == "F5_SWMP1_S49"|
                                                            Sample == "F5_SWMP3_S51"|
                                                            Sample == "F6_SWMP6_S65"|
                                                            Sample == "F1_SWMP2_S2" |
                                                            Sample == "F1_SWMP4_S4" |
                                                            Sample == "F1_SWMP11_S11" |
                                                            Sample == "F3_SWMP1_S25" ))

# Visualize pie chart
above_1_plot$plot_pie(facet_nrow = 3, add_label = FALSE, legend_text_italic = TRUE)


t1 <- trans_abund$new(dataset = test_neg_ctrl, taxrank = "Genus", ntaxa = 10)
t1$plot_bar(others_color = "grey70", facet = c("Field_day", "IDL"), xtext_keep = FALSE, 
            legend_text_italic = FALSE, barwidth = 1.5)

# Note that none of the contaminants are part of the Cyanobacteriia class, 
# so there should be minimal impact on conclusions drawn from
# the data with regard to the cyanobacteria community whether the contaminants are removed or not
