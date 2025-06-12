library(readxl)
library(writexl)
#### Read in previously joined WQ + taxonomy table ####
dir<-getwd()
wq_genus_data<-read_excel(paste0(dir,"/genus_cyano_df.xlsx"))

#### Subset potential microcystin producing taxa #####
#subset taxa
mc <- wq_genus_data[]
mc$tax_table <- subset(mc$tax_table, Class == "c__Cyanobacteriia")

mc$tax_table <- subset(mc$tax_table, (Genus == "g__Aphanizomenon_NIES81"|
                                        Genus == "g__Gloeotrichia_SAG_32.84"|
                                        Genus ==  "g__Microcystis_PCC-7914" |
                                        Genus ==  "g__Phormidium_SAG_37.90" |
                                        Genus == "g__Planktothrix_NIVA-CYA_15"| 
                                        Genus == "g__Pseudanabaena_PCC-7429"))

mc$tidy_dataset()
#60 samples with 0 abundance removed from otu_table
mc
#microtable-class object:
#sample_table have 11 rows and 40 columns
#otu_table have 12 rows and 11 columns
#tax_table have 12 rows and 7 columns
#phylo_tree have 12 tips
#rep_fasta have 12 sequences
mc$sample_names()
#[1] "F1_SWMP1_S1"   "F1_SWMP9_S9"   "F1_SWMP10_S10" "F1_SWMP11_S11" "F2_SWMP5_S17"  "F2_SWMP6_S18" 
#[7] "F2_SWMP8_S20"  "F2_SWMP9_S21"  "F2_SWMP10_S22" "F3_SWMP9_S33"  "F5_SWMP4_S52"

mc$cal_betadiv()

#incude the family Microcystaceae
mc_f <- clone(table1)
mc_f$tax_table <- subset(mc_f$tax_table, Class == "c__Cyanobacteriia")

mc_f$tax_table <- subset(mc_f$tax_table, (Genus == "g__Aphanizomenon_NIES81"|
                                            Genus == "g__Gloeotrichia_SAG_32.84"|
                                            Genus ==  "g__Microcystis_PCC-7914" |
                                            Genus ==  "g__Phormidium_SAG_37.90" |
                                            Genus == "g__Planktothrix_NIVA-CYA_15"| 
                                            Genus == "g__Pseudanabaena_PCC-7429"|
                                            Family == "f__Microcystaceae"))

mc_f$tidy_dataset()
#57 samples with 0 abundance removed from otu_table
mc_f
#microtable-class object:
#sample_table have 14 rows and 44 columns
#otu_table have 18 rows and 14 columns
#tax_table have 18 rows and 7 columns
#phylo_tree have 18 tips
#rep_fasta have 18 sequences
mc_f$sample_names()
#[1] "F1_SWMP1_S1"   "F1_SWMP9_S9"   "F1_SWMP10_S10" "F1_SWMP11_S11" "F2_SWMP5_S17"  "F2_SWMP6_S18" "F2_SWMP6_S19"
#[8] "F2_SWMP8_S20"  "F2_SWMP9_S21"  "F2_SWMP10_S22" "F3_SWMP9_S33"  "F3_SWMP9_S34"  "F3_SWMP9_S46"  "F5_SWMP4_S52"

mc_f$cal_betadiv()

#### MC taxa analysis #####
library(dbplyr)
library(tidyr)
library(ggplot2)
#boxplots
t1 <- trans_abund$new(dataset = mc, taxrank = "Genus", ntaxa = 12) 
t1$plot_box(group = "Field_day", xtext_angle = 30)

write.csv(t1$data_abund, file = "C:/Users/kaitl/Documents/July_Analysis/merge_concat.genus.mc.abund.csv", 
          row.names = TRUE)

#only one bar (mean) for each group, in this case field day
t1 <- trans_abund$new(dataset = mc, taxrank = "Genus", ntaxa = 32, groupmean = "Field_day")
t1$plot_bar(others_color = "grey70", xtext_keep = TRUE, legend_text_italic = FALSE, barwidth = 0.9)



#Try using all cyano data and calling out specific taxa


#potential toxin producers
my_colours_seq<-c("Microcystis_PCC-7914"="indianred4", 
                  "Phormidium_SAG_37.90"= "firebrick4"  , 
                  "Planktothrix_NIVA-CYA_15"="darkred", 
                  
                  #potential nitrogen fixers
                  "Calothrix_KVSF5"="lightsteelblue1", 
                  "Cuspidothrix_LMECYA_163"= "lightskyblue2", 
                  "Cyanothece_PCC-8801"= "lightskyblue",
                  "Planktothricoides_SR001" = "mediumturquoise",
                  "Richelia_HH01"="turquoise4", 
                  "Rivularia_PCC-7116"="dodgerblue4",
                  "Schizothrix_LEGE_07164" = "midnightblue",
                  
                  #potential toxin producers and nitrogen fixers
                  "Aphanizomenon_NIES81"="mediumorchid4",
                  "Gloeotrichia_SAG_32.84" = "mediumpurple4",
                  "Pseudanabaena_PCC-7429" = "purple4")

my_colours_seq<-c("Microcystis spp."="indianred4", 
                  "Phormidium spp."= "firebrick4"  , 
                  "Planktothrix spp."="darkred", 
                  
                  #potential nitrogen fixers
                  "Calothrix spp."="lightsteelblue1", 
                  "Cuspidothrix spp."= "lightskyblue2", 
                  "Cyanothece spp."= "lightskyblue",
                  "Planktothricoides spp." = "mediumturquoise",
                  "Richelia spp."="turquoise4", 
                  "Rivularia spp."="dodgerblue4",
                  "Schizothrix spp." = "midnightblue",
                  
                  #potential toxin producers and nitrogen fixers
                  "Aphanizomenon spp."="mediumorchid4",
                  "Gloeotrichia spp." = "mediumpurple4",
                  "Pseudanabaena spp." = "purple4",
                  
                  "Other Cyanobacteria"="grey")

x<-c("June 23rd",
     "July 20th",
     "Aug 3rd",
     "Aug 23rd",
     "Aug 31st",
     "Sept 27th")

t1 <- trans_abund$new(dataset = cyano, taxrank = "Genus") 

#specify order of bars (from top to bottom)
t1$data_abund$Taxonomy<-factor(t1$data_abund$Taxonomy,
                               levels=c("Microcystis_PCC-7914", "Phormidium_SAG_37.90" , "Planktothrix_NIVA-CYA_15", 
                                        
                                        "Calothrix_KVSF5",  "Cuspidothrix_LMECYA_163",  "Cyanothece_PCC-8801", "Planktothricoides_SR001",
                                        "Richelia_HH01", "Rivularia_PCC-7116","Schizothrix_LEGE_07164",
                                        
                                        "Aphanizomenon_NIES81", "Gloeotrichia_SAG_32.84","Pseudanabaena_PCC-7429"),
                               labels = c("Microcystis spp.", "Phormidium spp." , "Planktothrix spp.", 
                                          
                                          "Calothrix spp.",  "Cuspidothrix spp.",  "Cyanothece spp.", "Planktothricoides spp.",
                                          "Richelia spp.", "Rivularia spp.","Schizothrix spp.",
                                          
                                          "Aphanizomenon spp.", "Gloeotrichia spp.","Pseudanabaena spp."))
t1$data_abund$Taxonomy <- as.character(t1$data_abund$Taxonomy)
t1$data_abund$Taxonomy<-replace_na(t1$data_abund$Taxonomy, "Other Cyanobacteria")
t1$data_abund$Taxonomy <- factor(t1$data_abund$Taxonomy)

t1$data_abund$Taxonomy<-factor(t1$data_abund$Taxonomy,
                               levels=c("Other Cyanobacteria",
                                        "Microcystis spp.", "Phormidium spp." , "Planktothrix spp.", 
                                        
                                        "Calothrix spp.",  "Cuspidothrix spp.",  "Cyanothece spp.", "Planktothricoides spp.",
                                        "Richelia spp.", "Rivularia spp.","Schizothrix spp.",
                                        
                                        "Aphanizomenon spp.", "Gloeotrichia spp.","Pseudanabaena spp."))

write.csv(t1$data_abund, file = "C:/Users/kaitl/Documents/July_Analysis/merge_concat.genus.mc.nf.abund.csv", 
          row.names = TRUE)


#uses mean values, not accurate representation due to samples w/ no hits
ggplot(t1$data_abund, aes(x=factor(t1$data_abund$Date, level=x),
                          y=t1$data_abund$Abundance, fill = t1$data_abund$Taxonomy)) + 
  geom_bar(stat="summary", fun=mean) + 
  scale_fill_manual(values=my_colours_seq, name = "Genus") + 
  xlab("Field Day") +
  ylab("Relative Abundance") +
  my_theme_xlg

#try to assess relative abundance w/ mean values by removing the samples w/ no hits
t2<-clone(cyano)
t2$taxa_abund$Genus[t2$taxa_abund$Genus==0] <- NA
t2$taxa_abund<-drop_na(t2$taxa_abund, Genus)
t2 <- trans_abund$new(dataset = t2, taxrank = "Genus") 


ggplot(t2$data_abund, aes(x=factor(Date, level=x),
                          y=Abundance, fill = Taxonomy)) + 
  geom_bar(stat="summary", fun=mean) + 
  scale_fill_manual(values=my_colours_seq, name = "Genus") + 
  xlab("Field Day") +
  ylab("Relative Abundance") +
  my_theme_xlg

#shows each observation independently - probably best representation of the data
ggplot(t1$data_abund, aes(x=t1$data_abund$IDL,y=t1$data_abund$Abundance, 
                          fill = t1$data_abund$Taxonomy)) + 
  geom_bar(stat="sum", show.legend=c(size=FALSE)) + 
  scale_fill_manual(values=my_colours_seq, name = "Genus") + 
  xlab("SWMP") +
  ylab("Relative Abundance") +
  facet_wrap(factor(t1$data_abund$Date, level=x))+
  my_theme_xlg

#now lets look at abundance with all bacteria
t3<- clone(no_ctrl)
t3$tax_table <- subset(t3$tax_table, Kingdom == "k__Bacteria")
t3$tidy_dataset()
t3
#microtable-class object:
#sample_table have 69 rows and 44 columns
#otu_table have 8544 rows and 69 columns
#tax_table have 8544 rows and 7 columns
#phylo_tree have 8544 tips
#rep_fasta have 8544 sequences
t3$cal_betadiv()

t3 <- trans_abund$new(dataset = t3, taxrank = "Genus") 

t3$data_abund$Taxonomy<-factor(t3$data_abund$Taxonomy,
                               levels=c("Microcystis_PCC-7914", "Phormidium_SAG_37.90" , "Planktothrix_NIVA-CYA_15", 
                                        
                                        "Calothrix_KVSF5",  "Cuspidothrix_LMECYA_163",  "Cyanothece_PCC-8801", "Planktothricoides_SR001",
                                        "Richelia_HH01", "Rivularia_PCC-7116","Schizothrix_LEGE_07164",
                                        
                                        "Aphanizomenon_NIES81", "Gloeotrichia_SAG_32.84","Pseudanabaena_PCC-7429"),
                               labels = c("Microcystis spp.", "Phormidium spp." , "Planktothrix spp.", 
                                          
                                          "Calothrix spp.",  "Cuspidothrix spp.",  "Cyanothece spp.", "Planktothricoides spp.",
                                          "Richelia spp.", "Rivularia spp.","Schizothrix spp.",
                                          
                                          "Aphanizomenon spp.", "Gloeotrichia spp.","Pseudanabaena spp."))
t3$data_abund$Taxonomy <- as.character(t3$data_abund$Taxonomy)
t3$data_abund$Taxonomy<-replace_na(t3$data_abund$Taxonomy, "Other Bacteria")
t3$data_abund$Taxonomy <- factor(t3$data_abund$Taxonomy)

t3$data_abund$Taxonomy<-factor(t3$data_abund$Taxonomy,
                               levels=c("Other Bacteria",
                                        "Microcystis spp.", "Phormidium spp." , "Planktothrix spp.", 
                                        
                                        "Calothrix spp.",  "Cuspidothrix spp.",  "Cyanothece spp.", "Planktothricoides spp.",
                                        "Richelia spp.", "Rivularia spp.","Schizothrix spp.",
                                        
                                        "Aphanizomenon spp.", "Gloeotrichia spp.","Pseudanabaena spp."))

my_colours_seq2<-c("Microcystis spp."="indianred4", 
                   "Phormidium spp."= "firebrick4"  , 
                   "Planktothrix spp."="darkred", 
                   
                   #potential nitrogen fixers
                   "Calothrix spp."="lightsteelblue1", 
                   "Cuspidothrix spp."= "lightskyblue2", 
                   "Cyanothece spp."= "lightskyblue",
                   "Planktothricoides spp." = "mediumturquoise",
                   "Richelia spp."="turquoise4", 
                   "Rivularia spp."="dodgerblue4",
                   "Schizothrix spp." = "midnightblue",
                   
                   #potential toxin producers and nitrogen fixers
                   "Aphanizomenon spp."="mediumorchid4",
                   "Gloeotrichia spp." = "mediumpurple4",
                   "Pseudanabaena spp." = "purple4",
                   
                   "Other Bacteria"="grey")

ggplot(t3$data_abund, aes(x=IDL,y=Abundance, 
                          fill = Taxonomy)) + 
  geom_bar(stat="summary", fun=mean, show.legend=c(size=FALSE)) + 
  scale_fill_manual(values=my_colours_seq2, name = "Genus") + 
  xlab("SWMP") +
  ylab("Relative Abundance") +
  my_theme_xlg


ggplot(t3$data_abund, aes(x=IDL,y=Abundance, 
                          fill = Taxonomy)) + 
  geom_bar(stat="sum", show.legend=c(size=FALSE)) + 
  scale_fill_manual(values=my_colours_seq2, name = "Genus") + 
  xlab("SWMP") +
  ylab("Relative Abundance") +
  facet_wrap(factor(t3$data_abund$Date, level=x))+
  my_theme_xlg


ggplot(t3$data_abund, aes(x=IDL,y=Abundance, 
                          fill = Taxonomy)) + 
  geom_bar(stat="sum", show.legend=c(size=FALSE)) + 
  scale_fill_manual(values=my_colours_seq2, name = "Genus") + 
  xlab("SWMP") +
  ylab("Relative Abundance") +
  facet_wrap(factor(t3$data_abund$Date, level=x))+
  coord_cartesian(ylim=c(0,1.5))+
  my_theme_xlg




ggplot(t1$data_abund) + 
  geom_bar(data=t1$data_abund, aes(x=t1$data_abund$Field_day,
                                   y=t1$data_abund$Abundance, fill = t1$data_abund$Taxonomy),
           position='stack', stat='sum') + 
  geom_line(data=t1$data_abund, aes(x=t1$data_abund$Field_day,
                                    y=t1$data_abund$copies_mcyE_L), 
            linewidth=1.5, color="red")+
  scale_fill_manual(values=my_colours_seq, name = "Genus") + 
  xlab("Field Day") +
  ylab("Abundance") +
  my_theme_xlg


ggplot(t1$data_abund) + 
  geom_line(data=t1$data_abund, aes(x=t1$data_abund$Field_day,
                                    y=t1$data_abund$copies_mcyE_L), 
            linewidth=1.5, color="red")+
  scale_fill_manual(values=my_colours_seq, name = "Genus") + 
  xlab("Field Day") +
  ylab("Abundance") +
  my_theme_xlg


