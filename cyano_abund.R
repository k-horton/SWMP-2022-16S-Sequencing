library(decontam)
library(file2meco)
library(microeco)
library(dplyr)
library(rstatix)
library(ggplot2)
library(phyloseq)
library(EnvStats)
library(ggpubr)

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


ggsave("all_bact_abund_bar.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

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

ggsave("zoom_all_bact_abund_bar.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")


#### Check normality of the data #####
# Subset just cyanobacteria data
cyano<-subset(all_bact_abund$data_abund, all_bact_abund$data_abund$Taxonomy=="Cyanobacteria")


# first, aggregate abundance data 
# (multiple entries exist for each sample due to one entry being created 
#      for each taxa, even when it is absent)
# I already simplified to just Cyanobacteria overall, so we can condense these 
#    entries by summing Cyanobacteria abundance for each observation
cyano_sum<-aggregate(cyano$Abundance, by=list(Category=cyano$Sample), FUN=sum)

cyano_sum$SampleID<-cyano_sum$Category
cyano_sum$Abundance<-cyano_sum$x

cyano_sum<-subset(cyano_sum, select=-c(Category, x))

# Separate sample metadata so it can be joined to cyanobacteria abundance data
metadata<-no_ctrl$sample_table
# convert rownames to first column
metadata<-tibble::rownames_to_column(metadata, "SampleID")

# merge both dataframes
cyano_merge<-merge(cyano_sum, metadata, by="SampleID")

# check normality
res_aov <- aov(Abundance ~ Date,
               data = cyano_merge)
par(mfrow = c(1, 2)) 
# histogram
hist(res_aov$residuals)
# QQ-plot
qqPlot(res_aov$residuals, add.line = TRUE)
# not completely normal, more normal than data that was not summed

# try log transformation

res_aov <- aov(log10(Abundance+0.1) ~ Date,
               data = cyano_merge)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# This is much more normal - use log transformation moving forward
#      for parametric tests


# Now check normality for all other variables (already determined
#      previously for most water quality variables, see repository
#      SWMP-2022-Water-Quality https://github.com/k-horton/SWMP-2022-Water-Quality)

# since we don't need abundance data for this, we can just load in our basic dataset
#   without the abundance data
library(readxl)

swmp<-read_excel(paste0(dir,"/metadata_decontam.xlsx"))
#remove control samples
swmp<-subset(swmp, Control == FALSE)

# Chlorophyll-a #
res_aov <- aov(Chla_gL ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

  # try log transformation
res_aov <- aov(log10(Chla_gL+0.01) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

  # This is better - use log transformation 

# Unicellular cyanobacteria cells #
res_aov <- aov(Unicellularcells_L ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

  # try log transformation
res_aov <- aov(log10(Unicellularcells_L+0.1) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
swmp_sub<-swmp
swmp_sub$Unicellularcells_L[swmp_sub$Unicellularcells_L==0] <- NA
swmp_sub<-swmp_sub[complete.cases(swmp_sub),]

res_aov <- aov(Unicellularcells_L ~ Date,
               data = swmp_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Unicellularcells_L+0.01) ~ Date,
               data = swmp_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

#this looks good! Use log transformation


# Colonial cyanobacteria cells #
res_aov <- aov(Colonialcells_L ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

  # try log transformation
res_aov <- aov(log10(Colonialcells_L+0.01) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
  # This is better - use log transformation 


# Filamentous cyanobacteria cells #
res_aov <- aov(Filamentouscells_L ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Filamentouscells_L+0.01) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# This is questionable, try some other transformations:

  # root transformation -> x^(1/2)
res_aov <- aov((Filamentouscells_L+0.01)^(1/2) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

  # reciprocal transformation -> 1/x
res_aov <- aov(1/(Filamentouscells_L+0.01) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

  # squared transformation -> x^2
res_aov <- aov((Filamentouscells_L+0.01)^(2) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
swmp_sub<-swmp
swmp_sub$Filamentouscells_L[swmp_sub$Filamentouscells_L==0] <- NA
swmp_sub<-swmp_sub[complete.cases(swmp_sub),]

res_aov <- aov(Filamentouscells_L ~ Date,
               data = swmp_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

  # try log transformation
res_aov <- aov(log10(Filamentouscells_L+0.01) ~ Date,
               data = swmp_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

  #this looks good! Use log transformation

# Total cyanobacteria cells #
res_aov <- aov(Totalcells_L ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(Totalcells_L+0.01) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)
# This is better - use log transformation 

# 16S copies / L #
res_aov <- aov(copies_16S_L ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# try log transformation
res_aov <- aov(log10(copies_16S_L+0.01) ~ Date,
               data = swmp)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

# drop zero values and re-check
swmp_sub<-swmp
swmp_sub$copies_16S_L[swmp_sub$copies_16S_L==0] <- NA
swmp_sub<-swmp_sub[complete.cases(swmp_sub),]

res_aov <- aov(copies_16S_L ~ Date,
               data = swmp_sub)
par(mfrow = c(1, 2)) 
hist(res_aov$residuals)
qqPlot(res_aov$residuals,   add.line = TRUE)

res_aov <- aov(log10(copies_16S_L+0.01) ~ Date,
               data = swmp_sub)
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

##### Cyano Abundance - Repeated-measures ANOVA#####
# set date as a factor and order chronologically
cyano_merge$Date<-factor(cyano_merge$Date,
                         levels=c("June 23rd",
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
                                    "Sept 27th")
                         )


# Abundance summary statistics by field day - remember this is relative abundance of 
#   cyanobacteria compared to all bacteria detected with 16S sequencing
cyano_merge %>%
  group_by(Date) %>%
  get_summary_stats(Abundance, type = "mean_sd")
#   Date       variable      n  mean    sd
#    <chr>     <fct>     <dbl> <dbl> <dbl>
# 1 June 23rd   Abundance    12 2.03  3.72 
# 2 July 20th   Abundance    12 6.07  6.37   <- second highest abundance
# 3 Aug 3rd     Abundance    11 6.18  9.45   <- highest abundance
# 4 Aug 23rd    Abundance    11 1.07  0.923
# 5 Aug 31st    Abundance    12 0.236 0.472  <-lowest abundance
# 6 Sept 27th   Abundance    11 0.578 1.92  

# summary statistics by SWMP
cyano_merge %>%
  group_by(IDL) %>%
  get_summary_stats(Abundance, type = "mean_sd")
#  IDL   variable      n  mean     sd
#  <chr> <fct>     <dbl> <dbl>  <dbl>
# 1 A     Abundance     6 3.6    7.72 
# 2 B     Abundance     6 1.59   3.38 
# 3 C     Abundance     6 0.117  0.287 <- Lowest abundance
# 4 D     Abundance     6 1.85   2.92 
# 5 E     Abundance     5 4.29   6.81 
# 6 F     Abundance     6 2.34   4.73 
# 7 G     Abundance     5 1.00   1.34 
# 8 H     Abundance     5 3.44   6.54 
# 9 I     Abundance     6 9.84   10.8  <- Highest abundance
# 10 J    Abundance     6 1.07   1.36 
# 11 K    Abundance     6 1.06   0.992
# 12 L    Abundance     6 2.30   3.77 

# create new column w/ log transformed abundance data
cyano_merge$logAbund<-log10(cyano_merge$Abundance + 0.1)

# Repeated measures ANOVA (Within-subject):
res.aov <- anova_test(data = cyano_merge, dv = logAbund, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)
#ANOVA Table (type III tests)
#   Effect   DFn DFd     F     p      p<.05     ges
#1    Date   5  40  8.639 1.28e-05    *     0.403

# Abundance of cyanobacteria is significantly different between field days

# determine where the significant differences are:
t.test<-pairwise.t.test(cyano_merge$logAbund, cyano_merge$Date, 
                p.adjust.method = "bonferroni")
t.test
#Pairwise comparisons using t tests with pooled SD 

#data:  cyano_merge$logAbund and cyano_merge$Date 

#             June 23rd   July 20th   Aug 3rd   Aug 23rd  Aug 31st   
#  July 20th  0.0534     -           -         -         -     
#  Aug 3rd    0.5958     1.0000      -         -         -     
#  Aug 23rd   1.0000     0.4621      1.0000    -         -     
#  Aug 31st   0.8811     8.6e-05     0.0027    0.1584    -     
#  Sept 27th  0.4616     3.8e-05     0.0012    0.0769    1.0000

# P value adjustment method: bonferroni 
# Bonferroni correction multiplies p-values by the number of comparisons, this is the most 
#  conservative option


# create a boxplot to show this data
# Define theme for boxplot
box_theme <- theme(strip.background = element_rect(fill="white", color="black"),
                   panel.background = element_rect(fill="white", color = "black"),
                   panel.margin = unit(0.3, "lines"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.margin=margin(0.5,0.5,0.5,0.5,unit="cm"),
                   axis.title.x = element_text(size = 25, face = "bold"),
                   axis.title.y = element_text(size = 25, face = "bold"),
                   axis.text.x = element_text(size = 20, vjust = 0.5),
                   axis.text.y = element_text(size = 20),
                   legend.title = element_text(size = 22, face = "bold"),
                   legend.text = element_text(size = 20)
                   )


# Boxplot total filamentous cyanobacterial cells across field days
cyano_box<-ggplot(cyano_merge, aes(x=Date,
                             y=Abundance)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Relative Abundance of Cyanobacteria (%)") +
  box_theme
cyano_box


# add significance to plot


group_cyano<-cyano_merge %>%
  group_by(Date)%>%
  mutate(max_value = max(Abundance))%>%
  mutate(min_value = min(Abundance)) %>%
  mutate(abund.mean = mean(Abundance))

#subset just min max data
min_max<-unique(group_cyano[, c("Date", "max_value", "min_value", "abund.mean")]
                )

# generate p-values for labels
t.test.label<-cyano_merge %>%
  pairwise_t_test(logAbund~Date, p.adjust.method = "bonferroni") 
head(t.test.label,5)

# add matching grouping columns to min_max dataframe
min_max$group1<-min_max$Date
min_max$group2<-min_max$Date

# add group 1 min max values
by<-join_by(group1)
t.test.merge<-left_join(t.test.label, 
                        min_max[,c("group1", "min_value", "max_value", "abund.mean")], 
                        by)

#add group 1 prefix
colnames(t.test.merge)[c(11:13)] <- paste("g1",  colnames(t.test.merge)[c(11:13)], 
                                                 sep = '.')

# add group 2 min max values
by<-join_by(group2)
t.test.merge<-left_join(t.test.merge, 
                        min_max[,c("group2", "min_value", "max_value", "abund.mean")], 
                        by)

#add group 2 prefix
colnames(t.test.merge)[c(14:16)] <- paste("g2",  colnames(t.test.merge)[c(14:16)], 
                                          sep = '.')

# add unique column id w/ both dates
t.test.merge$comp<-paste(t.test.merge$group1, sep="-", t.test.merge$group2)


# define overall max, and max mean between dates
t.test.merge <- t.test.merge %>%
  mutate(over_max = pmax(g1.max_value, g2.max_value,  na.rm = TRUE))%>%
  mutate(max_mean = pmax(g1.abund.mean, g2.abund.mean,  na.rm = TRUE)
         )

#subset comparisons we want to include (p < 0.05)
label_sub<-subset(t.test.merge, c(comp == "June 23rd-July 20th" |
                                    comp =="July 20th-Aug 31st"|
                                    comp =="July 20th-Sept 27th"|
                                    comp =="Aug 3rd-Aug 31st"|
                                    comp =="Aug 3rd-Sept 27th")
                  ) 

# add full label column
label_sub$label<-ifelse((label_sub$p.adj > 0.05),
                       paste("p =", label_sub$p.adj, "(NS)"),
                       ifelse((label_sub$p.adj < 0.001),
                              paste("p =", label_sub$p.adj, "***"),
                              ifelse((label_sub$p.adj< 0.01),
                                     paste("p =", label_sub$p.adj, "**"),
                                     ifelse((label_sub$p.adj < 0.05),
                                            paste("p =", label_sub$p.adj, "*"),
                                            paste("error")
                                            )
                                     )
                              )
                       )
                       
# plot with significance
cyano_box+
  stat_pvalue_manual(label_sub, label = "label", y.position = 20, step.increase=0.1)

# Customize y position for each comparison
y <- c(22.5, 24, 25.5, 21, 19.5) 
cyano_box+
  stat_pvalue_manual(label_sub, label = "label", y.position=y)

ggsave("cyano_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")
##### Total cyano cells Abundance - Repeated-measures ANOVA#####
t_cyano<-swmp

# set date as a factor and order chronologically
t_cyano$Date<-factor(t_cyano$Date,
                   levels=c("June 23rd","July 20th","Aug 3rd",
                            "Aug 23rd", "Aug 31st",  "Sept 27th"),
                   labels = c("June 23rd", "July 20th", "Aug 3rd",
                              "Aug 23rd","Aug 31st",  "Sept 27th")
)

# Total cells summary statistics by field day
t_cyano %>%
  group_by(Date) %>%
  get_summary_stats(Totalcells_L, type = "mean_sd")
#   Date      variable         n       mean         sd
#   <fct>     <fct>        <dbl>      <dbl>      <dbl>
# 1 June 23rd Totalcells_L    12 115440374. 137931107.
# 2 July 20th Totalcells_L    12  78041504.  45295198.
# 3 Aug 3rd   Totalcells_L    12 103781047.  67129278.
# 4 Aug 23rd  Totalcells_L    12 146064602. 127764901. <-max
# 5 Aug 31st  Totalcells_L    12  77338927.  55161716. <-min
# 6 Sept 27th Totalcells_L    11 142148172. 101659588.

# summary statistics by SWMP
t_cyano %>%
  group_by(IDL) %>%
  get_summary_stats(Totalcells_L, type = "mean_sd")
#   IDL   variable         n       mean         sd
#   <chr> <fct>        <dbl>      <dbl>      <dbl>
# 1 A     Totalcells_L     6 136701598. 190961026. <- max
# 2 B     Totalcells_L     6  81476587.  46119442.
# 3 C     Totalcells_L     6 119359512.  37879971.
# 4 D     Totalcells_L     6 186862171. 117123086.
# 5 E     Totalcells_L     5 100952569.  82060920.
# 6 F     Totalcells_L     6  75362820.  46087871.
# 7 G     Totalcells_L     6 194436029. 130561297.  
# 8 H     Totalcells_L     6  72695230.  77612654.
# 9 I     Totalcells_L     6 109679626.  37579628.
# 10 J    Totalcells_L     6  46613235.  19613248.
# 11 K    Totalcells_L     6  38518157.  20458315.  <- min
# 12 L    Totalcells_L     6 156105783.  96485177.

# create new column w/ log transformed total cyanobcteria cell data
t_cyano$logT_cyano<-log10(t_cyano$Totalcells_L+0.01)

# Repeated measures ANOVA (Within-subject):
res.aov <- anova_test(data = t_cyano, dv = logT_cyano, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)
#ANOVA Table (type III tests)

#   Effect DFn DFd     F     p p<.05   ges
#1   Date   5  50 1.565 0.187       0.091

# determine where the significant differences are:
t.test<-pairwise.t.test(t_cyano$logT_cyano, t_cyano$Date, 
                        p.adjust.method = "bonferroni")
t.test
#Pairwise comparisons using t tests with pooled SD 
# data:  t_cyano$logT_cyano and t_cyano$Date 


#           June 23rd July 20th Aug 3rd Aug 23rd Aug 31st
#  July 20th 1.00      -         -       -        -       
#  Aug 3rd   1.00      1.00      -       -        -       
#  Aug 23rd  1.00      1.00      1.00    -        -       
#  Aug 31st  1.00      1.00      1.00    0.81     -       
#  Sept 27th 1.00      1.00      1.00    1.00     1.00    

# P value adjustment method: bonferroni 

# create a boxplot to show this data
# Define theme for boxplot
box_theme <- theme(strip.background = element_rect(fill="white", color="black"),
                   panel.background = element_rect(fill="white", color = "black"),
                   panel.margin = unit(0.3, "lines"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.margin=margin(0.5,0.5,0.5,0.5,unit="cm"),
                   axis.title.x = element_text(size = 25, face = "bold"),
                   axis.title.y = element_text(size = 25, face = "bold"),
                   axis.text.x = element_text(size = 20, vjust = 0.5),
                   axis.text.y = element_text(size = 20),
                   legend.title = element_text(size = 22, face = "bold"),
                   legend.text = element_text(size = 20)
)


# Boxplot total cyanobacterial cells across field days
t_cyano_box<-ggplot(t_cyano, aes(x=Date,
                             y=Totalcells_L)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Total cyanobacteria cells / L") +
  box_theme
t_cyano_box

# No significant comparisons, do not need to add significance to plot

ggsave("t_cyano_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

##### qPCR 16S copies - Repeated-measures ANOVA#####
q16S<-swmp

# set date as a factor and order chronologically
q16S$Date<-factor(q16S$Date,
                     levels=c("June 23rd","July 20th","Aug 3rd",
                              "Aug 23rd", "Aug 31st",  "Sept 27th"),
                     labels = c("June 23rd", "July 20th", "Aug 3rd",
                                "Aug 23rd","Aug 31st",  "Sept 27th")
)

# Total cells summary statistics by field day
q16S %>%
  group_by(Date) %>%
  get_summary_stats(copies_16S_L, type = "mean_sd")
#   Date      variable         n        mean          sd
#   <fct>     <fct>        <dbl>       <dbl>       <dbl>
# 1 June 23rd copies_16S_L    11 1035480626. 1417774672.
# 2 July 20th copies_16S_L    12 6277918231. 5719850288. <- max
# 3 Aug 3rd   copies_16S_L    12 2282486071. 3020406547.
# 4 Aug 23rd  copies_16S_L    12  786531677.  644141534.
# 5 Aug 31st  copies_16S_L    12  603297800.  734441378.  <- min
# 6 Sept 27th copies_16S_L    10 1724015839. 1485525455.

# summary statistics by SWMP
q16S %>%
  group_by(IDL) %>%
  get_summary_stats(copies_16S_L, type = "mean_sd")
#   IDL   variable      n        mean          sd
# <chr> <fct>        <dbl>       <dbl>       <dbl>
# 1 A     copies_16S_L     5 1918502569. 1824001268.
# 2 B     copies_16S_L     6 1120276430. 1723042676.
# 3 C     copies_16S_L     6 1280231716. 1897867023.
# 4 D     copies_16S_L     6 2446178230. 3168222608.
# 5 E     copies_16S_L     4  180147509.  161281756. <- min
# 6 F     copies_16S_L     6 1355851183.  820341352.
# 7 G     copies_16S_L     6 4605813454. 7674876302.
# 8 H     copies_16S_L     6 1772995405. 3912588206.
# 9 I     copies_16S_L     6 4606169650. 4613503945. <- max
# 10 J    copies_16S_L     6 2462426485. 3296992451.
# 11 K    copies_16S_L     6 2323092586  2380461313.
# 12 L    copies_16S_L     6  980322816. 1319581171.

# create new column w/ log transformed total cyanobcteria cell data
q16S$log16S<-log10(q16S$copies_16S_L+0.01)

q16S_drop<-drop_na(q16S, copies_16S_L)
# Repeated measures ANOVA (Within-subject):
res.aov <- anova_test(data = q16S_drop, dv = log16S, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)
#ANOVA Table (type III tests)

#   Effect  DFn   DFd     F     p p<.05   ges
#1   Date   2.02 18.19  1.294   0.298       0.109

# determine where the significant differences are:
t.test<-pairwise.t.test(q16S$log16S, q16S$Date, 
                        p.adjust.method = "bonferroni")
t.test
#Pairwise comparisons using t tests with pooled SD 
# data:  q16S$log16S and q16S$Date 


#           June 23rd July 20th Aug 3rd Aug 23rd Aug 31st
#  July 20th 1.00      -         -       -        -       
#  Aug 3rd   1.00      1.00      -       -        -       
#  Aug 23rd  1.00      1.00      1.00    -        -       
#  Aug 31st  1.00      0.62      1.00    1.00     -       
#  Sept 27th 1.00      1.00      1.00    1.00     1.00  
#
# P value adjustment method: bonferroni 

# create a boxplot to show this data
# Define theme for boxplot
box_theme <- theme(strip.background = element_rect(fill="white", color="black"),
                   panel.background = element_rect(fill="white", color = "black"),
                   panel.margin = unit(0.3, "lines"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.margin=margin(0.5,0.5,0.5,0.5,unit="cm"),
                   axis.title.x = element_text(size = 25, face = "bold"),
                   axis.title.y = element_text(size = 25, face = "bold"),
                   axis.text.x = element_text(size = 20, vjust = 0.5),
                   axis.text.y = element_text(size = 20),
                   legend.title = element_text(size = 22, face = "bold"),
                   legend.text = element_text(size = 20)
)


# Boxplot total cyanobacterial cells across field days
q16S_box<-ggplot(q16S, aes(x=Date,
                                 y=copies_16S_L)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("16S copies / L") +
  box_theme
q16S_box

# No significant comparisons, do not need to add significance to plot

ggsave("q16S_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

##### Chlorophyll a - Repeated-measures ANOVA#####
chla<-swmp

# set date as a factor and order chronologically
chla$Date<-factor(chla$Date,
                     levels=c("June 23rd","July 20th","Aug 3rd",
                              "Aug 23rd", "Aug 31st",  "Sept 27th"),
                     labels = c("June 23rd", "July 20th", "Aug 3rd",
                                "Aug 23rd","Aug 31st",  "Sept 27th")
)

# Total cells summary statistics by field day
chla %>%
  group_by(Date) %>%
  get_summary_stats(Chla_gL, type = "mean_sd")
#   Date      variable     n  mean    sd
#   <fct>     <fct>    <dbl> <dbl> <dbl>
# 1 June 23rd Chla_gL     12 0.015 0.015
# 2 July 20th Chla_gL     12 0.002 0.003  <- min
# 3 Aug 3rd   Chla_gL     12 0.026 0.044  <- max
# 4 Aug 23rd  Chla_gL     12 0.016 0.024
# 5 Aug 31st  Chla_gL     12 0.01  0.012
# 6 Sept 27th Chla_gL     11 0.019 0.016

# summary statistics by SWMP
chla %>%
  group_by(IDL) %>%
  get_summary_stats(Chla_gL, type = "mean_sd")
#   IDL   variable     n  mean    sd
# <chr> <fct>    <dbl> <dbl> <dbl>
# 1 A     Chla_gL      6 0.009 0.006
# 2 B     Chla_gL      6 0     0.001  <- min
# 3 C     Chla_gL      6 0.04  0.056  <- max
# 4 D     Chla_gL      6 0.036 0.021
# 5 E     Chla_gL      5 0.003 0.003
# 6 F     Chla_gL      6 0.001 0.001
# 7 G     Chla_gL      6 0.027 0.015
# 8 H     Chla_gL      6 0.008 0.007
# 9 I     Chla_gL      6 0.03  0.028
# 10 J    Chla_gL      6 0.005 0.006
# 11 K    Chla_gL      6 0.009 0.011
# 12 L    Chla_gL      6 0.004 0.006

# create new column w/ log transformed total cyanobcteria cell data
chla$logChla<-log10(chla$Chla_gL+0.01)

# Repeated measures ANOVA (Within-subject):
res.aov <- anova_test(data = chla, dv = logChla, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)
# ANOVA Table (type III tests)

# Effect DFn DFd     F     p p<.05   ges
# 1   Date   5  50 4.612 0.002     * 0.154

# determine where the significant differences are:
t.test<-pairwise.t.test(chla$logChla, chla$Date, 
                        p.adjust.method = "bonferroni")
t.test
#Pairwise comparisons using t tests with pooled SD 
# data:  chla$logChla and chla$Date 

#          June 23rd July 20th Aug 3rd Aug 23rd Aug 31st
# July 20th 0.314     -         -       -        -       
# Aug 3rd   1.000     0.234     -       -        -       
# Aug 23rd  1.000     0.521     1.000   -        -       
# Aug 31st  1.000     1.000     1.000   1.000    -       
# Sept 27th 1.000     0.064     1.000   1.000    1.000   

# P value adjustment method: bonferroni 

# create a boxplot to show this data
# Define theme for boxplot
box_theme <- theme(strip.background = element_rect(fill="white", color="black"),
                   panel.background = element_rect(fill="white", color = "black"),
                   panel.margin = unit(0.3, "lines"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.margin=margin(0.5,0.5,0.5,0.5,unit="cm"),
                   axis.title.x = element_text(size = 25, face = "bold"),
                   axis.title.y = element_text(size = 25, face = "bold"),
                   axis.text.x = element_text(size = 20, vjust = 0.5),
                   axis.text.y = element_text(size = 20),
                   legend.title = element_text(size = 22, face = "bold"),
                   legend.text = element_text(size = 20)
)


# Boxplot total cyanobacterial cells across field days
chla_box<-ggplot(chla, aes(x=Date,
                                 y=Chla_gL)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Chlorophyll a (g/L)") +
  box_theme
chla_box

# No significant comparisons, do not need to add significance to plot

ggsave("chla_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

