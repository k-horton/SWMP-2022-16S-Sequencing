# Assess cyanobacterial community - specific taxa

# First, load sequence data and remove contaminant samples as determined previously 
# (see files decontam_prev.R, decontam_freq.R)
# If you already have the ps.prev.noncontam object from previous scripts, skip next 
# section and continue to "If you already have ps.prev.noncontam object"
library(decontam)
library(file2meco)
library(microeco)
library(dplyr)
library(rstatix)
library(ggplot2)
library(phyloseq)
library(EnvStats)

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



#### Continue on here: subset to only cyanobacteria #####

#subset taxa
cyano <- clone(no_ctrl)
cyano$tax_table <- subset(cyano$tax_table, Class == "c__Cyanobacteriia")

# remove non-cyano associated data 
cyano$tidy_dataset()
# 23 samples with 0 abundance are removed from the otu_table ...

# How many observations are remaining?
cyano
#microtable-class object:
#sample_table have 46 rows and 40 columns
#otu_table have 128 rows and 46 columns
#tax_table have 128 rows and 7 columns
#phylo_tree have 128 tips
#rep_fasta have 128 sequences

# There are 128 unique cyanobacteria OTU in the sample

# Which samples have these taxa?
cyano$sample_names()

#[1] "F1_SWMP1_S1"   "F1_SWMP4_S4"   "F1_SWMP6_S6"   "F1_SWMP7_S7"  
#[5] "F1_SWMP8_S8"   "F1_SWMP9_S9"   "F1_SWMP10_S10" "F1_SWMP11_S11"
#[9] "F1_SWMP12_S12" "F2_SWMP1_S13"  "F2_SWMP2_S14"  "F2_SWMP3_S15" 
#[13] "F2_SWMP4_S16"  "F2_SWMP5_S17"  "F2_SWMP6_S18"  "F2_SWMP7_S19" 
#[17] "F2_SWMP8_S20"  "F2_SWMP9_S21"  "F2_SWMP10_S22" "F2_SWMP11_S23"
#[21] "F2_SWMP12_S24" "F3_SWMP1_S25"  "F3_SWMP2_S26"  "F3_SWMP4_S28" 
#[25] "F3_SWMP5_S29"  "F3_SWMP6_S30"  "F3_SWMP8_S32"  "F3_SWMP9_S33" 
#[29] "F3_SWMP10_S34" "F3_SWMP11_S35" "F3_SWMP12_S36" "F4_SWMP1_S37" 
#[33] "F4_SWMP2_S38"  "F4_SWMP4_S40"  "F4_SWMP5_S41"  "F4_SWMP6_S42" 
#[37] "F4_SWMP7_S43"  "F4_SWMP9_S45"  "F4_SWMP10_S46" "F4_SWMP11_S47"
#[41] "F4_SWMP12_S48" "F5_SWMP4_S52"  "F5_SWMP7_S55"  "F5_SWMP10_S58"
#[45] "F5_SWMP12_S60" "F6_SWMP9_S68" 

#### Cyanobacteria family level taxonomy plots #####

# Plot boxplot of cyanobacteria taxa at the Family level
t1 <- trans_abund$new(dataset = cyano, taxrank = "Family", ntaxa = 10) #10 taxa present
t1$plot_box(group = "Field_day", xtext_angle = 30)

# Export .png file of relative cyanobacteria abundance boxplot at family level
ggsave("family_boxplot.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

# Export .csv file of relative cyanobacteria abundance data at the family level
write.csv(t1$data_abund, file = paste0(dir, "/Abundance_plots/family_boxplot.csv"), 
          row.names = TRUE)

# Bar-chart of relative abundance for each pond, grouped by field day
t1$plot_bar(others_color = "grey70", facet = "Field_day", xtext_keep = FALSE, legend_text_italic = FALSE)

# Bar chart of mean relative abundance for each field day
t1 <- trans_abund$new(dataset = cyano, taxrank = "Family", ntaxa = 10, groupmean = "Field_day")
t1$plot_bar(others_color = "grey70", xtext_keep = TRUE, legend_text_italic = FALSE, barwidth = 0.9)

# Export .png file of mean relative cyanobacteria abundance bar-chart at family level
ggsave("mean_family_barchart.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

#### Cyanobacteria genus level taxonomy plots #####

# Plot boxplot of cyanobacteria taxa at the Genus level
t1 <- trans_abund$new(dataset = cyano, taxrank = "Genus", ntaxa = 32) #32 taxa present
t1$plot_box(group = "Field_day", xtext_angle = 30)

#we'll save the .csv file with all taxa at genus level and then make a boxplot with only the top 15 taxa
write.csv(t1$data_abund, file = paste0(dir, "/Abundance_plots/genus_boxplot.csv"), 
          row.names = TRUE)

# Cut down to top 15 taxa (genus level)
t1 <- trans_abund$new(dataset = cyano, taxrank = "Genus", ntaxa = 15)
t1$plot_box(group = "Field_day", xtext_angle = 35)

# Export .png file of relative cyanobacteria abundance boxplot at genus level
ggsave("genus_boxplot.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

# Bar-chart of relative abundance for each pond, grouped by field day
t1$plot_bar(others_color = "grey70", facet = "Field_day", xtext_keep = FALSE, legend_text_italic = FALSE)

# Bar chart of mean relative abundance for each field day
t1 <- trans_abund$new(dataset = cyano, taxrank = "Genus", groupmean = "Field_day", ntaxa=32)
t1$plot_bar(others_color = "grey70", xtext_keep = TRUE, legend_text_italic = FALSE, barwidth = 0.9)

# Export .png file of mean relative cyanobacteria abundance bar-chart at genus level
ggsave("mean_genus_barchart.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")
#### Cyanobacteria diversity metrics #####

##### Filamentous cells Abundance - Repeated-measures ANOVA#####
filam<-swmp
#drop zeroes
filam_drop<-swmp
filam_drop$Filamentouscells_L[filam_drop$Filamentouscells_L==0] <- NA
filam_drop<-filam_drop[complete.cases(filam_drop),]

# set date as a factor and order chronologically
filam$Date<-factor(filam$Date,
                   levels=c("June 23rd","July 20th","Aug 3rd",
                            "Aug 23rd", "Aug 31st",  "Sept 27th"),
                   labels = c("June 23rd", "July 20th", "Aug 3rd",
                              "Aug 23rd","Aug 31st",  "Sept 27th")
)

filam_drop$Date<-factor(filam_drop$Date,
                        levels=c("June 23rd","July 20th", "Aug 3rd",
                                 "Aug 23rd","Aug 31st", "Sept 27th"),
                        labels = c("June 23rd", "July 20th", "Aug 3rd",
                                   "Aug 23rd","Aug 31st", "Sept 27th")
)

# Filamentous cells summary statistics by field day
filam %>%
  group_by(Date) %>%
  get_summary_stats(Filamentouscells_L, type = "mean_sd")
#   Date      variable               n     mean       sd
#   <fct>     <fct>              <dbl>    <dbl>    <dbl>
# 1 June 23rd Filamentouscells_L    12 1104685. 2381645.
# 2 July 20th Filamentouscells_L    12  488611. 1211249.
# 3 Aug 3rd   Filamentouscells_L    12  335243.  786422.<- min
# 4 Aug 23rd  Filamentouscells_L    12  607328.  969738.
# 5 Aug 31st  Filamentouscells_L    12 1409259. 3040185.
# 6 Sept 27th Filamentouscells_L    11 2718912. 6823813.<-max

filam_drop %>%
  group_by(Date) %>%
  get_summary_stats(Filamentouscells_L, type = "mean_sd")
#   Date      variable               n     mean       sd
#  <fct>     <fct>              <dbl>    <dbl>    <dbl>
# 1 June 23rd Filamentouscells_L     3  290099.  133496. <- min
# 2 July 20th Filamentouscells_L     5 1172665. 1741230.
# 3 Aug 3rd   Filamentouscells_L     3 1340971. 1174116.
# 4 Aug 23rd  Filamentouscells_L     5 1425900. 1069980.
# 5 Aug 31st  Filamentouscells_L     6 2818519. 3945576.
# 6 Sept 27th Filamentouscells_L     7 4253176. 8369753.  <- max

# summary statistics by SWMP
filam %>%
  group_by(IDL) %>%
  get_summary_stats(Filamentouscells_L, type = "mean_sd")
# IDL   variable               n     mean       sd
# <chr> <fct>              <dbl>    <dbl>    <dbl>
# 1 A     Filamentouscells_L     6   63858.  105030.
# 2 B     Filamentouscells_L     6   73289.  179520.
# 3 C     Filamentouscells_L     6 1419325. 2081588.
# 4 D     Filamentouscells_L     6  169858.  303144.
# 5 E     Filamentouscells_L     5   40740.   91097. <- min
# 6 F     Filamentouscells_L     6  841205. 1983920.
# 7 G     Filamentouscells_L     6 1237289. 2378532.
# 8 H     Filamentouscells_L     6  292347.  545572.
# 9 I     Filamentouscells_L     6 7488948. 8318941. <- max
# 10 J    Filamentouscells_L     6  468694.  441283.
# 11 K    Filamentouscells_L     6  623595. 1062153.
# 12 L    Filamentouscells_L     6  162564.  274041.

filam_drop %>%
  group_by(IDL) %>%
  get_summary_stats(Filamentouscells_L, type = "mean_sd")
#    IDL   variable               n     mean       sd
#<chr> <fct>              <dbl>    <dbl>    <dbl>
# 1 A     Filamentouscells_L     1  247350.      NA   <- not enough observations
# 2 B     Filamentouscells_L     1  439733.      NA   <- not enough observations
# 3 C     Filamentouscells_L     3 2785840. 2285185.
# 4 D     Filamentouscells_L     2  509573.  336545.
# 5 E     Filamentouscells_L     1  203700.      NA   <- min, not enough observations
# 6 F     Filamentouscells_L     2 2523616. 3344874.
# 7 G     Filamentouscells_L     2  690317.   16004.
# 8 H     Filamentouscells_L     2  877042.  680178.  <- max
# 9 I     Filamentouscells_L     5 7718174. 9279650.
# 10 J    Filamentouscells_L     4  703041.  323848.
# 11 K    Filamentouscells_L     3 1247191. 1286005.
# 12 L    Filamentouscells_L     3  325127.  329340.

# create new column w/ log transformed filamentous cell data
filam$logFilam<-log10(filam$Filamentouscells_L+0.01)

filam_drop$logFilam<-log10(filam_drop$Filamentouscells_L+0.01)

# Repeated measures ANOVA (Within-subject):
res.aov <- anova_test(data = filam, dv = logFilam, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)
#ANOVA Table (type III tests)

#   Effect DFn DFd     F     p p<.05   ges
# 1   Date   5  50 1.168 0.338       0.081

# the filam_drop dataframe is missing some SWMP observations, 
#    so cannot do a within-subject analysis

t.test<-pairwise.t.test(filam$logFilam, filam$Date, 
                        p.adjust.method = "bonferroni")
t.test

#Pairwise comparisons using t tests with pooled SD 

#data:  filam$logFilam and filam$Date 

#             June 23rd July 20th Aug 3rd Aug 23rd Aug 31st
#  July 20th 1.00      -         -       -        -       
#  Aug 3rd   1.00      1.00      -       -        -       
#  Aug 23rd  1.00      1.00      1.00    -        -       
#  Aug 31st  1.00      1.00      1.00    1.00     -       
#  Sept 27th 1.00      1.00      0.46    1.00     1.00    

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


# Boxplot total filamentous cyanobacterial cells across field days
filam_box<-ggplot(filam, aes(x=Date,
                             y=Filamentouscells_L)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Filamentous cyanobacteria cells / L") +
  box_theme
filam_box


# No significant comparisons, do not need to add significance to plot

ggsave("filam_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

##### Filamentous cells - Rosner Outlier test #####
ros_fil<-filam

#plot filamentous cells without any grouping
ros_filam_box<-ggplot(filam, aes(y=Filamentouscells_L)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("") +
  ylab("Filamentous cyanobacteria cells / L") +
  box_theme
ros_filam_box
# There looks to be 1, 2, 4, 5, or 7 outliers based on the boxplot

ros.test<-rosnerTest(ros_fil$Filamentouscells_L, k = 12, alpha = 0.0001, warn = TRUE)
# k = number of suspect outliers
ros.test$all.stats

# make list of outliers
outliers<-subset(ros.test$all.stats, Outlier == TRUE)
outliers<-outliers$Obs.Num

# remove outliers from dataframe
ros_fil<-ros_fil%>%
  slice(-c(outliers))

#check statistics
res.aov <- anova_test(data = ros_fil, dv = logFilam, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)

t.test<-pairwise.t.test(ros_fil$logFilam, ros_fil$Date, 
                        p.adjust.method = "bonferroni")
t.test
# there is not much change in the reusults even with the outliers removed,
#     still no signicant differences in filamentous cells between field days

##### Colonial cells Abundance - Repeated-measures ANOVA#####
colony<-swmp

# set date as a factor and order chronologically
colony$Date<-factor(colony$Date,
                    levels=c("June 23rd","July 20th","Aug 3rd",
                             "Aug 23rd", "Aug 31st",  "Sept 27th"),
                    labels = c("June 23rd", "July 20th", "Aug 3rd",
                               "Aug 23rd","Aug 31st",  "Sept 27th")
)

# Colonial cells summary statistics by field day
colony %>%
  group_by(Date) %>%
  get_summary_stats(Colonialcells_L, type = "mean_sd")
#   Date      variable            n       mean         sd
#   <fct>     <fct>           <dbl>      <dbl>      <dbl>
# 1 June 23rd Colonialcells_L    12 114007119. 137578246.
# 2 July 20th Colonialcells_L    12  77067839.  44969107.
# 3 Aug 3rd   Colonialcells_L    12 103367854.  66803026.
# 4 Aug 23rd  Colonialcells_L    12 145333648. 128030418. <- max
# 5 Aug 31st  Colonialcells_L    12  75820922.  54972509. <- min
# 6 Sept 27th Colonialcells_L    11 139085515. 102159236.

# summary statistics by SWMP
colony %>%
  group_by(IDL) %>%
  get_summary_stats(Colonialcells_L, type = "mean_sd")
#   IDL   variable            n       mean         sd
# <chr> <fct>           <dbl>      <dbl>      <dbl>
# 1 A     Colonialcells_L     6 136322256. 190886259.
# 2 B     Colonialcells_L     6  81278986.  45991664.
# 3 C     Colonialcells_L     6 117708734.  37999183.
# 4 D     Colonialcells_L     6 186495782. 117400929.  
# 5 E     Colonialcells_L     5 100752290.  82118992.
# 6 F     Colonialcells_L     6  74298139.  44666186.
# 7 G     Colonialcells_L     6 193066477. 130234432.  <- max
# 8 H     Colonialcells_L     6  72119852.  77262971.
# 9 I     Colonialcells_L     6 101584569.  36547963.
# 10 J    Colonialcells_L     6  45938744.  19736047.
# 11 K    Colonialcells_L     6  37587833.  20479031.  <- min
# 12 L    Colonialcells_L     6 155823263.  96443217.

# create new column w/ log transformed colony forming cell data
colony$logColony<-log10(colony$Colonialcells_L+0.1)

# Repeated measures ANOVA (Within-subject):
res.aov <- anova_test(data = colony, dv = logColony, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)
# ANOVA Table (type III tests)

# Effect DFn DFd     F     p p<.05   ges
# 1     Date   5  50 1.578 0.183       0.092

# determine where the significant differences are:
t.test<-pairwise.t.test(colony$logColony, colony$Date, 
                        p.adjust.method = "bonferroni")
t.test
#Pairwise comparisons using t tests with pooled SD 
# data:  colony$logColony and colony$Date 

#           June 23rd July 20th Aug 3rd Aug 23rd Aug 31st
# July 20th 1.00      -         -       -        -       
# Aug 3rd   1.00      1.00      -       -        -       
# Aug 23rd  1.00      1.00      1.00    -        -       
# Aug 31st  1.00      1.00      1.00    0.78     -       
# Sept 27th 1.00      1.00      1.00    1.00     1.00      

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


# Boxplot total colony forming cyanobacterial cells across field days
colony_box<-ggplot(colony, aes(x=Date,
                               y=Colonialcells_L)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Colonial cyanobacteria cells / L") +
  box_theme
colony_box


# No significant comparisons, do not need to add significance to plot

ggsave("colony_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")
##### Unicellular cells Abundance - Repeated-measures ANOVA#####
uni<-swmp

# set date as a factor and order chronologically
uni$Date<-factor(uni$Date,
                 levels=c("June 23rd","July 20th","Aug 3rd",
                          "Aug 23rd", "Aug 31st",  "Sept 27th"),
                 labels = c("June 23rd", "July 20th", "Aug 3rd",
                            "Aug 23rd","Aug 31st",  "Sept 27th")
)

# Unicellular cells summary statistics by field day
uni %>%
  group_by(Date) %>%
  get_summary_stats(Unicellularcells_L, type = "mean_sd")
#   Date      variable               n    mean      sd
#   <fct>     <fct>              <dbl>   <dbl>   <dbl>
# 1 June 23rd Unicellularcells_L    12 328570. 215330.
# 2 July 20th Unicellularcells_L    12 485054. 668097. <- max
# 3 Aug 3rd   Unicellularcells_L    12  77950. 126637. <- min
# 4 Aug 23rd  Unicellularcells_L    12 123626. 109560.
# 5 Aug 31st  Unicellularcells_L    12 108745.  82483.
# 6 Sept 27th Unicellularcells_L    11 343745. 192936.

# summary statistics by SWMP
uni %>%
  group_by(IDL) %>%
  get_summary_stats(Unicellularcells_L, type = "mean_sd")
#  IDL   variable               n    mean      sd
#<chr> <fct>              <dbl>   <dbl>   <dbl>
# 1 A     Unicellularcells_L     6 315484. 214332.
# 2 B     Unicellularcells_L     6 124312. 110261.
# 3 C     Unicellularcells_L     6 231453. 229504.
# 4 D     Unicellularcells_L     6 196531. 143932.
# 5 E     Unicellularcells_L     5 159539. 147078.
# 6 F     Unicellularcells_L     6 223476. 160040.
# 7 G     Unicellularcells_L     6 132264. 146762.
# 8 H     Unicellularcells_L     6 283031. 377753.
# 9 I     Unicellularcells_L     6 606110. 924876. <- max
# 10 J    Unicellularcells_L     6 205797. 234395.
# 11 K    Unicellularcells_L     6 306729. 267853.
# 12 L    Unicellularcells_L     6 119957.  69171. <- min


# create new column w/ log transformed unicellular cell data
uni$logUni<-log10(uni$Unicellularcells_L+0.1)

# Repeated measures ANOVA (Within-subject):
res.aov <- anova_test(data = uni, dv = logUni, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)
#ANOVA Table (type III tests)

#   Effect  DFn   DFd    F     p p<.05   ges
#1   Date 1.71 17.13 4.02 0.042     * 0.258

# determine where the significant differences are:
t.test<-pairwise.t.test(uni$logUni, uni$Date, 
                        p.adjust.method = "bonferroni")
t.test
#Pairwise comparisons using t tests with pooled SD 
# data:  uni$logUni and uni$Date 

#           June 23rd July 20th Aug 3rd Aug 23rd Aug 31st
#  July 20th 1.000     -         -       -        -       
#  Aug 3rd   0.015     0.016     -       -        -       
#  Aug 23rd  0.357     0.377     1.000   -        -       
#  Aug 31st  1.000     1.000     0.216   1.000    -       
#  Sept 27th 1.000     1.000     0.017   0.370    1.000   

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


# Boxplot total unicellular cyanobacterial cells across field days
uni_box<-ggplot(uni, aes(x=Date,
                         y=Unicellularcells_L)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Unicellular cyanobacteria cells / L") +
  box_theme
uni_box


# add significance to plot

group_uni<-uni %>%
  group_by(Date)%>%
  mutate(max_value = max(Unicellularcells_L))%>%
  mutate(min_value = min(Unicellularcells_L)) %>%
  mutate(uni.mean = mean(Unicellularcells_L))

#subset just min max data
min_max<-unique(group_uni[, c("Date", "max_value", "min_value", "uni.mean")]
)

# generate p-values for labels
t.test.label<-uni %>%
  pairwise_t_test(logUni~Date, p.adjust.method = "bonferroni") 
head(t.test.label,5)

# add matching grouping columns to min_max dataframe
min_max$group1<-min_max$Date
min_max$group2<-min_max$Date

# add group 1 min max values
by<-join_by(group1)
t.test.merge<-left_join(t.test.label, 
                        min_max[,c("group1", "min_value", "max_value", "uni.mean")], 
                        by)

#add group 1 prefix
colnames(t.test.merge)[c(11:13)] <- paste("g1",  colnames(t.test.merge)[c(11:13)], 
                                          sep = '.')

# add group 2 min max values
by<-join_by(group2)
t.test.merge<-left_join(t.test.merge, 
                        min_max[,c("group2", "min_value", "max_value", "uni.mean")], 
                        by)

#add group 2 prefix
colnames(t.test.merge)[c(14:16)] <- paste("g2",  colnames(t.test.merge)[c(14:16)], 
                                          sep = '.')

# add unique column id w/ both dates
t.test.merge$comp<-paste(t.test.merge$group1, sep="-", t.test.merge$group2)


# define overall max, and max mean between dates
t.test.merge <- t.test.merge %>%
  mutate(over_max = pmax(g1.max_value, g2.max_value,  na.rm = TRUE))%>%
  mutate(max_mean = pmax(g1.uni.mean, g2.uni.mean,  na.rm = TRUE)
  )

#subset comparisons we want to include (p < 0.05)
label_sub<-subset(t.test.merge, c(comp == "June 23rd-Aug 3rd" |
                                    comp =="July 20th-Aug 3rd"|
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

# plot with significance - there is no significant differences
uni_box+
  stat_pvalue_manual(label_sub, label = "label", y.position = 1e6, step.increase=0.1)

# Customize y position for each comparison
y <- c(1.3e6, 1.1e6, 1.5e6) 
uni_box+
  stat_pvalue_manual(label_sub, label = "label", y.position=y)

ggsave("uni_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")
