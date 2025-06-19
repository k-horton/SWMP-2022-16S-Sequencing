library(decontam)
library(file2meco)
library(microeco)
library(dplyr)
library(rstatix)
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(writexl)
library(ggplot2)

#### Prepare microtable as before ####
#create microtable 
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
  # drop the positive control sample from the microtable
  table2<-clone(tcontam)
  table2$sample_table<-subset(table2$sample_table, table2$sample_table$Control!="Pos")
  # Remove pollution (mitochondrial and chloroplast sequences) from the microtable
  table2$filter_pollution(taxa = c("mitochondria", "chloroplast"))
  # convert microtable to phyloseq object
  physeq <- meco2phyloseq(table2)
  physeq
}
#and run decontamination 
{sample_data(physeq)$is.neg <- sample_data(physeq)$Control == "Neg"
  contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg", threshold = 0.10)
  # Subset only the sequences designated as contaminants
  prev.contam<-subset(contamdf.prev, contamdf.prev$contaminant == TRUE)
  # Make a list of contaminant sequence IDs
  prev.contam.list<-row.names(prev.contam)
  # Check what the matching taxa are for the contaminant sequence IDs
  prev.contam.taxa<-clone(table2)
  prev.contam.taxa$otu_table<-prev.contam.taxa$otu_table[row.names(prev.contam.taxa$otu_table) %in% c(prev.contam.list),]
  prev.contam.taxa$tidy_dataset()
  prev.taxa.name<-trans_abund$new(dataset = prev.contam.taxa, taxrank="Genus")
  # Remove contaminants from dataframe
  ps.prev.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)
  ps.prev.noncontam
}
#convert phyloseq object back to a microtable and remove control samples
{df <- phyloseq2meco(ps.prev.noncontam)
  no_ctrl <- clone(df)
  no_ctrl$sample_table <- subset(no_ctrl$sample_table, Control == "FALSE")
  no_ctrl$tidy_dataset()}
#### Plot abundance of Cyanobacteria vs all bacteria using microtable ####
# make a copy of the microtable to work with
bact<-clone(no_ctrl)

# Calculate abundance for all bacteria at the Genus level
bact$tax_table <- subset(bact$tax_table, Kingdom == "k__Bacteria")
all_bact_abund <- trans_abund$new(dataset = bact, taxrank = "Genus") 
# Adjust abundance based on initial DNA concentration of sample
all_bact_abund$Abundance<-all_bact_abund$Abundance / all_bact_abund$q_conc

# Make another copy of the microtable, and subset to the Cyanobacteriia class
cyano_mt<-clone(no_ctrl)
cyano_mt$tax_table <- subset(cyano_mt$tax_table, Class == "c__Cyanobacteriia")

# Calculate abundance for all Cyanobacteria at the Genus level
cyano_abund <- trans_abund$new(dataset = cyano_mt, taxrank = "Genus") 
# Adjust abundance based on initial DNA concentration of sample
cyano_abund$Abundance<-cyano_abund$Abundance / cyano_abund$q_conc

# Output Cyanobacteria Genus names to a list
cyano_list<-cyano_abund$data_abund$Taxonomy #1518 rows
un_cyano_list<-unique(cyano_list) #33 unique names
un_cyano_list
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

# Output Bacteria Genus names w/o Cyanobacteria included
bact_list<-clone(no_ctrl)
bact_list$tax_table <- subset(bact_list$tax_table, c(Kingdom == "k__Bacteria" &
                                                       Class != "c__Cyanobacteriia"))
bact_list$tidy_dataset()

bact_abund <- trans_abund$new(dataset = bact_list, taxrank = "Genus") 
bact_abund$Abundance<-bact_abund$Abundance / bact_abund$q_conc
bact_list<-bact_abund$data_abund$Taxonomy #56097 rows
un_bact_list<-unique(bact_list) #813 unique names

# Create lists of new names for cyanobacteria and bacteria in the original dataset
cyano_label<-rep(list("Cyanobacteria"), 33) #33 unique names
bact_label<-rep(list("Other Bacteria"), 813) #813 unique names
# Re-name with lists
all_bact_abund$data_abund$Taxonomy<-factor(all_bact_abund$data_abund$Taxonomy,
                                       levels=c(un_bact_list, un_cyano_list),
                                       labels = c(bact_label, cyano_label))
table(all_bact_abund$data_abund$Taxonomy)
# Other Bacteria  Cyanobacteria 
# 56097           2208

# set colours for each group
my_colours_cyano<-c("Cyanobacteria" = "turquoise",
                   "Other Bacteria"="grey40")

# Define the 'Date' level for the facet-wrap
date<-c("June 23rd", "July 20th","Aug 3rd",
        "Aug 23rd", "Aug 31st","Sept 27th")

# Plot abundance data with each sample represented as a unique bar
abund_bar<- ggplot(all_bact_abund$data_abund, aes(x= IDL,y= Abundance, 
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
abund_bar

ggsave("all_bact_abund_bar.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

# Trim the plot to 50% abundance so it's easier to see abundance of Cyanobacteria
abund_bar+
  coord_cartesian(ylim=c(0,0.40))

ggsave("zoom_all_bact_abund_bar.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

#### Read in joined WQ and micro data dataframe ####
dir<-getwd()
join_data<-read_excel(paste0(dir,"/total_cyano_df.xlsx"))
#### Prepare dataframe and subset to variables of interest ####
# set date as a factor and order chronologically
join_data$Date<-factor(join_data$Date,
                         levels=c("June 23rd", "July 20th","Aug 3rd",
                                  "Aug 23rd", "Aug 31st", "Sept 27th"),
                         labels = c("June 23rd",    "July 20th",  "Aug 3rd",
                                    "Aug 23rd", "Aug 31st", "Sept 27th") )

dataset<-join_data[c(3:5,14,40:42, 50,55,94:96,101)]

#### Summ. Stats - Cyano abund. for DATE and SWMP ID ####
# DATE
summary_stats_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(6:13)])[
  !names(dataset[c(6:13)]) %in% "Date"]

dataset$Date <- as.factor(dataset$Date)

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  summ_stats <- (dataset %>% 
                   group_by(Date) %>%
                   get_summary_stats(variable_name, type = "mean_sd"))

  stats_table <- as.data.frame((summ_stats))
  stats_table$Cat_Variable <- rownames(stats_table)
  stats_table <- stats_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  stats_table$Dependent_Variable <- variable_name
  summary_stats_list[[variable_name]] <- stats_table
}

final_summ_stats <- bind_rows(summary_stats_list)

dir.create(file.path(paste0(dir,"/stats"), "/overall_cyano_stats"), showWarnings = FALSE)
excel_output_path <- paste0(dir,"/stats/overall_cyano_stats/o_cyano_DATE_summ_stats.xlsx")
write_xlsx(final_summ_stats, path = excel_output_path)

# SWMP ID
summary_stats_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(6:13)])[
  !names(dataset[c(6:13)]) %in% "IDL"]

dataset$Date <- as.factor(dataset$Date)

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  summ_stats <- (dataset %>% 
                   group_by(IDL) %>%
                   get_summary_stats(variable_name, type = "mean_sd"))
  
  stats_table <- as.data.frame((summ_stats))
  stats_table$Cat_Variable <- rownames(stats_table)
  stats_table <- stats_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  stats_table$Dependent_Variable <- variable_name
  summary_stats_list[[variable_name]] <- stats_table
}

final_summ_stats <- bind_rows(summary_stats_list)

excel_output_path <- paste0(dir,"/stats/overall_cyano_stats/o_cyano_SWMPID_summ_stats.xlsx")
write_xlsx(final_summ_stats, path = excel_output_path)

#### ANOVA - Cyano abund. for DATE ####
# Initialize an empty list to store ANOVA results
anova_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(6:13)])[
  !names(dataset[c(6:13)]) %in% "Date"]

dataset$Date <- as.factor(dataset$Date)
# set Date label order
date_level<-c("June 23rd","July 20th","Aug 3rd",
              "Aug 23rd", "Aug 31st",  "Sept 27th")

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  # Create the formula for ANOVA
  formula_anova <- as.formula(paste(variable_name, "~ Date"))
  
  # Perform ANOVA
  aov_model <- aov(formula_anova, data = dataset)
  
  anova_table <- as.data.frame(summary(aov_model)[[1]])
  anova_table$Cat_Variable <- rownames(anova_table)
  anova_table <- anova_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  anova_table$Dependent_Variable <- variable_name
  anova_results_list[[variable_name]] <- anova_table
}

final_anova_df <- bind_rows(anova_results_list)

# remove extra blank spaces in the variable names
final_anova_df$Cat_Variable<-factor(final_anova_df$Cat_Variable,
                                    levels = c("Date       ","Residuals  "),
                                    labels = c("Date","Residuals"))

# drop observations for residuals
final_anova_df_2<-final_anova_df[final_anova_df$Cat_Variable == "Date",]


excel_output_path <- paste0(dir,"/stats/overall_cyano_stats/o_cyano_DATE_anova.xlsx")
write_xlsx(final_anova_df_2, path = excel_output_path)

#### T-test - Cyano abund. for FLOW ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(6:13)])[
  !names(dataset[c(6:13)]) %in% "Storm_Base"]

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  # Create the formula for T-test
  formula_t.test <- as.formula(paste(variable_name, "~ Storm_Base"))
  
  # Perform T-test
  t.test<- t.test(formula_t.test, data = dataset)
  
  t.test_table <- as.data.frame((t.test[1:4]))
  t.test_table$Cat_Variable <- "Storm_Base"
  t.test_table <- t.test_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  t.test_table$Dependent_Variable <- variable_name
  t.test_results_list[[variable_name]] <- t.test_table
}

final_t.test_df <- bind_rows(t.test_results_list)

# drop observations for residuals
final_t.test_df_2<-final_t.test_df[final_t.test_df$Cat_Variable == "Storm_Base",]

excel_output_path <- paste0(dir,"/stats/overall_cyano_stats/o_cyano_FLOW_t.test.xlsx")
write_xlsx(final_t.test_df_2, path = excel_output_path)

#### T-test - Cyano abund. for DREDGE ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(6:13)])[
  !names(dataset[c(6:13)]) %in% "Dredge_Cat"]

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  # Create the formula for T-test
  formula_t.test <- as.formula(paste(variable_name, "~ Dredge_Cat"))
  
  # Perform T-test
  t.test<- t.test(formula_t.test, data = dataset)
  
  t.test_table <- as.data.frame((t.test[1:4]))
  t.test_table$Cat_Variable <- "Dredge_Cat"
  t.test_table <- t.test_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  t.test_table$Dependent_Variable <- variable_name
  t.test_results_list[[variable_name]] <- t.test_table
}

final_t.test_df <- bind_rows(t.test_results_list)

# drop observations for residuals
final_t.test_df_2<-final_t.test_df[final_t.test_df$Cat_Variable == "Dredge_Cat",]

excel_output_path <- paste0(dir,"/stats/overall_cyano_stats/o_cyano_DREDGE_t.test.xlsx")
write_xlsx(final_t.test_df_2, path = excel_output_path)
#### T-test - Cyano abund. for PRE_2003 ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(6:13)])[
  !names(dataset[c(6:13)]) %in% "Pre_2003"]

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  # Create the formula for T-test
  formula_t.test <- as.formula(paste(variable_name, "~ Pre_2003"))
  
  # Perform T-test
  t.test<- t.test(formula_t.test, data = dataset)
  
  t.test_table <- as.data.frame((t.test[1:4]))
  t.test_table$Cat_Variable <- "Pre_2003"
  t.test_table <- t.test_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  t.test_table$Dependent_Variable <- variable_name
  t.test_results_list[[variable_name]] <- t.test_table
}

final_t.test_df <- bind_rows(t.test_results_list)

# drop observations for residuals
final_t.test_df_2<-final_t.test_df[final_t.test_df$Cat_Variable == "Pre_2003",]

excel_output_path <- paste0(dir,"/stats/overall_cyano_stats/o_cyano_2003_t.test.xlsx")
write_xlsx(final_t.test_df_2, path = excel_output_path)

#### Define aesthetics for plots #####
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
                   legend.text = element_text(size = 20),
                   legend.key=element_blank()
)

# set Date label order
date_level<-c("June 23rd","July 20th","Aug 3rd",
              "Aug 23rd", "Aug 31st",  "Sept 27th")

my_colours_storm<-c("Stormflow"= "steelblue4",
                    "Low Flow" = "slategray1")
#### Plot - Abundance from 16S sequencing #####
cyano_box<-ggplot(dataset, aes(x=Date,
                             y=Abundance)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Relative Abundance of Cyanobacteria (%)") +
  box_theme
cyano_box

# Define significance labels for plot
{group_cyano<-dataset %>%
  group_by(Date)%>%
  mutate(max_value = max(Abundance))%>%
  mutate(min_value = min(Abundance)) %>%
  mutate(abund.mean = mean(Abundance))

#subset just min max data
min_max<-unique(group_cyano[, c("Date", "max_value", "min_value", "abund.mean")])

# generate p-values for labels
t.test.label<-dataset %>%
  pairwise_t_test(logAbund~Date, p.adjust.method = "holm") 
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
t.test.merge<-left_join(t.test.merge, min_max[,c("group2", "min_value", "max_value", "abund.mean")], 
                        by)

#add group 2 prefix
colnames(t.test.merge)[c(14:16)] <- paste("g2",  colnames(t.test.merge)[c(14:16)], 
                                          sep = '.')
# add unique column id w/ both dates
t.test.merge$comp<-paste(t.test.merge$group1, sep="-", t.test.merge$group2)

# define overall max, and max mean between dates
t.test.merge <- t.test.merge %>%
  mutate(over_max = pmax(g1.max_value, g2.max_value,  na.rm = TRUE))%>%
  mutate(max_mean = pmax(g1.abund.mean, g2.abund.mean,  na.rm = TRUE) )

#subset comparisons we want to include (p < 0.05)
label_sub<-subset(t.test.merge, p.adj < 0.05)              

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
                                            ) ) )  )}
                       
# Plot with significance
cyano_box+
  stat_pvalue_manual(label_sub, label = "label", y.position = 80, step.increase=0.1)

# Customize y position for each comparison
y <- c(195, 75, 60, 210, 130, 115, 100) 
cyano_box+
  stat_pvalue_manual(label_sub, label = "label", y.position=y)

ggsave("cyano_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")
#### Plot - Abundance from microscopy#####
t_cyano_box<-ggplot(dataset, aes(x=Date,
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

#### Plot - Abundance from qPCR targeting 16S ####
q16S_box<-ggplot(dataset, aes(x=Date,
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

#### Plot - Abundance of all algae from Chl-a#####
chla_box<-ggplot(dataset, aes(x=Date,
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
