# Assess cyanobacterial community - specific taxa
library(dplyr)
library(rstatix)
library(ggplot2)
library(EnvStats)
library(tidyr)
library(readxl)
library(writexl)

#### Read in genus data  ####
dir<-getwd()

genus_data<-read_excel(paste0(dir,"/wq_genus_cyano_ASV_abundance.xlsx"))
cyano_data<-read_excel(paste0(dir,"/wq_overall_cyano_ASV_abundance.xlsx"))

genus_data$Date<-factor(genus_data$Date,
                     levels=c("June 23rd","July 20th","Aug 3rd",
                              "Aug 23rd", "Aug 31st",  "Sept 27th"),
                     labels = c("June 23rd", "July 20th", "Aug 3rd",
                                "Aug 23rd","Aug 31st",  "Sept 27th"))

cyano_data$Date<-factor(cyano_data$Date,
                        levels=c("June 23rd","July 20th","Aug 3rd",
                                 "Aug 23rd", "Aug 31st",  "Sept 27th"),
                        labels = c("June 23rd", "July 20th", "Aug 3rd",
                                   "Aug 23rd","Aug 31st",  "Sept 27th"))
#### Relative abundance tables for genus data ####

genus_pivot<-pivot_wider(genus_data[c(2,4,13)], 
                   id_cols = c("Genus"), 
                   names_from = "Date", names_sort = TRUE,
                   values_from = "genus_prop_abund",
                   values_fn=mean)

excel_output_path <- paste0(dir,"/stats/cyano_abund_DATE_16S_sequencing.xlsx")
write_xlsx(genus_pivot, path = excel_output_path)


genus_pivot<-pivot_wider(genus_data[c(3,4,13)], 
                         id_cols = c("Genus"), 
                         names_from = "IDL", names_sort = TRUE,
                         values_from = "genus_prop_abund",
                         values_fn=mean)

excel_output_path <- paste0(dir,"/stats/cyano_abund_SWMPID_16S_sequencing.xlsx")
write_xlsx(genus_pivot, path = excel_output_path)

#### subset data ####
dataset <-cyano_data[c(2:3,12,21,47,54:57,101:104)]
# set date as a factor and order chronologically
dataset$Date<-factor(dataset$Date,
                   levels=c("June 23rd","July 20th","Aug 3rd",
                            "Aug 23rd", "Aug 31st",  "Sept 27th"),
                   labels = c("June 23rd", "July 20th", "Aug 3rd",
                              "Aug 23rd","Aug 31st",  "Sept 27th"))
#### Define aesthetics and labels ####
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

my_colours_celltype<- c("Unicellular" = "aquamarine", 
                        "Colonial" = "cyan3", 
                        "Filamentous" ="turquoise4")

bar_theme <- theme(strip.background = element_rect(fill="white", color="black"),
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

#### Summ. Stats - abundance for each cell type for DATE and SWMP ID ####
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

dir<-getwd()
dir.create(file.path(paste0(dir,"/stats"), "/m_community_stats"), showWarnings = FALSE)
excel_output_path <- paste0(dir,"/stats/m_community_stats/m_community_DATE_summ_stats.xlsx")
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

excel_output_path <- paste0(dir,"/stats/m_community_stats/m_community_SWMPID_summ_stats.xlsx")
write_xlsx(final_summ_stats, path = excel_output_path)

#### ANOVA - Each Cell-type for DATE ####
#subset dataframe
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


excel_output_path <- paste0(dir,"/stats/m_community_stats/m_community_DATE_anova.xlsx")
write_xlsx(final_anova_df_2, path = excel_output_path)

#### Plots - Filamentous cells #####
# Boxplot total filamentous cyanobacterial cells across field days
filam_box<-ggplot(dataset, aes(x=Date,
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

#### Rosner Outlier test - -Filamentous cells #####
ros_fil<-dataset

#plot filamentous cells without any grouping
ros_filam_box<-ggplot(dataset, aes(y=Filamentouscells_L)) + 
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
res.aov <- anova_test(data = ros_fil, dv = logF_cyano, 
                      wid = IDL, within = Date)
get_anova_table(res.aov)

t.test<-pairwise.t.test(ros_fil$logF_cyano, ros_fil$Date, 
                        p.adjust.method = "holm")
t.test
# there is not much change in the reusults even with the outliers removed,
#     still no signicant differences in filamentous cells between field days

#### Plots - Colonial cells #####
# Boxplot total colony forming cyanobacterial cells across field days
colony_box<-ggplot(dataset, aes(x=Date,
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
#### Plots - Unicellular cells #####
# Boxplot total unicellular cyanobacterial cells across field days
uni_box<-ggplot(dataset, aes(x=Date,
                         y=Unicellularcells_L)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Unicellular cyanobacteria cells / L") +
  box_theme
uni_box

ggsave("uni_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")

#### Convert cell type to categorical variable #####
cyano_long<-dataset%>% 
  pivot_longer(
    cols = c("Unicellularcells_L", "Colonialcells_L", "Filamentouscells_L"), 
    names_to = "Cyano_cell_type", 
    values_to = "Cyano_cell_abund"
  )

cyano_long$logCyano_cell_abund<-log10(cyano_long$Cyano_cell_abund+0.001)

#### ANOVA and Pairwise T-test - for cell type as cat. variable across field days ####
# Repeated measures ANOVA (Within-subject):
res.aov <- anova_test(data = cyano_long, dv = logCyano_cell_abund, 
                      wid = IDL, within = c(Date, Cyano_cell_type))
get_anova_table(res.aov)
#ANOVA Table (type III tests)
#                Effect  DFn  DFd      F        p p<.05   ges
# 1                 Date  5.0 40.0  1.889 1.18e-01       0.056
# 2      Cyano_cell_type  1.1  8.8 50.227 5.25e-05     * 0.496
# 3 Date:Cyano_cell_type 10.0 80.0  1.332 2.28e-01       0.074

anova_table <- as.data.frame((res.aov$ANOVA))

# determine where the significant differences are:
t.test<-cyano_long %>%
  group_by(Date) %>%
  pairwise_t_test(data =., logCyano_cell_abund ~Cyano_cell_type) %>%
  adjust_pvalue(method="holm")
t.test
# Date      .y.    group1 group2    n1    n2       p p.signif   p.adj p.adj.signif
# <chr>     <chr>  <chr>  <chr>  <int> <int>   <dbl> <chr>      <dbl> <chr>       
# 1 Aug 23rd  logCy… Colon… Filam…    11    11 1.28e-4 ***      1.54e-3 ***         
# 2 Aug 23rd  logCy… Colon… Unice…    11    11 3.85e-3 **       3.47e-2 **          
# 3 Aug 23rd  logCy… Filam… Unice…    11    11 2.17e-1 ns       2.17e-1 ns          
# 4 Aug 31st  logCy… Colon… Filam…    12    12 5.09e-6 ****     7.63e-5 ****        
# 5 Aug 31st  logCy… Colon… Unice…    12    12 1.75e-2 *        1.05e-1 *           
# 6 Aug 31st  logCy… Filam… Unice…    12    12 6.03e-3 **       4.22e-2 *           
# 7 Aug 3rd   logCy… Colon… Filam…    11    11 2.95e-7 ****     4.72e-6 ****        
# 8 Aug 3rd   logCy… Colon… Unice…    11    11 1.44e-3 **       1.44e-2 **          
# 9 Aug 3rd   logCy… Filam… Unice…    11    11 4.73e-3 **       3.78e-2 **          
# 10 July 20th logCy… Colon… Filam…    12    12 1.41e-7 ****     2.54e-6 ****        
# 11 July 20th logCy… Colon… Unice…    12    12 3.38e-2 *        1.35e-1 *           
# 12 July 20th logCy… Filam… Unice…    12    12 9.43e-5 ****     1.32e-3 ***         
# 13 June 23rd logCy… Colon… Filam…    12    12 2.36e-7 ****     4.01e-6 ****        
# 14 June 23rd logCy… Colon… Unice…    12    12 3.44e-2 *        1.35e-1 *           
# 15 June 23rd logCy… Filam… Unice…    12    12 1.54e-4 ***      1.69e-3 ***         
# 16 Sept 27th logCy… Colon… Filam…    11    11 1.08e-4 ***      1.40e-3 ***         
# 17 Sept 27th logCy… Colon… Unice…    11    11 2.07e-2 *        1.05e-1 *           
# 18 Sept 27th logCy… Filam… Unice…    11    11 5.31e-2 ns       1.35e-1 ns  

# output T-test and ANOVA results
excel_output_path <- paste0(dir,"/stats/m_community_stats/m_community_CELLTYPE_t.test.xlsx")
write_xlsx(t.test, path = excel_output_path)

excel_output_path <- paste0(dir,"/stats/m_community_stats/m_community_CELLTYPE_anova.xlsx")
write_xlsx(anova_table, path = excel_output_path)

# re-name levels for cell type
cyano_long$Cyano_cell_type<-factor(cyano_long$Cyano_cell_type ,
                        levels = c("Colonialcells_L","Unicellularcells_L",  "Filamentouscells_L"),
                        labels = c("Colonial","Unicellular",  "Filamentous"))

#### T-test - Cell-type for FLOW ####
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

excel_output_path <- paste0(dir,"/stats/m_community_stats/m_community_FLOW_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)

#### T-test - Cell-type for DREDGE ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(6:13)])[
  !names(dataset[c(6:13)]) %in% "Date"]

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

excel_output_path <- paste0(dir,"/stats/m_community_stats/m_community_DREDGE_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)

#### T-test - Cell-type for PRE_2003 ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(6:13)])[
  !names(dataset[c(6:13)]) %in% "Date"]

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

excel_output_path <- paste0(dir,"/stats/m_community_stats/m_community_2003_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)
