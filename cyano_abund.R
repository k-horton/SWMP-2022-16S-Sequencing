library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(writexl)
library(ggplot2)

#### Read in data ####
dir<-getwd()
genus_data<-read_excel(paste0(dir,"/wq_genus_cyano_ASV_abundance.xlsx"))
family_data<-read_excel(paste0(dir,"/wq_family_cyano_ASV_abundance.xlsx"))
cyano_data<-read_excel(paste0(dir,"/wq_cyano_ASV_abundance.xlsx"))

#### Plot abundance of Cyanobacteria vs all bacteria ####
# set colours for each group
my_colours_cyano<-c("Cyanobacteria" = "turquoise",
                    "Other Bacteria"="grey40")

# Define the 'Date' level for the facet-wrap
date<-c("June 23rd", "July 20th","Aug 3rd",
        "Aug 23rd", "Aug 31st","Sept 27th")

# pivot data
cyano_pivot<-pivot_longer(cyano_data, cols=c(7,8), 
                          names_to = "Group",
                          values_to="Abundance")

cyano_pivot$Group<-factor(cyano_pivot$Group,
                          levels = c("Prop_abund_bact", "Prop_abund_cyano"),
                          labels = c("Other Bacteria", "Cyanobacteria"))
# Plot abundance data with each sample represented as a unique bar
abund_bar<- ggplot(cyano_pivot, aes(x= IDL,y= Abundance*100, 
                                                  fill = Group)) + 
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=my_colours_cyano, name = "Group") + 
  xlab("SWMP") +
  ylab("Relative Abundance (%)") +
  facet_wrap(factor(cyano_pivot$Date, level=date))+
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

#### Summ. Stats - Cyano abund. for DATE and SWMP ID ####
# Subset to variables of interest 
dataset<-cyano_data[c(2:8,11,20,46:48,56,97,98,103)]

# DATE
summary_stats_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(3:7,11:16)])[
  !names(dataset[c(3:7,11:16)]) %in% "Date"]

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
numerical_cols <- names(dataset[c(3:7,11:16)])[
  !names(dataset[c(3:7,11:16)]) %in% "IDL"]

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
# re-name Other Bacteria column
dataset<-rename(dataset, Other_Bacteria = `Other Bacteria`)
# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(3:7,11:16)])[
  !names(dataset[c(3:7,11:16)]) %in% "Date"]

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
numerical_cols <- names(dataset[c(3:7,11:16)])[
  !names(dataset[c(3:7,11:16)]) %in% "Storm_Base"]

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

excel_output_path <- paste0(dir,"/stats/overall_cyano_stats/o_cyano_FLOW_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)

#### T-test - Cyano abund. for DREDGE ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(3:7,11:16)])[
  !names(dataset[c(3:7,11:16)]) %in% "Dredge_Cat"]

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
write_xlsx(final_t.test_df, path = excel_output_path)
#### T-test - Cyano abund. for PRE_2003 ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset[c(3:7,11:16)])[
  !names(dataset[c(3:7,11:16)]) %in% "Pre_2003"]

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

excel_output_path <- paste0(dir,"/stats/overall_cyano_stats/o_cyano_2003_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)

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
dataset<-cyano_data[c(2:4,7:11,20,46:48,56,97:98,103)]

cyano_box<-ggplot(dataset, aes(x=factor(x=Date, levels=date_level),
                               y=Prop_abund_cyano*100)) + 
  geom_boxplot(fill = "lightgray") + 
  xlab("Field Day") +
  ylab("Relative Abundance of Cyanobacteria (%)") +
  box_theme
cyano_box

# Define significance labels for plot
{# drop na's
  group_drop<-drop_na(dataset)
  
  group_cyano<-group_drop %>%
    group_by(Date)%>%
    mutate(max_value = max(Prop_abund_cyano*100))%>%
    mutate(min_value = min(Prop_abund_cyano*100)) %>%
    mutate(abund.mean = mean(Prop_abund_cyano*100))
  
  #subset just min max data
  min_max<-unique(group_cyano[, c("Date", "max_value", "min_value", "abund.mean")])
  
  # generate p-values for labels
  t.test.label<-dataset %>%
    pairwise_t_test(logAbund_cyano~Date, p.adjust.method = "holm") 
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
  colnames(t.test.merge)[c(10:12)] <- paste("g1",  colnames(t.test.merge)[c(10:12)], 
                                            sep = '.')
  # add group 2 min max values
  by<-join_by(group2)
  t.test.merge<-left_join(t.test.merge, min_max[,c("group2", "min_value", "max_value", "abund.mean")], 
                          by)
  
  #add group 2 prefix
  colnames(t.test.merge)[c(13:15)] <- paste("g2",  colnames(t.test.merge)[c(13:15)], 
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
  stat_pvalue_manual(label_sub, label = "label", y.position = 60, step.increase=0.1)

# Customize y position for each comparison
y <- c(48, 69, 74, 53, 79) 
cyano_box+
  stat_pvalue_manual(label_sub, label = "label", y.position=y)

ggsave("cyano_abund_box.png",  plot = last_plot(), 
       path = paste0(dir, "/Abundance_plots"),
       width = 14, height = 8, units = "in")
#### Plot - Abundance from microscopy#####
t_cyano_box<-ggplot(dataset, aes(x=factor(x=Date, levels=date_level),
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
q16S_box<-ggplot(dataset, aes(x=factor(x=Date, levels=date_level),
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
chla_box<-ggplot(dataset, aes(x=factor(x=Date, levels=date_level),
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
