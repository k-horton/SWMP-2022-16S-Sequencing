library(dplyr)
library(EnvStats)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(readxl)
library(rstatix)
library(writexl)

#### Read in dataframes ####
dir<-getwd()
genus_data<-read_excel(paste0(dir,"/wq_genus_cyano_ASV_abundance.xlsx"))
allcyano_data<-read_excel(paste0(dir,"/wq_overall_cyano_ASV_abundance.xlsx"))
#### Prepare dataframes and subset to variables of interest ####
# set date as a factor and order chronologically
allcyano_data$Date<-factor(allcyano_data$Date,
                       levels=c("June 23rd", "July 20th","Aug 3rd",
                                "Aug 23rd", "Aug 31st", "Sept 27th"),
                       labels = c("June 23rd",    "July 20th",  "Aug 3rd",
                                  "Aug 23rd", "Aug 31st", "Sept 27th") )

allcyano<-allcyano_data[c(1:3,7:12,21:38,47:49,51:58,62:104)]


#### Determine # of potentially toxic reads per sample ####
toxic_reads<-genus_data[c(1:14)]

mcyE_only<-drop_na(toxic_reads[toxic_reads$seq_mcyE_producer=="Yes",])
mcyE_only_sum<-aggregate(genus_prop_abund ~ Sample, data = mcyE_only, FUN = sum)
mcyE_only_sum<-rename(mcyE_only_sum, "seq_pot_mcyE"="genus_prop_abund")

oth_tox<-drop_na(toxic_reads[toxic_reads$seq_oth_tox=="Yes",])
oth_tox_sum<-aggregate(genus_prop_abund ~ Sample, data = oth_tox, FUN = sum)
oth_tox_sum<-rename(oth_tox_sum, "seq_pot_oth_tox"="genus_prop_abund")

all_tox<-drop_na(toxic_reads[(toxic_reads$seq_oth_tox=="Yes" | 
                                toxic_reads$seq_mcyE_producer=="Yes" ),])
all_tox_sum<-aggregate(genus_prop_abund ~ Sample, data = all_tox, FUN = sum)
all_tox_sum<-rename(all_tox_sum, "seq_pot_all_tox"="genus_prop_abund")

merge_table<-merge(mcyE_only_sum, oth_tox_sum, by="Sample", all=TRUE)
merge_table2<-merge(merge_table, all_tox_sum)

excel_output_path <- paste0(dir,"/seq_tox_prod_sum.xlsx")
write_xlsx(merge_table2, path = excel_output_path)

#### Summ. Stats - mcyE copies for DATE and SWMP ID ####
# DATE
summary_stats_list <- list()

# Get the names of the numerical variables you want to include
numerical_cols <- names(allcyano[c(31,33,76,81)])[
  !names(allcyano[c(31,33,76,81)]) %in% "Date"]

allcyano$Date <- as.factor(allcyano$Date)

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  summ_stats <- (allcyano %>% 
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

dir.create(file.path(paste0(dir,"/stats"), "/toxin_cyano_stats"), showWarnings = FALSE)
excel_output_path <- paste0(dir,"/stats/toxin_cyano_stats/mcyE_DATE_summ_stats.xlsx")
write_xlsx(final_summ_stats, path = excel_output_path)

# SWMP ID
summary_stats_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(allcyano[c(31,33,76,81)])[
  !names(allcyano[c(31,33,76,81)]) %in% "IDL"]

allcyano$Date <- as.factor(allcyano$Date)

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  summ_stats <- (allcyano %>% 
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

excel_output_path <- paste0(dir,"/stats/toxin_cyano_stats/mcyE_SWMPID_summ_stats.xlsx")
write_xlsx(final_summ_stats, path = excel_output_path)


#### ANOVA - mcyE copies/L + potential toxin producers (sequencing) for DATE ####
# Initialize an empty list to store ANOVA results
dataset<-merge(allcyano, merge_table2, by="Sample", all=TRUE)

dataset2<-dataset %>%
  mutate(seq_pot_mcyE = case_when(is.na(seq_pot_mcyE) & is.na(cyano_library_size) ~ NA,
                                  is.na(seq_pot_mcyE) & cyano_library_size >0  ~ 0,
                                  seq_pot_mcyE > 0 & cyano_library_size >0  ~ seq_pot_mcyE))%>%
  mutate(seq_pot_oth_tox = case_when(is.na(seq_pot_oth_tox) & is.na(cyano_library_size) ~ NA,
                                     is.na(seq_pot_oth_tox) & cyano_library_size >0  ~ 0,
                                     seq_pot_oth_tox > 0 & cyano_library_size >0  ~ seq_pot_oth_tox))%>%
  mutate(seq_pot_all_tox = case_when(is.na(seq_pot_all_tox) & is.na(cyano_library_size) ~ NA,
                                     is.na(seq_pot_all_tox) & cyano_library_size >0  ~ 0,
                                     seq_pot_all_tox > 0 & cyano_library_size >0  ~ seq_pot_all_tox))

dataset2$logseq_pot_mcyE <- log10(dataset2$seq_pot_mcyE+0.001)
dataset2$logseq_pot_oth_tox <- log10(dataset2$seq_pot_oth_tox+0.001)
dataset2$logseq_pot_all_tox <- log10(dataset2$seq_pot_all_tox+0.001)

anova_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset2[c(31,33,76,81:87)])[
  !names(dataset2[c(31,33,76,81:87)]) %in% "Date"]

dataset2$Date <- as.factor(dataset2$Date)
# set Date label order
date_level<-c("June 23rd","July 20th","Aug 3rd",
              "Aug 23rd", "Aug 31st",  "Sept 27th")


# Loop through each numerical variable
for (variable_name in numerical_cols) {
  # Create the formula for ANOVA
  formula_anova <- as.formula(paste(variable_name, "~ Date"))
  
  # Perform ANOVA
  aov_model <- aov(formula_anova, data = dataset2)
  
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

excel_output_path <- paste0(dir,"/stats/toxin_cyano_stats/mcyE_DATE_anova.xlsx")
write_xlsx(final_anova_df_2, path = excel_output_path)

#### T-test - all quantitative data for mcyE pres/abs  ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(allcyano[c(11:27,29,30, 34:37,39,41:80)])[
  !names(allcyano[c(11:27,29,30, 34:37,39,41:80)]) %in% "mcyE_Pres"]

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  # Create the formula for T-test
  formula_t.test <- as.formula(paste(variable_name, "~ mcyE_Pres"))
  
  # Perform T-test
  t.test<- t.test(formula_t.test, data = allcyano)
  
  t.test_table <- as.data.frame((t.test[1:4]))
  t.test_table$Cat_Variable <- "mcyE_Pres"
  t.test_table <- t.test_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  t.test_table$Dependent_Variable <- variable_name
  t.test_results_list[[variable_name]] <- t.test_table
}

final_t.test_df <- bind_rows(t.test_results_list)

excel_output_path <- paste0(dir,"/stats/toxin_cyano_stats/mcyE_PRES_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)
#### T-test - mcyE copies and all categorical variables ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the categorical variables you want to test
categorical_cols <- names(allcyano[c(9:10,28,38)]
                          )[!names(allcyano[c(9:10,28,38)]) %in% "logmcyE"]

allcyano$Storm_Base<-as.factor(allcyano$Storm_Base)
allcyano$Dredge_Cat<-as.factor(allcyano$Dredge_Cat)
allcyano$Pre_2003<-as.factor(allcyano$Pre_2003)
allcyano$All_Susp_Toxic<-as.factor(allcyano$All_Susp_Toxic)

# Loop through each categorical variable
for (variable_name in categorical_cols) {
  # Create the formula for T-test
  formula_t.test <- as.formula(paste("logmcyE ~", variable_name))
  
  # Perform T-test
  t.test<- t.test(formula_t.test, data = allcyano)
  
  t.test_table <- as.data.frame((t.test[1:4]))
  t.test_table$Num_Variable <- "logmcyE"
  t.test_table <- t.test_table %>%
    select(Num_Variable, everything()) # Reorder columns
  
  t.test_table$Dependent_Variable <- variable_name
  t.test_results_list[[variable_name]] <- t.test_table
}

final_t.test_df <- bind_rows(t.test_results_list)

excel_output_path <- paste0(dir,"/stats/toxin_cyano_stats/mcyE_allCAT_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)
#### Pearson - mcyE copies + toxic reads x WQ and abundance data ####
# Subset to variables of interest 
dataset<-allcyano_data[c(1:3,7:12,21:38,47:49,51,53:57,62:104,106)]
dataset<-merge(dataset, merge_table2, by="Sample", all=TRUE)

dataset2<-dataset %>%
  mutate(seq_pot_mcyE = case_when(is.na(seq_pot_mcyE) & is.na(cyano_library_size) ~ NA,
                         is.na(seq_pot_mcyE) & cyano_library_size >0  ~ 0,
                         seq_pot_mcyE > 0 & cyano_library_size >0  ~ seq_pot_mcyE))%>%
  mutate(seq_pot_oth_tox = case_when(is.na(seq_pot_oth_tox) & is.na(cyano_library_size) ~ NA,
                                  is.na(seq_pot_oth_tox) & cyano_library_size >0  ~ 0,
                                  seq_pot_oth_tox > 0 & cyano_library_size >0  ~ seq_pot_oth_tox))%>%
  mutate(seq_pot_all_tox = case_when(is.na(seq_pot_all_tox) & is.na(cyano_library_size) ~ NA,
                                  is.na(seq_pot_all_tox) & cyano_library_size >0  ~ 0,
                                  seq_pot_all_tox > 0 & cyano_library_size >0  ~ seq_pot_all_tox))

dataset2$logseq_pot_mcyE <- log10(dataset2$seq_pot_mcyE+0.001)
dataset2$logseq_pot_oth_tox <- log10(dataset2$seq_pot_oth_tox+0.001)
dataset2$logseq_pot_all_tox <- log10(dataset2$seq_pot_all_tox+0.001)

# Get the names of the toxin and WQ variables being tested
toxin_cols <- colnames(dataset2)[c(31,32,74,79,81:86)]
toxin_cols <- toxin_cols[!toxin_cols %in% "Sample"]

wq_cols <- colnames(dataset2)[c(4:8,11:27,29,30,33:73,75:78,80)]
wq_cols <- wq_cols[!wq_cols %in% "Sample"]

pearson_results_list <- list()

for (toxin_variable in toxin_cols) {
  pearson_table <- data.frame(
    wq_variable=numeric(),
    cor=numeric(),
    p=numeric()
  )
  for (wq_variable in wq_cols) {
    # Perform Pearson correlation
    pearson <- cor_test(toxin_variable,wq_variable, data=dataset2, method="pearson",
                        alternative="two.sided")
    pearson<-as.data.frame(pearson)
    pearson_table[nrow(pearson_table) + 1,] <- list(pearson$var2, pearson$cor, pearson$p)    
    pearson_results_list[[toxin_variable]] <- pearson_table
  }
}


excel_output_path <- paste0(dir,"/stats/toxin_cyano_stats/toxin_PEARSON.xlsx")
write_xlsx(pearson_results_list, path = excel_output_path)


