library(dplyr)
library(rstatix)
library(ggplot2)
library(phyloseq)
library(EnvStats)
library(ggpubr)
library(writexl)

#### Read in dataframes ####
dir<-getwd()
genus_data<-read_excel(paste0(dir,"/genus_cyano_df.xlsx"))
allcyano_data<-read_excel(paste0(dir,"/total_cyano_df.xlsx"))
#### Prepare dataframes and subset to variables of interest ####
# set date as a factor and order chronologically
allcyano_data$Date<-factor(allcyano_data$Date,
                       levels=c("June 23rd", "July 20th","Aug 3rd",
                                "Aug 23rd", "Aug 31st", "Sept 27th"),
                       labels = c("June 23rd",    "July 20th",  "Aug 3rd",
                                  "Aug 23rd", "Aug 31st", "Sept 27th") )

allcyano<-allcyano_data[c(3:5,14:31,40:42,44:51,55,56,59:103)]

#### Summ. Stats - mcyE copies for DATE and SWMP ID ####
# DATE
summary_stats_list <- list()

# Get the names of the numerical variables you want to include
numerical_cols <- names(allcyano[c(25,27,74,79)])[
  !names(allcyano[c(25,27,74,79)]) %in% "Date"]

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
numerical_cols <- names(allcyano[c(25,27,74,79)])[
  !names(allcyano[c(25,27,74,79)]) %in% "IDL"]

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


#### ANOVA - mcyE copies/L for DATE ####
# Initialize an empty list to store ANOVA results
anova_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(allcyano[c(25,27,74,79)])[
  !names(allcyano[c(25,27,74,79)]) %in% "Date"]

allcyano$Date <- as.factor(allcyano$Date)
# set Date label order
date_level<-c("June 23rd","July 20th","Aug 3rd",
              "Aug 23rd", "Aug 31st",  "Sept 27th")


# Loop through each numerical variable
for (variable_name in numerical_cols) {
  # Create the formula for ANOVA
  formula_anova <- as.formula(paste(variable_name, "~ Date"))
  
  # Perform ANOVA
  aov_model <- aov(formula_anova, data = allcyano)
  
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

write_xlsx(final_anova_df_2, path = excel_output_path)
excel_output_path <- paste0(dir,"/stats/toxin_cyano_stats/mcyE_DATE_anova.xlsx")

#### T-test - all quantitative data for mcyE pres/abs  ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(allcyano[c(5:21,23,24, 28:31,33,35:73,75:78)])[
  !names(allcyano[c(5:21,23,24, 28:31,33,35:73,75:78)]) %in% "mcyE_Pres"]

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
categorical_cols <- names(allcyano[c(3:4,22,32,34)]
                          )[!names(allcyano[c(3:4,22,32,34)]) %in% "logmcyE"]

allcyano$Storm_Base<-as.factor(allcyano$Storm_Base)
allcyano$Dredge_Cat<-as.factor(allcyano$Dredge_Cat)
allcyano$Pre_2003<-as.factor(allcyano$Pre_2003)
allcyano$Seq_MC_pres<-as.factor(allcyano$Seq_MC_pres)
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