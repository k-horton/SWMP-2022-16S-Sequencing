library(dplyr)
library(EnvStats)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(readxl)
library(writexl)
library(rstatix)

#### Read in dataframes ####
dir<-getwd()
genus_data<-read_excel(paste0(dir,"/wq_genus_cyano_ASV_abundance.xlsx"))
cyano_data<-read_excel(paste0(dir,"/wq_cyano_ASV_abundance.xlsx"))
#### Prepare dataframes and subset to variables of interest ####
# re-name Other Bacteria column
cyano_data<-rename(cyano_data, Other_Bacteria = `Other Bacteria`)

# set date as a factor and order chronologically
cyano_data$Date<-factor(cyano_data$Date,
                           levels=c("June 23rd", "July 20th","Aug 3rd",
                                    "Aug 23rd", "Aug 31st", "Sept 27th"),
                           labels = c("June 23rd",    "July 20th",  "Aug 3rd",
                                      "Aug 23rd", "Aug 31st", "Sept 27th") )

cyano_sub<-cyano_data[c(2,3,7:11,20,46,58:60,105)]

genus_data$Date<-factor(genus_data$Date,
                        levels=c("June 23rd", "July 20th","Aug 3rd",
                                 "Aug 23rd", "Aug 31st", "Sept 27th"),
                        labels = c("June 23rd",    "July 20th",  "Aug 3rd",
                                   "Aug 23rd", "Aug 31st", "Sept 27th") )

genus_sub<-genus_data[c(2:7,13:15,24,50,62:64,109)]

#### Determine # of potential N-fixers per sample ####
n_fix_reads<-genus_data[c(1:14)]

t1_only<-drop_na(n_fix_reads[n_fix_reads$seq_N_fixer_type=="I",])
t1_only_sum<-aggregate(genus_prop_abund ~ Sample, data = t1_only, FUN = sum)
t1_only_sum<-rename(t1_only_sum, "seq_pot_nfix_t1"="genus_prop_abund")

t2_t3<-drop_na(n_fix_reads[(n_fix_reads$seq_N_fixer_type=="II"|
                              n_fix_reads$seq_N_fixer_type=="III"),])
t2_t3_sum<-aggregate(genus_prop_abund ~ Sample, data = t2_t3, FUN = sum)
t2_t3_sum<-rename(t2_t3_sum, "seq_pot_nfix_t2t3"="genus_prop_abund")

all_nfix<-drop_na(n_fix_reads[n_fix_reads$seq_N_fixer=="Yes",])
all_nfix_sum<-aggregate(genus_prop_abund ~ Sample, data = all_nfix, FUN = sum)
all_nfix_sum<-rename(all_nfix_sum, "seq_pot_all_nfix"="genus_prop_abund")

merge_table<-merge(t1_only_sum, t2_t3_sum, by="Sample", all=TRUE)
merge_table2<-merge(merge_table, all_nfix_sum)

excel_output_path <- paste0(dir,"/seq_nfix_sum.xlsx")
write_xlsx(merge_table2, path = excel_output_path)
#### Summ. Stats - Total heterocysts /L for DATE and SWMP ID ####
# DATE
summary_stats_list <- list()

# Get the names of the numerical variables you want to include
numerical_cols <- names(cyano_sub[c(10,11,13)])[
  !names(cyano_sub[c(10,11,13)]) %in% "Date"]

cyano_sub$Date <- as.factor(cyano_sub$Date)

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  summ_stats <- (cyano_sub%>% 
                   group_by(Date) %>%
                   get_summary_stats(variable_name, type = "mean_sd"))
  
  stats_table <- as.data.frame((summ_stats))
  stats_table$Cat_Variable <- rownames(stats_table)
  stats_table <- stats_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  stats_table$Dependent_Variable <- variable_name
  summary_stats_list[[variable_name]] <- stats_table
}

final_summ_stats <- (bind_rows(summary_stats_list))[c(2,4:7)]

dir.create(file.path(paste0(dir,"/stats"), "/n_lim_cyano_stats"), showWarnings = FALSE)
excel_output_path <- paste0(dir,"/stats/n_lim_cyano_stats/het_DATE_summ_stats.xlsx")
write_xlsx(final_summ_stats, path = excel_output_path)

# SWMP ID
summary_stats_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(cyano_sub[c(10,11,13)])[
  !names(cyano_sub[c(10,11,13)]) %in% "IDL"]

cyano_sub$Date <- as.factor(cyano_sub$Date)

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  summ_stats <- (cyano_sub %>% 
                   group_by(IDL) %>%
                   get_summary_stats(variable_name, type = "mean_sd"))
  
  stats_table <- as.data.frame((summ_stats))
  stats_table$Cat_Variable <- rownames(stats_table)
  stats_table <- stats_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  stats_table$Dependent_Variable <- variable_name
  summary_stats_list[[variable_name]] <- stats_table
}

final_summ_stats <- (bind_rows(summary_stats_list))[c(2,4:7)]

excel_output_path <- paste0(dir,"/stats/n_lim_cyano_stats/het_SWMPID_summ_stats.xlsx")
write_xlsx(final_summ_stats, path = excel_output_path)

#### ANOVA - Heterocysts /L + potential N-fixers (sequencing) for DATE ####
dataset<-cyano_data[c(1:4,7:11,12,20,46,58:60,105)]

dataset<-merge(dataset, merge_table2, by="Sample", all=TRUE)

dataset2<-dataset %>%
  rename("cyano_library_size"="Cyanobacteria")%>%
  mutate(seq_pot_nfix_t1 = case_when(is.na(seq_pot_nfix_t1) & is.na(cyano_library_size) ~ NA,
                                  is.na(seq_pot_nfix_t1) & cyano_library_size >0  ~ 0,
                                  seq_pot_nfix_t1 > 0 & cyano_library_size >0  ~ seq_pot_nfix_t1))%>%
  mutate(seq_pot_nfix_t2t3 = case_when(is.na(seq_pot_nfix_t2t3) & is.na(cyano_library_size) ~ NA,
                                     is.na(seq_pot_nfix_t2t3) & cyano_library_size >0  ~ 0,
                                     seq_pot_nfix_t2t3 > 0 & cyano_library_size >0  ~ seq_pot_nfix_t2t3))%>%
  mutate(seq_pot_all_nfix = case_when(is.na(seq_pot_all_nfix) & is.na(cyano_library_size) ~ NA,
                                     is.na(seq_pot_all_nfix) & cyano_library_size >0  ~ 0,
                                     seq_pot_all_nfix > 0 & cyano_library_size >0  ~ seq_pot_all_nfix))

dataset2$logseq_pot_nfix_t1<- log10(dataset2$seq_pot_nfix_t1+0.001)
dataset2$logseq_pot_nfix_t2t3 <- log10(dataset2$seq_pot_nfix_t2t3+0.001)
dataset2$logseq_pot_all_nfix <- log10(dataset2$seq_pot_all_nfix+0.001)

# Initialize an empty list to store ANOVA results
anova_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(dataset2[c(13:14,16:22)])[
  !names(dataset2[c(13:14,16:22)])  %in% "Date"]

dataset2$Date <- as.factor(dataset2$Date)
# set Date label order
date_level<-c("June 23rd","July 20th","Aug 3rd",
              "Aug 23rd", "Aug 31st",  "Sept 27th")

# Loop through each variable
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

excel_output_path <- paste0(dir,"/stats/n_lim_cyano_stats/n_lim_DATE_anova.xlsx")
write_xlsx(final_anova_df_2, path = excel_output_path)

#### T-test - all quantitative data for Heterocyst pres/abs  ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the numerical variables you want to test
numerical_cols <- names(cyano_data[c(4:10,21:37,47,48,53:56,61:98,100:103)])[
  !names(cyano_data[c(4:10,21:37,53:56,61:98,100:103)]) %in% "Het_Pres"]

# Loop through each numerical variable
for (variable_name in numerical_cols) {
  # Create the formula for T-test
  formula_t.test <- as.formula(paste(variable_name, "~ Het_Pres"))
  
  # Perform T-test
  t.test<- t.test(formula_t.test, data = cyano_data)
  
  t.test_table <- as.data.frame((t.test[1:4]))
  t.test_table$Cat_Variable <- "Het_Pres"
  t.test_table <- t.test_table %>%
    select(Cat_Variable, everything()) # Reorder columns
  
  t.test_table$Dependent_Variable <- variable_name
  t.test_results_list[[variable_name]] <- t.test_table
}

final_t.test_df <- bind_rows(t.test_results_list)

excel_output_path <- paste0(dir,"/stats/n_lim_cyano_stats/het_PRES_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)
#### T-test - Heterocysts / L and all categorical variables ####
# Initialize an empty list to store T-test results
t.test_results_list <- list()

# Get the names of the categorical variables you want to test
categorical_cols <- names(cyano_sub[c(7:9)])[
  !names(cyano_sub[c(7:9)]) %in% "logHet"]

cyano_sub$Storm_Base<-as.factor(cyano_sub$Storm_Base)
cyano_sub$Dredge_Cat<-as.factor(cyano_sub$Dredge_Cat)
cyano_sub$Pre_2003<-as.factor(cyano_sub$Pre_2003)

# Loop through each categorical variable
for (variable_name in categorical_cols) {
  # Create the formula for T-test
  formula_t.test <- as.formula(paste("logHet ~", variable_name))
  
  # Perform T-test
  t.test<- t.test(formula_t.test, data = cyano_sub)
  
  t.test_table <- as.data.frame((t.test[1:4]))
  t.test_table$Num_Variable <- "logHet"
  t.test_table <- t.test_table %>%
    select(Num_Variable, everything()) # Reorder columns
  
  t.test_table$Dependent_Variable <- variable_name
  t.test_results_list[[variable_name]] <- t.test_table
}

final_t.test_df <- bind_rows(t.test_results_list)

excel_output_path <- paste0(dir,"/stats/n_lim_cyano_stats/het_allCAT_t.test.xlsx")
write_xlsx(final_t.test_df, path = excel_output_path)
#### Pearson - heterocysts + N-fix x WQ and abundance data ####
# Subset to variables of interest 
dataset<-cyano_data[c(1:4,7:12,20:37,46:48,50,52:56,58:59,61:105)]
dataset<-merge(dataset, merge_table2, by="Sample", all=TRUE)

dataset2<-dataset %>%
  rename("cyano_library_size"="Cyanobacteria")%>%
  mutate(seq_pot_nfix_t1 = case_when(is.na(seq_pot_nfix_t1) & is.na(cyano_library_size) ~ NA,
                                     is.na(seq_pot_nfix_t1) & cyano_library_size >0  ~ 0,
                                     seq_pot_nfix_t1 > 0 & cyano_library_size >0  ~ seq_pot_nfix_t1))%>%
  mutate(seq_pot_nfix_t2t3 = case_when(is.na(seq_pot_nfix_t2t3) & is.na(cyano_library_size) ~ NA,
                                       is.na(seq_pot_nfix_t2t3) & cyano_library_size >0  ~ 0,
                                       seq_pot_nfix_t2t3 > 0 & cyano_library_size >0  ~ seq_pot_nfix_t2t3))%>%
  mutate(seq_pot_all_nfix = case_when(is.na(seq_pot_all_nfix) & is.na(cyano_library_size) ~ NA,
                                      is.na(seq_pot_all_nfix) & cyano_library_size >0  ~ 0,
                                      seq_pot_all_nfix > 0 & cyano_library_size >0  ~ seq_pot_all_nfix))

dataset2$logseq_pot_nfix_t1<- log10(dataset2$seq_pot_nfix_t1+0.001)
dataset2$logseq_pot_nfix_t2t3 <- log10(dataset2$seq_pot_nfix_t2t3+0.001)
dataset2$logseq_pot_all_nfix <- log10(dataset2$seq_pot_all_nfix+0.001)

# Get the names of the n-lim and WQ variables being tested
nlim_cols <- colnames(dataset2)[c(38,39,84:90)]
nlim_cols <- nlim_cols[!nlim_cols %in% "Sample"]

wq_cols <- colnames(dataset2)[c(4:5,7,12:28,30,31,34:37,40:77,79:82)]
wq_cols <- wq_cols[!wq_cols %in% "Sample"]

pearson_results_list <- list()

for (nlim_variable in nlim_cols) {
  pearson_table <- data.frame(
    wq_variable=numeric(),
    cor=numeric(),
    p=numeric()
  )
  for (wq_variable in wq_cols) {
    # Perform Pearson correlation
    pearson <- cor_test(nlim_variable,wq_variable, data=dataset2, method="pearson",
                        alternative="two.sided")
    pearson<-as.data.frame(pearson)
    pearson_table[nrow(pearson_table) + 1,] <- list(pearson$var2, pearson$cor, pearson$p)    
    pearson_results_list[[nlim_variable]] <- pearson_table
  }
}

excel_output_path <- paste0(dir,"/stats/n_lim_cyano_stats/nlim_PEARSON.xlsx")
write_xlsx(pearson_results_list, path = excel_output_path)
