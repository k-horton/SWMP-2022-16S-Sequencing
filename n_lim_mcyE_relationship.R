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
cyano_data<-read_excel(paste0(dir,"/wq_cyano_ASV_abundance.xlsx"))

nfix_data<-read_excel(paste0(dir,"/seq_nfix_sum.xlsx"))
tox_data<-read_excel(paste0(dir,"/seq_tox_prod_sum.xlsx"))

#### Prepare dataframes####
# re-name Other Bacteria column
cyano_data<-rename(cyano_data, Other_Bacteria = `Other Bacteria`)

# set date as a factor and order chronologically
cyano_data$Date<-factor(cyano_data$Date,
                        levels=c("June 23rd", "July 20th","Aug 3rd",
                                 "Aug 23rd", "Aug 31st", "Sept 27th"),
                        labels = c("June 23rd",    "July 20th",  "Aug 3rd",
                                   "Aug 23rd", "Aug 31st", "Sept 27th") )

genus_data$Date<-factor(genus_data$Date,
                        levels=c("June 23rd", "July 20th","Aug 3rd",
                                 "Aug 23rd", "Aug 31st", "Sept 27th"),
                        labels = c("June 23rd",    "July 20th",  "Aug 3rd",
                                   "Aug 23rd", "Aug 31st", "Sept 27th") )
#### Pearson - heterocysts + N-fix x toxin data ####
# Subset to variables of interest 
dataset<-cyano_data[c(1,4,50,52,55,58,59,99,102,104,105)]
dataset1<-merge(dataset, nfix_data, by="Sample", all=TRUE)
dataset2<-merge(dataset1, tox_data, by="Sample", all=TRUE)

dataset3<-dataset2 %>%
  rename("cyano_library_size"="Cyanobacteria")%>%
  mutate(seq_pot_nfix_t1 = case_when(is.na(seq_pot_nfix_t1) & is.na(cyano_library_size) ~ NA,
                                     is.na(seq_pot_nfix_t1) & cyano_library_size >0  ~ 0,
                                     seq_pot_nfix_t1 > 0 & cyano_library_size >0  ~ seq_pot_nfix_t1))%>%
  mutate(seq_pot_nfix_t2t3 = case_when(is.na(seq_pot_nfix_t2t3) & is.na(cyano_library_size) ~ NA,
                                       is.na(seq_pot_nfix_t2t3) & cyano_library_size >0  ~ 0,
                                       seq_pot_nfix_t2t3 > 0 & cyano_library_size >0  ~ seq_pot_nfix_t2t3))%>%
  mutate(seq_pot_all_nfix = case_when(is.na(seq_pot_all_nfix) & is.na(cyano_library_size) ~ NA,
                                      is.na(seq_pot_all_nfix) & cyano_library_size >0  ~ 0,
                                      seq_pot_all_nfix > 0 & cyano_library_size >0  ~ seq_pot_all_nfix))%>%
  mutate(seq_pot_mcyE = case_when(is.na(seq_pot_mcyE) & is.na(cyano_library_size) ~ NA,
                                  is.na(seq_pot_mcyE) & cyano_library_size >0  ~ 0,
                                  seq_pot_mcyE > 0 & cyano_library_size >0  ~ seq_pot_mcyE))%>%
  mutate(seq_pot_oth_tox = case_when(is.na(seq_pot_oth_tox) & is.na(cyano_library_size) ~ NA,
                                     is.na(seq_pot_oth_tox) & cyano_library_size >0  ~ 0,
                                     seq_pot_oth_tox > 0 & cyano_library_size >0  ~ seq_pot_oth_tox))%>%
  mutate(seq_pot_all_tox = case_when(is.na(seq_pot_all_tox) & is.na(cyano_library_size) ~ NA,
                                     is.na(seq_pot_all_tox) & cyano_library_size >0  ~ 0,
                                     seq_pot_all_tox > 0 & cyano_library_size >0  ~ seq_pot_all_tox))


dataset3$logseq_pot_nfix_t1<- log10(dataset3$seq_pot_nfix_t1+0.001)
dataset3$logseq_pot_nfix_t2t3 <- log10(dataset3$seq_pot_nfix_t2t3+0.001)
dataset3$logseq_pot_all_nfix <- log10(dataset3$seq_pot_all_nfix+0.001)

dataset3$logseq_pot_mcyE <- log10(dataset3$seq_pot_mcyE+0.001)
dataset3$logseq_pot_oth_tox <- log10(dataset3$seq_pot_oth_tox+0.001)
dataset3$logseq_pot_all_tox <- log10(dataset3$seq_pot_all_tox+0.001)

# Get the names of the n-lim and WQ variables being tested
nlim_cols <- colnames(dataset3)[c(5:7,9,11:14,18:20)]
nlim_cols <- nlim_cols[!nlim_cols %in% "Sample"]

tox_cols <- colnames(dataset3)[c(1,2,8,10,15:17,21:23)]
tox_cols <- tox_cols[!tox_cols %in% "Sample"]

pearson_results_list <- list()

for (nlim_variable in nlim_cols) {
  pearson_table <- data.frame(
    wq_variable=numeric(),
    cor=numeric(),
    p=numeric()
  )
  for (tox_variable in tox_cols) {
    # Perform Pearson correlation
    pearson <- cor_test(nlim_variable,tox_variable, data=dataset3, method="pearson",
                        alternative="two.sided")
    pearson<-as.data.frame(pearson)
    pearson_table[nrow(pearson_table) + 1,] <- list(pearson$var2, pearson$cor, pearson$p)    
    pearson_results_list[[nlim_variable]] <- pearson_table
  }
}

excel_output_path <- paste0(dir,"/stats//nlim_tox_PEARSON.xlsx")
write_xlsx(pearson_results_list, path = excel_output_path)
