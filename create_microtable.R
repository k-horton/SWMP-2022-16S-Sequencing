# The qiime2meco() function from the file2meco library is designed to create the microtable object using files from QIIME2
# Load in file2meco library
library(file2meco)

# Assign current working directory to 'dir'
dir<-getwd()

# alternatively, manually set the working directory:
# dir<- "C:/Path/To/Directory"

# Define the path to each of the files needed to create the microtable object. Replace [Path to file] with the actual path
abund_file_path <- paste0(dir, "/merge_concat.seqs_count_table.qza")
sample_file_path <- paste0(dir, "/metadata_decontam.xlsx")
taxonomy_file_path <- paste0(dir, "/merge_concat.silva_taxonomy.qza")
tree_data <- paste0(dir, "/merge_concat.rooted_tree.qza")
rep_data <- paste0(dir, "/merge_concat.representative_seqs.qza")

# construct microtable object
table1 <- qiime2meco(abund_file_path, sample_table = sample_file_path, 
                     taxonomy_table = taxonomy_file_path, phylo_tree = tree_data, 
                     rep_fasta = rep_data, auto_tidy = TRUE)
# 2 samples with 0 abundance removed from otu_table

table1
#microtable-class object:
  #sample_table have 71 rows and 45 columns
  #otu_table have 9040 rows and 71 columns
  #tax_table have 9040 rows and 7 columns

# Verify all samples were imported correctly
table1$sample_names()
# [1] "F1_SWMP1_S1"          "F1_SWMP2_S2"          "F1_SWMP3_S3"          "F1_SWMP4_S4"         
# [5] "F1_SWMP5_S5"          "F1_SWMP6_S6"          "F1_SWMP7_S7"          "F1_SWMP8_S8"         
# [9] "F1_SWMP9_S9"          "F1_SWMP10_S10"        "F1_SWMP11_S11"        "F1_SWMP12_S12"       
# [13] "F2_SWMP1_S13"         "F2_SWMP2_S14"         "F2_SWMP3_S15"         "F2_SWMP4_S16"        
# [17] "F2_SWMP5_S17"         "F2_SWMP6_S18"         "F2_SWMP7_S19"         "F2_SWMP8_S20"        
# [21] "F2_SWMP9_S21"         "F2_SWMP10_S22"        "F2_SWMP11_S23"        "F2_SWMP12_S24"       
# [25] "F3_SWMP1_S25"         "F3_SWMP2_S26"         "F3_SWMP3_S27"         "F3_SWMP4_S28"        
# [29] "F3_SWMP5_S29"         "F3_SWMP6_S30"         "F3_SWMP8_S32"         "F3_SWMP9_S33"        
# [33] "F3_SWMP10_S34"        "F3_SWMP11_S35"        "F3_SWMP12_S36"        "F4_SWMP1_S37"        
# [37] "F4_SWMP2_S38"         "F4_SWMP3_S39"         "F4_SWMP4_S40"         "F4_SWMP5_S41"        
# [41] "F4_SWMP6_S42"         "F4_SWMP7_S43"         "F4_SWMP9_S45"         "F4_SWMP10_S46"       
# [45] "F4_SWMP11_S47"        "F4_SWMP12_S48"        "F5_SWMP1_S49"         "F5_SWMP2_S50"        
# [49] "F5_SWMP3_S51"         "F5_SWMP4_S52"         "F5_SWMP5_S53"         "F5_SWMP6_S54"        
# [53] "F5_SWMP7_S55"         "F5_SWMP8_S56"         "F5_SWMP9_S57"         "F5_SWMP10_S58"       
# [57] "F5_SWMP11_S59"        "F5_SWMP12_S60"        "F6_SWMP1_S61"         "F6_SWMP2_S62"        
# [61] "F6_SWMP3_S63"         "F6_SWMP4_S64"         "F6_SWMP6_S65"         "F6_SWMP7_S66"        
# [65] "F6_SWMP8_S67"         "F6_SWMP9_S68"         "F6_SWMP10_S69"        "F6_SWMP11_S70"       
# [69] "F6_SWMP12_S71"        "POS_CTRL_S72"         "Negative_Control_S73"

# Remove pollution (mitochondrial and chloroplast sequences) from the microtable
table1$filter_pollution(taxa = c("mitochondria", "chloroplast"))
# Total 366 features are removed from tax_table ...
