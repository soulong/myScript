
# devtools::install_github("stjude/ChIPseqSpikeInFree")
library(ChIPseqSpikeInFree)
library(tidyverse)


# bam_dir <- "F:/workspace/IMR/result/03_alignment"
bam_dir <- "F:/workspace/MEF/result/03_alignment"

# gsize_file <- "F:/software/R_scripts/GRCh38.primary_assembly.genome.fa.fai" # hs
gsize_file <- "F:/software/R_scripts/GRCm39.primary_assembly.genome.fa.fai" # mm


setwd(bam_dir)
sample_sheet <- "../../samplesheet_CutTag.csv"
sample_info <- read_csv(sample_sheet) %>% #glimpse()
  mutate(ID=str_c(sample, ".filtered.bam")) %>% 
  glimpse()

save_prefix="SpikeFree"
meta_file <- str_c(save_prefix, "_metaFile.txt")
tibble(ID=sample_info$ID, ANTIBODY="H3K27me3", GROUP=sample_info$group) %>% 
  write_tsv(meta_file)

# run
ChIPseqSpikeInFree(bamFiles=sample_info$ID, 
                   chromFile=gsize_file,
                   metaFile=meta_file,
                   prefix=save_prefix,
                   cutoff_QC=1.2,
                   maxLastTurn=0.99,
                   ncores=6)


# # transform sf for bigwig/bedgraph normalization
# statistics <- "../statistics_CPM.csv" %>% 
#   read_csv() %>% #glimpse()
#   # transmute(sample, filtered_reads) %>% 
#   left_join(sample_info) %>%
#   print()
# 
# scale_factor <- str_c(save_prefix, "_SF.txt") %>% 
#   read_table() %>%
#   left_join(statistics) %>% 
#   mutate(scale_factor=1000000 / (filtered_reads * SF)) %>% 
#   glimpse()
# 
# # save
# scale_factor %>% 
#   write_csv(str_c(save_prefix, ".csv"))
# 
# scale_factor %>% 
#   dplyr::select(sample, scale_factor=scale_factor, SF) %>% 
#   write_csv(str_c(save_prefix, "_scale_factors.csv"))


