
library(tidyverse)
library(Biostrings)

rstudioapi::getActiveDocumentContext()$path %>% 
  dirname() %>% setwd()
# setwd("C:/Users/haohe/Desktop")


## load vector sequence
lentiCRISPRv2 <- "./CRISPRv2.txt" %>% 
  read_delim(delim="/t", col_names=F) %>% 
  pull()
# replace this sequence with sgRNA sequence
lentiCRISPRv2_replacement <- "gagacggttgtaaatgagcacacaaaatacacatgctaaaatattatattctatgacctttataaaatcaaccaaaatcttctttttaataactttagtatcaataattagaatttttatgttcctttttgcaaacttttaataaaaatgagcaaaataaaaaaacgctagttttagtaactcgcgttgttttcttcacctttaataatagctactccaccacttgttcctaagcggtcagctcctgcttcaatcattttttgagcatcttcaaatgttctaactccaccagctgctttaactaaagcattgtctttaacaactgacttcattagtttaacatcttcaaatgttgcacctgattttgaaaatcctgttgatgttttaacaaattctaatccagcttcaacagctatttcacaagctttcatgatttcttcttttgttaataaacaattttccataatacatttaacaacatgtgatccagctgctttttttacagctttcatgtcttctaaaactaattcataatttttgtcttttaatgcaccaatatttaataccatatcaatttctgttgcaccatctttaattgcttcagaaacttcgaatgcttttgtagctgttgtgcatgcacctagaggaaaacctacaacatttgttattcctacatttgtgccttttaataattctttacaatagcttgttcaatatgaattaacacaaactgttgcaaaatcaaattcaattgcttcatcacataattgtttaatttcagctttcgtagcatcttgttttaataatgtgtgatctatatatttgtttagtttcattttttctcctatatattcatttttaattttaattctttaataatttcgtctactttaactttagcgttttgaacagattcaccaacacctataaaataaatttttagtttaggttcagttccacttgggcgaacagcaaatcatgacttatcttctaaataaaattttagtaagtcttgtcctggcatattatacattccatcgatgtagtcttcaacattaacaactttaagtccagcaatttgagttaagggtgttgctctcaatgatttcattaatggttcaatttttaatttcttttcttctggtttaaaattcaagtttaaagtgaaagtgtaatatgcacccatttctttaaataaatcttctaaatagtctactaatgttttattttgttttttataaaatcaagcagcctctgctattaatatagaagcttgtattccatctttatctctagctgagtcatcaattacatatccataactttcttcataagcaaaaacaaaatttaatccgttatcttcttctttagcaatttctctacccattcatttaaatccagttaaagtttttacaatattaactccatatttttcatgagcgattctatcacccaaatcacttgttacaaaacttgaatatagagccggattttttggaatgctatttaagcgttttagatttgataattttcaatcaattaaaattggtcctgtttgatttccatctaatcttacaaaatgaccatcatgttttattgccattccaaatctgtcagcatctgggtcattcataataataatatctgcatcatgtttaataccatattcaagcggtatttttcatgcaggatcaaattctggatttggatttacaacatttttaaatgtttcatcttcaaatgcatgctcttcaacctcaataacgttatatcctgattcacgtaatatttttggggtaaatttagttcctgttccattaactgcgctaaaaataatttttaaatcttttttagcttcttgctcttttttgtacgtctct"


if(T) {
  ## load sgRNA design table
  # for CRISPR KO, download from https://portals.broadinstitute.org/gppx/crispick/public
  crisprPick_file <- "C:/Users/haohe/Downloads/b093d1f6-0bc8-4ba2-aab7-cc899e149e90-sgrna-designs.txt"
  crispr <- read_delim(crisprPick_file) %>%
    dplyr::select(symbol=`Target Gene Symbol`,
                  sgRNA=`sgRNA Sequence`)
  
  crispr <- readxl::read_excel(
    file.path("C:/Users/haohe/Desktop/2025-03-31_rra_selected_hits.xlsx")) %>% 
    arrange(group) %>% 
    mutate(symbol=str_c(group, symbol, sep="_")) %>% 
    dplyr::select(symbol=symbol,
                  sgRNA=sequence)
  
  ## or manually defined
  crispr <- tibble(
    symbol=c("PSEN1", "PSEN2"),
    sgRNA=c("CCCTGAGCCATTATCTAATGGAC","CCCTAATGTCGGCTGAGAGCCCC")
  )
}


## get primer and map
crispr_res <- crispr %>% 
  mutate(rep=row_number(), .by=c(symbol), .after=symbol) %>% 
  mutate(sgRNA_revCom=DNAStringSet(sgRNA) %>% 
           reverseComplement() %>% 
           as.character()) %>% 
  mutate(
    f_name=str_glue("sg{symbol}#{rep}_F"),
    f_primer=str_c("caccg", sgRNA),
    r_name=str_glue("sg{symbol}#{rep}_R"),
    r_primer=str_c("aaac", sgRNA_revCom, "c")) %>% 
  mutate(map=map_chr(sgRNA, \(x) str_replace(
    lentiCRISPRv2, lentiCRISPRv2_replacement, x))) %>% 
  print()


## write
write_csv(crispr_res, str_glue("{Sys.Date()}_crispr.csv"))
if(!dir.exists("map")) dir.create("map")
for(i in seq_len(nrow(crispr_res))) {
  Biostrings::writeXStringSet(
    DNAStringSet(crispr_res[[i, "map"]]),
    str_glue('map/P{i+253}_sg_{crispr_res[[i, "symbol"]]}#{crispr_res[[i, "rep"]]}.fasta')
  )
}



