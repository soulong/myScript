
if(F) {
  rstudioapi::getActiveDocumentContext()$path |>
    dirname() |> setwd()
} else {
  setwd("/media/hao/Data/Project_kinase_substrate_design/2025-05-12_ALK_sub")
}

library(tidyverse)
library(readxl)
library(writexl)
# devtools::install_github("greg-botwin/reversetranslate")
library(reversetranslate)
library(Biostrings)

## * design primer (codon optimize)  -------------
aa_seq <- "2025-06-16_result_merged_top2_selected.xlsx" %>% 
  # readxl::read_excel(skip=1, col_names=F) %>% 
  readxl::read_excel() %>% 
  dplyr::select(2:3) %>% 
  set_names(c("id","aa")) %>% 
  filter(!is.na(id)) %>% 
  mutate(aa=str_to_upper(aa)) %>% 
  print()


# rev translate to forward
# codon usage table: reversetranslate::hsapien_tbl
set.seed(42)
aa_seq <- aa_seq %>% 
  # mutate(f=seq_rev_translate(as_aa(aa))) %>% 
  mutate(f=map_chr(aa, \(x) suppressMessages(
    reverse_translate(
      amino_acid_seq=x,
      codon_tbl=hsapien_tbl,
      limit=0, model="proportional")))) %>%
  print()

# get reverse sequence
aa_seq <- aa_seq %>% 
  mutate(r=reverseComplement(DNAStringSet(f)) %>% 
           as.character()) %>% 
  print()

# get feature
aa_seq <- aa_seq %>% 
  mutate(gc=letterFrequency(DNAStringSet(f), letters="CG", as.prob=T) %>% 
           magrittr::multiply_by(100) %>% round() %>% as.vector()) %>% 
  print()

# add cut site sequence
aa_seq <- aa_seq %>% 
  mutate(f_primer=str_c("ctaga", f, "g"),
         r_primer=str_c("gatcc", r, "t")) %>% 
  # primer name
  mutate(f_name=((row_number()-1)*1 + 291), .before=f_primer) %>% 
  mutate(f_name=str_glue("#{f_name}_f")) %>% 
  mutate(r_name=((row_number()-1)*1 + 291), .before=r_primer) %>% 
  mutate(r_name=str_glue("#{r_name}_r")) %>% 
  print()


# write
writexl::write_xlsx(aa_seq, str_glue(
  "{Sys.Date()}_aa_seq.xlsx"))



