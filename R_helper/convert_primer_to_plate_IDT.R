
library(tidyverse)
library(readxl)
library(writexl)

wd <- 'Z:\\Data\\Others\\2025-08-15_TF_cloning'
setwd(wd)
getwd()


# load primer
f <- '2025-08-29_cloning_primer.csv'
primer <- read_csv(f) %>% 
  mutate(id=row_number()) %>% print()


# set plate
plate_nrow <- 8
plate_ncol <- 12
plate <- tibble(well=map2_chr(
  rep(LETTERS[1:plate_nrow], plate_ncol), 
  rep(1:plate_ncol, each = plate_nrow), # %>% str_pad(2, pad = "0"), 
  ~ str_c(.x, .y))) %>% 
  mutate(id=row_number()) %>% print()


# set position
plate_primer <- merge(plate, primer, by='id', all = T) %>% 
  as_tibble() %>% 
  filter(!is.na(name)) %>% glimpse()


# convert to long shape
plate_primer_long <- plate_primer %>% 
  dplyr::select(!c(id, primer_r)) %>% 
  dplyr::rename(`_f`=primer_f, `_r`=primer_r_revcomp) %>% 
  pivot_longer(`_f`:`_r`, names_to="order",values_to='sequence') %>% 
  unite('name', name:order, sep='', remove=T) %>% 
  print()
  

# save
write_csv(plate_primer_long, str_replace(f,'.csv','_plate.csv'))