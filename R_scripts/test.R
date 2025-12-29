
install.packages("writexl")
BiocManager::install("bambu")
BiocManager::install(c("rtracklayer","enrichplot","ChipSeeker","Rsamtools","Biostrings","GenomicFeatures","TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm39.knownGene","GSVA","GenomicAlignment","biomaRt","janitor"))

.libPaths()

library(tidyverse)
rstudioapi::getActiveDocumentContext()$path %>% 
  dirname() %>% setwd()
