
# processing gtf to bed
library(tidyverse)
library(GenomicRanges)

setwd("D:/Index/hs/v49")
# setwd("D:/Index/mm/vM38")

gtf <- list.files(pattern=".gtf.gz") %>% .[1] %>% 
  rtracklayer::import()


gtf$source %>% table()
gtf$level %>% table()
gtf$type %>% table()
gtf$gene_type %>% table()
gtf$tag %>% table()


genes <- 
  gtf[
    gtf$type == "gene" & 
      gtf$level %in% c(1, 2) &
      # gtf$source == "HAVANA" &
      str_detect(gtf$gene_name,'^MT[t]-', negate=T),] #%>%  .[1:100,] %>% view()

# tss <- promoters(genes, 0, 0)

table(genes$gene_type)
table(genes$source)



# save to UCSC style bed
genes[genes$gene_type == "protein_coding", "gene_name"] %>% 
  rtracklayer::export("genes_protein_coding.bed")

genes[genes$gene_type == "lncRNA", "gene_name"] %>% 
  rtracklayer::export("genes_lncRNA.bed")

genes[!(genes$gene_type %in% c("protein_coding", "lncRNA")), "gene_name"] %>% 
  rtracklayer::export("genes_others.bed")
