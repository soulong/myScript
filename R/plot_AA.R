library(Biostrings)
library(tidyverse)

seq <- Biostrings::readAAStringSet("C:/Users/haohe/Desktop/Q16288.fasta") %>%
  as.character() %>%
  set_names(nm = NULL) %>%
  str_split("", simplify = T) %>%
  as.character()
  

uni <- unique(seq) %>%
  set_names(., nm = .)

mat <- map_df(uni, ~ str_detect(seq, .x) %>% as.integer) %>%
  t() %>% 
  as.matrix()

pheatmap::pheatmap(mat, scale = "none", border_color = "red3",
                   cluster_cols = F, cluster_rows = T, 
                   color = colorRampPalette(c("white", "red3"))(50), 
                   filename = "C:/Users/haohe/Desktop/NTRK3.pdf", width = 5, height = 3)




seq <- Biostrings::readAAStringSet("C:/Users/haohe/Desktop/Q9Y3Y2.fasta") %>%
  as.character() %>%
  set_names(nm = NULL) %>%
  str_split("", simplify = T) %>%
  as.character()


uni <- unique(seq) %>%
  set_names(., nm = .)

mat <- map_df(uni, ~ str_detect(seq, .x) %>% as.integer) %>%
  t() %>% 
  as.matrix()
colnames(mat) <- seq_len(dim(mat)[2])

pheatmap::pheatmap(mat, scale = "none", border_color = "red3",
                   cluster_cols = F, cluster_rows = T, 
                   color = colorRampPalette(c("white", "red3"))(50), 
                   fontsize_col = 2,
                   filename = "C:/Users/haohe/Desktop/CHTOP1.pdf", width = 5, height = 3)
