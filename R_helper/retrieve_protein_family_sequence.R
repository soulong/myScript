
library(org.Hs.eg.db)
library(UniProt.ws)
library(tidyverse)
library(readxl)
library(writexl)
library(Biostrings)
library(furrr)
plan(multisession, workers=4)

# set data directory
dataset_dir <- "/media/hao/Data/Project_TF_PPI/2024-09-12_coactivators"

## retrieve items
library(biomaRt)
# listDatasets(useMart("ENSEMBL_MART_ENSEMBL")) %>% view()
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
# listFilters(mart) %>% view()
# listAttributes(mart) %>% view()
gene <- getBM(
  attributes=c("hgnc_symbol","ensembl_gene_id","hgnc_id","entrezgene_id","gene_biotype"),
  filters="go",
  values=c("GO:0003713"),
  mart=mart, uniqueRows=T, useCache=F) %>% as_tibble()
gene <- gene %>% 
  filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>% 
  distinct(hgnc_symbol, .keep_all=T) %>% 
  arrange(hgnc_symbol) %>% 
  print()

## get transcript id
transcript <- getBM(
  attributes=c("ensembl_gene_id","ensembl_transcript_id","uniprotswissprot", 
               "transcript_gencode_basic","transcript_is_canonical",
               "transcript_appris","transcript_biotype"),
  filters="ensembl_gene_id",
  values=unique(gene$ensembl_gene_id),
  mart=mart, uniqueRows=T, useCache=F) %>% as_tibble()
transcript <- transcript %>%
  filter(transcript_biotype == "protein_coding", transcript_is_canonical==1) %>%
  print()

## retrieve cdna & aa sequence
seq <- getSequence(
  id=transcript$uniprotswissprot,
  type="uniprotswissprot",
  seqType="coding",
  mart=mart, useCache=F) %>% as_tibble() %>% 
  dplyr::rename(cdna=coding)

  

## retrieve aa sequence and annotation with uniprot
up <- UniProt.ws(taxId=9606)
# keytypes(up)
# columns(up)
aa_up <- UniProt.ws::select(
  up, 
  keys=transcript$uniprotswissprot, 
  columns=columns(uniprot),
  # columns=c("id","protein_families","sequence","reviewed", 
  #           # str_subset(columns(up), "^ft_"),
  #           "ft_signal","ft_chain","ft_region","ft_motif","ft_binding","ft_topo_dom","ft_transmem","ft_dna_bind"
  #           ), 
  keytype="UniProtKB") %>% as_tibble() %>% 
  filter(Reviewed=="reviewed") %>% 
  glimpse()
# filter annotation
aa_up_filter <- aa_up %>% 
  dplyr::select(uniprotswissprot=From, uniprot_aa_sequence=Sequence)


## merge cdna and aa
library(Biostrings)
res <- left_join(aa_up_filter, seq) %>% 
  mutate(cdna_aa_sequence=map_chr(cdna, \(x) DNAString(x) %>% translate() %>% 
                                    as.character() %>% str_replace("[*]","")))

res2 <- res %>% 
  mutate(if_aa_same = ifelse(uniprot_aa_sequence==cdna_aa_sequence, 1, 0)) %>% 
  arrange(uniprotswissprot, desc(if_aa_same)) %>% 
  slice_head(n=1, by=uniprotswissprot) %>% 
  mutate(uniprot_aa_length=str_length(uniprot_aa_sequence), 
         .after=uniprot_aa_sequence)
# add symbol
res3 <- getBM(
  attributes=c("hgnc_symbol","uniprotswissprot"),
  filters="uniprotswissprot",
  values=unique(res2$uniprotswissprot),
  mart=mart, uniqueRows=T, useCache=F) %>% as_tibble() %>% 
  left_join(res2, .) %>% 
  arrange(hgnc_symbol) %>%
  relocate(hgnc_symbol, .before=1)



## save data
write_xlsx(res3, file.path(dataset_dir, "coactivators.xlsx"))


# write fasta
fasta <- data[, c("fasta_id","protein_sequence_dRTK")] %>% 
  deframe() %>% AAStringSet() %>% print()
writeXStringSet(fasta, file.path(dataset_dir, "family_dRTK.fasta"))


