
library(UniProt.ws)
library(tidyverse)
library(readxl)
library(writexl)
library(Biostrings)
library(biomaRt)
# library(furrr)
# plan(multisession, workers=4)


# set data directory
wd <- 'Z:\\Data\\Others\\2025-08-15_TF_cloning'
setwd(wd)
getwd()


## input gene list ----------------------------
complex <- read_excel('tf_coactivtor_match_ORFome.xlsx', sheet='coactivator') %>% 
  glimpse() 
# count(complex, type)
part <- complex %>% 
  dplyr::filter(!(type %in% c('swi_snf_complex', 'mediator'))) %>% 
  pull(symbol) %>% unique() %>% print()

# nymc <- read_excel('MycN interaction mapping.xlsx') %>% 
#   pull(1) %>% last() %>% 
#   str_split(' ', simplify=T) %>% 
#   str_squish() %>% unique() %>% print()
# 
# input_list <- c(part, nymc) %>% unique() %>% print()
input_list <- part

# # goids to gene
# # listDatasets(useMart("ENSEMBL_MART_ENSEMBL")) %>% view()
# mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host="https://useast.ensembl.org/")
# # listFilters(mart) %>% view()
# # listAttributes(mart) %>% view()
# goid_to_gene <- function(goids) {
#   gene <- getBM(
#     attributes=c("hgnc_symbol","ensembl_gene_id","hgnc_id","entrezgene_id","gene_biotype"),
#     filters="go",
#     values=goids,
#     mart=mart, uniqueRows=T, useCache=F) %>% as_tibble()
#   gene <- gene %>% 
#     filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>% 
#     distinct(hgnc_symbol, .keep_all=T) %>% 
#     arrange(hgnc_symbol)
#   return(gene)
# }




## add ensembl annotation ----------------------------
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
# listAttributes(mart, page='feature_page') %>% view()
# listAttributes(mart, page='sequences') %>% view()
# listFilters(mart) %>% view()    hgnc_symbol entrezgene_id go_id
# searchAttributes(mart=mart, 'entrez|MANE')

## id anno
mart_anno_id <- getBM(
  attributes=c('hgnc_symbol','entrezgene_id','uniprotswissprot',
               "ensembl_gene_id","ensembl_transcript_id",
               'transcript_mane_select','transcript_is_canonical',
               'transcript_gencode_basic','transcript_primary_basic'),
  filters="hgnc_symbol", values=input_list,
  mart=mart, uniqueRows=T, useCache=F) %>% 
  as_tibble() %>% glimpse()
mart_anno_id <- mart_anno_id %>%
  dplyr::rename(
    symbol=hgnc_symbol, entrez_id=entrezgene_id, uniprot=uniprotswissprot) %>% 
  group_by(symbol, entrez_id) %>% 
  mutate(across(everything(), \(x) ifelse(x=='',NA,x))) %>% 
  arrange(uniprot, transcript_mane_select,transcript_is_canonical,  
          transcript_gencode_basic, transcript_primary_basic, 
          .by_group=TRUE) %>% #view()
  slice_head(n=1) %>% 
  dplyr::select(!c(transcript_mane_select:transcript_primary_basic)) %>% 
  ungroup() %>% glimpse()

## sequence anno
mart_anno_cds <- getBM(
  attributes=c('ensembl_transcript_id','coding'),
  filters="ensembl_transcript_id", 
  values=mart_anno_id$ensembl_transcript_id,
  mart=mart, uniqueRows=T, useCache=F) %>% 
  as_tibble() %>% glimpse()
mart_anno_cds <- mart_anno_cds %>% 
  dplyr::rename(cds=coding) %>% 
  dplyr::mutate(cds=str_sub(cds, end=-4)) %>% # remove stop codon
  mutate(protein=translate(DNAStringSet(cds)) %>% as.character %>% 
           str_replace_all('[*]','')) %>% # remove stop AA *
  mutate(cds_length=str_length(cds), .before=cds) %>% glimpse()

## merge id and sequence
df <- left_join(mart_anno_id, mart_anno_cds) %>% 
  filter(cds_length < 10000) %>% 
  glimpse()

write_csv(df, str_glue('{Sys.Date()}_cloning.csv'))

# df <- read_csv('2025-08-18_cloning.csv')
ggplot(df, aes(x=cds_length)) +
  geom_histogram(binwidth = 100) +
  theme_bw() +
  labs(x='cDNA length', y='Count')



## add uniprot annotation ----------------------------
up <- UniProt.ws(taxId=9606)
# keytypes(up)  GeneID UniProtKB Gene_Name
# columns(up)
up_anno <- UniProt.ws::select(
  up, keys=df$uniprot, 
  columns=columns(up),
  keytype="UniProtKB") %>% 
  as_tibble() %>% dplyr::filter(Reviewed=="reviewed") %>% glimpse()
up_anno_selected <- up_anno %>% 
  dplyr::select(#query=From, 
                uniprot=Entry,
                # symbol='Gene.Names..primary.', alias='Gene.Names..synonym.',
                # protein_length=Length, 
                protein_uniprot=Sequence) %>%
  # mutate(across(entrez_id:ensembl_transcript_id, \(x) str_replace_all(x,";",""))) %>% 
  # mutate(ensembl_transcript_id=str_replace_all(ensembl_transcript_id,"[.]\\d?","")) %>% 
  glimpse()


df_uniprot <- left_join(df, up_anno_selected) %>% 
  glimpse()

# write_csv(df_uniprot, 'nmyc_related.csv')
















aa_up_dbd <- aa_up %>% 
  dplyr::select(uniprotswissprot=From, uniprot_aa_sequence=Sequence,
                region=Region, dbd=DNA.binding)
# aa_up_dbd$region[1:5]
aa_up_dbd <- aa_up_dbd %>% 
  mutate(DBD_start=str_extract(dbd, "DNA_BIND (\\d+)..(\\d+); /", group=1),
         DBD_stop=str_extract(dbd, "DNA_BIND (\\d+)..(\\d+); /", group=2)) %>% 
  print()



# ## fetch mobidb
# library(protti)
# mobi <- protti::fetch_mobidb(uniprot_ids=res3$uniprotswissprot)

## get activation domain and repression domain
# refer to: https://www.nature.com/articles/s41586-023-05906-y#Sec32
TF_AD_RD <- list(AD="Activation Domains", 
                 RD="Repression Domains") %>% 
  map(\(x) read_xlsx("41586_2023_5906_MOESM6_ESM.xlsx", x)) %>% 
  list_rbind(names_to="type") %>% 
  reframe(start=min(Start), stop=max(End),
          .by=c(type, `UniProt ID`)) %>% 
  pivot_wider(names_from=c(type), values_from=c(start, stop), 
              names_glue="{type}_{.value}", names_vary="slowest") %>% 
  print()

## get activation domain from paddle
# ref to: https://paddle.stanford.edu/
TF_AD <- list(high="human_TF_ADs_high-strength.tsv", 
              medium="human_TF_ADs_medium-strength.tsv") %>% 
  map(read_tsv) %>% 
  list_rbind(names_to="confidence") %>% 
  dplyr::select(!2) %>% 
  mutate(confidence=factor(confidence, c("high", "medium"))) %>% 
  slice_head(n=1, by=c(pID, prot_name)) %>% 
  separate(prot_name, c("symbol", "species"), sep="_") %>% 
  rename_with(\(x) str_c("AD2_", x), c(start, stop)) %>% 
  dplyr::select(pID:AD2_stop) %>% 
  print()




# merge with AD
res4 <- res3 %>% 
  left_join(TF_AD_RD, by=join_by(uniprotswissprot==`UniProt ID`)) %>% 
  left_join(TF_AD, by=join_by(uniprotswissprot==pID)) %>% 
  transmute(symbol=hgnc_symbol, 
            uniprot=uniprotswissprot,
            uniprot_aa=uniprot_aa_sequence, aa_length=uniprot_aa_length,
            cds=cdna, cdna_aa=cdna_aa_sequence, 
            if_uniprot_cdna_aa_same=if_aa_same,
            DBD_aa_start=DBD_start, DBD_aa_stop=DBD_stop,
            AD_aa_start=AD_start, AD_aa_stop=AD_stop,
            RD_aa_start=RD_start, RD_aa_stop=RD_stop,
            AD2_aa_start=AD2_start, AD2_aa_stop=AD2_stop,
            # AD_mean_Z=mean_Z, AD_confidence=confidence,
            uniprot_anno_region=region, uniprot_anno_dbd=dbd) %>% 
  glimpse()


# add tumor type
res5 <- tf %>% 
  reframe(
    tumor_type=str_c(unique(tumor_type), collapse="; "),
    tumor_subtype=str_c(unique(tumor_subtype), collapse="; "),
    .by=tf) %>% 
  dplyr::rename(symbol=tf) %>% 
  right_join(res4)


## save data
save.image("TF.RData")
write_xlsx(res5, file.path(dataset_dir, "TF.xlsx"))


res4 %>% 
  ggplot(aes(aa_length)) +
  geom_histogram(bins=100, color="black", fill="steelblue") +
  # geom_density(color="red3") +
  theme_bw()

count(res4, !is.na(DBD_aa_start)) # 56
count(res4, !is.na(AD_aa_start))  # 38
count(res4, !is.na(RD_aa_start)) # 102
count(res4, !is.na(AD2_aa_start)) # 48

# # write fasta
# fasta <- data[, c("fasta_id","protein_sequence_dRTK")] %>% 
#   deframe() %>% AAStringSet() %>% print()
# writeXStringSet(fasta, file.path(dataset_dir, "family_dRTK.fasta"))


