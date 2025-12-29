
# set path
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |>
  setwd()

# load packages
library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(GenomeInfoDb)
library(furrr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
options(future.globals.maxSize = 2 * 1024 * 1024^2) # 2 GB
BiocManager::install("ChIPseeker")

setwd("/mnt/f/workspace/IMR/result/")


## define input -------------------------------
peak_files <- list.files("04_peaks", pattern="broadPeak$", full.names=T, recursive=F)
print(peak_files)
df <- tibble(
  peak_file=peak_files,
  sample=map_chr(peak_files, basename) %>% 
    str_replace_all("_peaks.*Peak$",""),
  bam_file=str_c("03_alignment/", sample, ".filtered.bam"),
) %>% print()


## consensus peak -------------------------------
consensus <- rtracklayer::import("04_peaks/consensus.bed")
# consensus_annno <- ChIPseeker::annotatePeak(
#   consensus, 
#   tssRegion = c(-3000, 3000),
#   TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#   level = "gene")
# consensus <- consensus_annno@anno
consensus <- keepStandardChromosomes(
  consensus, pruning.mode="coarse")


## gene/tss -------------------------------
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) %>% 
  keepStandardChromosomes(pruning.mode="coarse") %>% 
  sort() %>% print()

# used in following functions
genes_reduced <- genes %>% 
  resize(width=width(genes) + 2000, fix="center") %>% 
  reduce() #%>% print()
tss_reduced <- promoters(genes, 1000, 1000) %>% 
  reduce() #%>% print()





compute_region <- function(gr, g=genes_reduced, t=tss_reduced) {
  in_gene <- subsetByOverlaps(gr, g, ignore.strand=FALSE) %>% 
    reduce() #%>% print()
  flank_gene <- c(
    flank(in_gene, width=10000, start=T, both=F),
    flank(in_gene, width=10000, start=F, both=F)
  ) %>% reduce()# %>% print()

  in_tss <- subsetByOverlaps(gr, t, ignore.strand=FALSE) %>% 
    reduce() #%>% print()
  flank_tss <- c(
    flank(in_tss, width=3000, start=T, both=F),
    flank(in_tss, width=3000, start=F, both=F)
  ) %>% reduce() #%>% print()
  
  return(list(in_gene, flank_gene, in_tss, flank_tss) %>% 
           set_names(nm=c("in_gene", "flank_gene", "in_tss", "flank_tss")))
}

do_sampling <- function(peak, peak_name="peak", sample_size=1000) {
  peak_sampled <- sample(consensus, size=sample_size, replace=F)
  peak_regions <- compute_region(peak_sampled)
  names(peak_regions) <- str_c(peak_name,"_", names(peak_regions))
  peak_regions[[peak_name]] <- peak_sampled
  return(peak_regions)
}
# do_sampling(consensus, "consensus")

get_stats <- function(bam_path, 
                      peak_path, 
                      g=genes_reduced,
                      t=tss_reduced,
                      consens=consensus,
                      sample_rep=NULL, 
                      sample_size=1000) {
  
  # bam_path <- df$bam_file[1]
  bam_gr <- GenomicAlignments::readGAlignmentPairs(bam_path) %>% 
    granges()
  
  if(sample_size <= 1) sample_size <- floor(length(bam_gr)*sample_size)
  if(sample_size > length(bam_gr)) sample_size <- length(sample_size)
  
  # peak_path <- df$peak_file[1]
  peak <- rtracklayer::import(peak_path) %>% reduce()
  
  # set rep
  if(is.null(sample_rep)) {
    sample_rep <- 1
    sample_size <- length(bam_gr)
  }
  
  # run here
  res <- list()
  for(i in seq_len(sample_rep)) {
    gene_sampled=sample(g, sample_size, replace=F)
    tss_sampled=sample(t, sample_size, replace=F)
    regions <- c(gene=list(gene_sampled), 
                 tss=list(tss_sampled),
                 gene_flank=c(
                   flank(gene_sampled, width=10000, start=T, both=F),
                   flank(gene_sampled, width=10000, start=F, both=F)
                 ) %>% reduce() %>% list(),
                 tss_flank=c(
                   flank(tss_sampled, width=3000, start=T, both=F),
                   flank(tss_sampled, width=3000, start=F, both=F)
                 ) %>% reduce() %>% list() )
    if(!is.null(peak_path)) regions <- c(regions, do_sampling(peak, "peak", sample_size))
    if(!is.null(consens)) regions <- c(regions, do_sampling(consens, "consensus", sample_size))
    # print(regions)
    
    gr_counts <- map_int(
      regions, \(x) countOverlaps(bam_gr, x, ignore.strand=T) %>% sum()) %>% 
      enframe() %>% 
      pivot_wider()
    
    gr_counts_list <- list(gr_counts) %>% set_names(i)
    
    res <- c(res, gr_counts_list)
  }
  
  data <- list_rbind(res, names_to="rep") %>% 
    mutate(total=length(bam_gr), .before=1)
  
  return(data)
}



# df2 <- df[c(1,3),]
df2 <- df

stats <- map2(df2$bam_file, df2$peak_file, 
              \(x, y) get_stats(x, y, genes_reduced, tss_reduced, 
                                consensus, sample_rep=100, sample_size=1000)) %>% 
  set_names(nm=df2$sample)


stats_tidy <- stats %>% 
  list_rbind(names_to="sample") %>% #glimpse()
  # dplyr::rename(gene_flank=genes_flank) %>% 
  mutate(peak_of_total=peak/total,
         consensus_of_total=consensus/total,
         gene_of_total=gene/total,
         tss_of_total=tss/total,
         gene_in_of_flank=gene/gene_flank,
         tss_in_of_flank=tss/tss_flank,
         consensus_gene_in_of_flank=consensus_in_gene/consensus_flank_gene,
         consensus_tss_in_of_flank=consensus_in_tss/consensus_flank_tss,
         peak_in_of_flank=peak_in_gene/peak_flank_gene,
         peak_tss_in_of_flank=peak_in_tss/peak_flank_tss
  ) %>% 
  # view()
  glimpse()

stats_tidy %>% 
  reframe(across(everything(), mean), .by=sample) %>%
  view()
writexl::write_xlsx(stats_tidy, str_glue("stats_tidy.xlsx"))



df <- read_csv("D:\\LGLab\\Others\\to_huang\\frap\\frap.csv") %>% 
  janitor::clean_names()
df1 <- df %>% 
  pivot_longer(!axis_s) %>% 
  separate_wider_delim(name, delim="_roi_", names=c("channel","roi")) %>% 
  reframe(axis_s=axis_s, value=value, ratio=value/max(value), .by=roi) %>% 
  filter(axis_s <= 60)

df1 %>% 
  ggplot(aes(axis_s, ratio)) +
  # stat_summary(geom="ribbon", size=0.1, linewidth=0.25, fill="grey80") +
  stat_summary(geom="pointrange", size=0.1, linewidth=0.25, color="steelblue") +
  stat_summary(geom="line", linewidth=0.5, color="steelblue") +
  # geom_line() +
  labs(x="Time (s)", y="Recovery Ratio") +
  theme_bw()
ggsave("frap2.pdf", width=3, height=2)
writexl::write_xlsx(df1, "frap.xlsx")




# plot

data <- readxl::read_excel("2025-12-23_stats_tidy.xlsx") %>% 
  glimpse()













