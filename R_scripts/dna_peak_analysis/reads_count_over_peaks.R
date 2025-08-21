
# set path
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |>
  setwd()

# load packages
library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(furrr)


if(F) {
  c("HUVEC_5mM","IMR_0.5mM","MEF_5mM")
  setwd("HUVEC_5mM")
  setwd("..")
  getwd()
}


#### define input data -------------------------------
# bam and peak files were co-coordinated in order
peak_files <- list.files(".", pattern=".narrowPeak$", full.names=F, recursive=F)
print(peak_files)

bam_files <- list.files(".", pattern=".bam$", full.names=F, recursive=F)
print(bam_files)

sample_names <- map_chr(peak_files, basename)



#### read bam and peak once -------------------------------
# collect result [list of vector]
stats <- list()

for(idx in seq_along(bam_files)) {
  # idx = 1
  
  # read the BED file
  # MACS2 bed is not standard bed format
  # use narrowPeak format instead if using rtracklayer import
  # bed <- rtracklayer::import(peak_files[idx])
  bed <- ChIPseeker::readPeakFile(peak_files[idx])
  
  # read the BAM file
  bam <- GenomicAlignments::readGAlignments(bam_files[idx])
  
  # count reads within genomic ranges
  reads_within_ranges <- GenomicRanges::countOverlaps(bam, bed, ignore.strand=T)
  # get sum of within peak reads
  reads_within_ranges <- sum(reads_within_ranges)
  
  # calculate total reads in the BAM file
  total_reads <- length(bam)
  
  # calculate reads outside of ranges
  reads_outside_ranges <- total_reads - reads_within_ranges
  
  # return result
  stats <- c(stats, list(c(within=reads_within_ranges, outside=reads_outside_ranges)))
  
  # print
  cat("Process sample:", sample_names[idx])
  cat("Reads within genomic ranges:", reads_within_ranges, "\n")
  cat("Reads outside of genomic ranges:", reads_outside_ranges, "\n")
  cat("within / outside:", round(100*(reads_within_ranges/reads_outside_ranges), 2), "%\n\n")
  
}




#### boostrap sampling -------------------------------

# use sampling to generate distribution

# single core
count_reads_with_sampling <- function(bam_file, 
                                      range_file,
                                      num_iterations = 100,
                                      sample_ratio = 0.5) {
  # read the bed file
  range_gr <- ChIPseeker::readPeakFile(range_file)
  
  # read the bam file
  bam <- GenomicAlignments::readGAlignments(bam_file)

  # calculate sample size
  sample_size <- round(sample_ratio*length(bam), 0)
  
  
  # using progress bar
  pb <- progress::progress_bar$new(total=num_iterations, 
                                   format="iteraction: [:bar] :percent in :elapsed",
                                   width=60)
  
  # initialize vector to store results for each iteration
  inside_counts <- numeric(num_iterations)
  outside_counts <- numeric(num_iterations)
  
  for(iteration in seq_len(num_iterations)) {
    # update progress bar
    pb$tick()
    
    # Generate random indices for sampling with replacement
    sample_indices <- sample(seq_along(bam), size = sample_size, replace = TRUE)
    # Subset the BAM file with sampled indices
    sampled_bam <- bam[sample_indices]
    
    # count reads within genomic ranges
    reads_within_ranges_sampled <- countOverlaps(sampled_bam, range_gr)
    
    # store the result for this iteration
    inside_counts[iteration] <- sum(reads_within_ranges_sampled)
    outside_counts[iteration] <- length(sampled_bam) - sum(reads_within_ranges_sampled)
  }
  
  result <- tibble(inside_counts=inside_counts, outside_counts=outside_counts)
  
  return(result)
}


# multiple cores
count_reads_with_sampling_mc <- function(bam_file, 
                                         range_file,
                                         num_iterations = 100,
                                         sample_ratio = 0.5) {
  # read the bed file
  range_gr <- ChIPseeker::readPeakFile(range_file)
  
  # read the bam file
  bam <- GenomicAlignments::readGAlignments(bam_file)

  # calculate sample size
  sample_size <- round(sample_ratio*length(bam), 0)
  
  
  # define single core function
  func <- function(iteraction, bam, range_gr, size) {
    # Generate random indices for sampling with replacement
    sample_indices <- sample(seq_along(bam), size = size, replace = TRUE)
    # Subset the BAM file with sampled indices
    sampled_bam <- bam[sample_indices]
    
    # count reads within genomic ranges
    reads_within_ranges_sampled <- countOverlaps(sampled_bam, range_gr)

    # return res
    tibble(inside_counts = sum(reads_within_ranges_sampled), 
           outside_counts = length(sampled_bam) - sum(reads_within_ranges_sampled))
  }
  
  
  # using multiple cores
  result_list <- future_map(seq_len(num_iterations), 
                            ~ func(.x, bam, range_gr, sample_size),
                            .options = furrr_options(seed=T),
                            .progress = T)
  
  # to tibble
  result <- bind_rows(result_list)
  
  return(result)
}



# run here

# for multiple cores
plan(multisession, workers=(parallelly::availableCores()-2))
# plan(sequential)
plan()

# collect all sample results [list of tibble]
stats_bs <- list()

for(idx in seq_along(bam_files)) {
  cat(">>>> index on:", idx, "\n")
  
  res <- count_reads_with_sampling_mc(bam_files[idx],
                                      peak_files[idx],
                                      num_iterations = 200,
                                      sample_ratio = 0.25)
  stats_bs[[sample_names[idx]]] <- res
}


# to tidy
stats_bs <- stats_bs %>%
  bind_rows(.id="sample") %>%
  mutate(inout_ratio = inside_counts / outside_counts)
# save
write_csv(stats_bs, str_glue("{Sys.Date()}_for_bs.csv"))




#### plot result -------------------------------

# by count ratio
stats_bs %>%
  # bind_rows(.id="sample") %>%
  # mutate(inout_ratio = inside_counts / outside_counts) %T>%
  # write_csv(str_glue("{Sys.Date()}_data_for_bs.csv")) %>%
  ggplot(aes(sample, inout_ratio, color=sample)) +
  # geom_jitter() +
  # geom_boxplot() +
  geom_violin(linewidth=0.7, width=0.5, show.legend=F) +
  labs(x="", y="Ratio of reads inside and outside of peaks") +
  # coord_cartesian(ylim=c(0.075, 0.105)) +
  theme_bw(12)
# save plot
cowplot::ggsave2(str_glue("{Sys.Date()}_for_bs.pdf"), 
                 width=5, height=5)
 


# by count ratio normalized by genomic length
# calculate peak length first
inside_peak_length <- peak_files %>%
  map(~ ChIPseeker::readPeakFile(.x) %>%
        GenomicRanges::reduce() %>%
        IRanges::width() %>%
        sum() %>%
        tibble(inside_peak_length=.)) %>%
  set_names(sample_names)
# hg38_legnth =  3096346000
peak_length <- bind_rows(inside_peak_length, .id="sample") %>%
  mutate(outside_peak_length = 3096346000 - inside_peak_length)
  
# range_gr <- ChIPseeker::readPeakFile(range_file)


# calculate normalized coverage
stats_bs %>%
#   bind_rows(.id="sample") %>%
#   mutate(inout_ratio = inside_counts / outside_counts) %>%
  # stats_bs %>%
  left_join(peak_length) %>%
    mutate(inside_counts_norm = inside_counts / inside_peak_length,
           outside_counts_norm = outside_counts / outside_peak_length,
           inout_norm_ratio = inside_counts_norm / outside_counts_norm) %T>%
  write_csv(str_glue("{Sys.Date()}_for_bs_norm.csv")) %>%
  ggplot(aes(sample, inout_norm_ratio, color=sample)) +
  # geom_jitter() +
  # geom_boxplot() +
  geom_violin(linewidth=0.7, width=0.5, show.legend=F) +
  labs(x="", y="Ratio of normalized inside and outside peak coverage") +
  # coord_cartesian(ylim=c(0.2, 0.4)) +
  theme_bw(12)
# save plot
cowplot::ggsave2(str_glue("{Sys.Date()}_for_bs_norm.pdf"), 
                 width=5, height=5)
