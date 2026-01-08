
library(tidyverse)
library(RSQLite)
library(ggpubr)

# rstudioapi::getActiveDocumentContext()$path %>%
#   dirname() %>% setwd()
# getwd()
setwd("/mnt/f/workspace/IF_H3K27me3_analysis/")

# use_db <- "hist_100_60000_quantile_50_logspace"
# use_db <-  "hist_100_8000_quantile_50_linear"
# use_db <- "hist_100_60000_quantile_100_linear"
# use_db <- "hist_150_65535_mean_50_linear"
# use_db <- "hist_150_65535_mean_50_sqrt"
# use_db <- "hist_150_65535_mean_50_linear"

for(use_db in c(#"hist_100_60000_quantile_50_logspace",
                "hist_100_6000_mean_50_geomspace",
                "hist_150_8000_mean_50_geomspace",
                "hist_150_8000_mean_50_linear",
                "hist_150_65535_mean_50_geomspace",
                "hist_200_6000_mean_50_geomspace")) {
  
  db_files <- list.files(".", pattern = ".db", recursive = TRUE, full.names = TRUE) %>%
    str_subset(use_db) %>%
    print()

  read_db <- function(f, table_name = "cell") {
    conn <- DBI::dbConnect(RSQLite::SQLite(), f)
    df <- tbl(conn, table_name) %>% collect()
    DBI::dbDisconnect(conn)
    return(df)
  }
  
  
  raw_list <- map(db_files, read_db, table_name="cell")
  
  raw <- list_rbind(raw_list) %>% glimpse()
  # count(raw, directory)
  
  hist <- raw %>% 
    mutate(
      # Extract clean directory name
      directory=basename(dirname(directory)) %>% str_remove_all("__20.*"),
      # Compute bin edge mean
      bin_edge_mean=(bin_edge_low + bin_edge_high) / 2,
      # Cell type
      celltype=ifelse(str_detect(directory, fixed("MEF")), "MEF", "IMR90"),
      # Assign group — NO reference to 'group' inside!
      group=case_when(
        directory == "20251114 MEF KO H3K27me3 IF" ~
          case_when(
            well == "B6" ~ "DMSO",
            well == "C6" ~ "NMN_5mM",
            well == "D6" ~ "NMN_10mM",
            TRUE ~ NA_character_
          ),
        directory == "LH 20231107 IMR90 NMN treatment H3K27me3" ~
          case_when(
            well == "B1" ~ "DMSO",
            well %in% c("B2", "B3", "B4", "B5", "B6", "C1", "C2", "C3", "C4", "C5", "C6", 
                        "D1", "D2", "D3", "D4", "D5", "D6") ~ "NMN",  # or just TRUE if all others are NMN
            TRUE ~ NA_character_
          ),
        directory == "LH 20231107 IMR90 young H3K27me3" ~
          ifelse(well == "D4", "Untreat", NA_character_),
        directory == "20251201 MEF WT H3K27me3 IF NMN 226" ~
          case_when(
            well == "C2" ~ "DMSO",
            well == "C3" ~ "NMN_5mM",
            well == "C4" ~ "NMN_10mM",
            TRUE ~ NA_character_
          ),
        directory == "20251211 MEF WT NMN EED226 H3K27me3 IF" ~
          case_when(
            well == "D2" ~ "DMSO",
            well == "D3" ~ "NMN_3mM",
            well == "D4" ~ "NMN_10mM",
            TRUE ~ NA_character_
          ),
        directory == "20251214 IMR90 H3K27me3 IF NMN 226" ~
          case_when(
            well == "B2"  ~ "DMSO",
            well == "B3"  ~ "EED226_1uM",
            well == "B4"  ~ "EED226_3uM",
            well == "B5"  ~ "EED226_10uM",
            well == "C1"  ~ "NMN_1uM",
            well %in% c("D2", "C3") ~ "NMN_5uM",
            TRUE ~ NA_character_
          ),
        directory == "20251222 old plate IMR90 NMN" ~
          case_when(
            well == "B4" ~ "DMSO",
            well == "A5" ~ "NMN_0.3mM",
            well == "B6" ~ "NMN_0.6mM",
            TRUE ~ NA_character_
          ),
        directory == "20251222 old plate2 IMR90 NMN" ~
          case_when(
            well == "C5" ~ "DMSO",
            well == "C3" ~ "NMN_0.5mM",
            well == "C4" ~ "NMN_1mM",
            TRUE ~ NA_character_
          ),
        directory == "20251223 MEF H3K27me3 IF NMN" ~
          case_when(
            well == "C2" ~ "DMSO",
            well == "C3" ~ "NMN_0.5mM",
            TRUE ~ NA_character_
          ),
        TRUE ~ NA_character_
      )
    ) %>%
    select(-timepoint, -stack)
  
  # count(hist, directory)
  # # filter(hist, is.na(group)) %>% pull(directory)
  # # hist$celltype %>% unique()
  # # hist$bin_index %>% unique()
  # 
  # cell <- hist %>%
  #   distinct(directory, well, field, cell_id, total_pixels, cell_mean_intensity) %>%
  #   glimpse()
  # 
  # cell$cell_mean_intensity %>% expm1() %>% 
  #   quantile(c(0.01, 0.02,0.05,0.25,0.5,0.75,0.95,0.98,0.99))
  # cell$total_pixels %>% 
  #   quantile(c(0.01, 0.02,0.05,0.25,0.5,0.75,0.95,0.98,0.99))
  # cell %>% 
  #   # mutate(cell_mean_intensity=expm1(cell_mean_intensity)) %>% 
  #   ggplot(aes(x=total_pixels, y=cell_mean_intensity, color=well)) +
  #   # geom_histogram(bins=200) +
  #   geom_point() +
  #   # geom_line() +
  #   # geom_density() +
  #   # geom_smooth() +
  #   # facet_wrap(vars(well), scales='fixed') +
  #   coord_cartesian(xlim=c(0, 15000), ylim=c(0, 8000)) +
  #   theme_bw()
  
  
  glimpse(hist)
  hist_filter <- hist %>% 
    filter(!is.na(group)) %>% 
    # filter(directory==unique(hist_selected$directory)) %>%
    # filter(between(cell_mean_intensity, 300, 6000)) %>%
    # filter(between(expm1(cell_mean_intensity), 300, 6000)) %>%
    # filter(directory=="20251114 MEF KO H3K27me3 IF") %>%
    # filter(well=="B6") %>% 
    # filter(field %in% c("1","2")) %>% 
    dplyr::filter(between(total_pixels, 4500, 9000))
  
  
  
  
  hist_cdf <- hist_filter %>% 
    reframe(bin_mean=mean(bin_mean),
            bin_edge_mean=mean(bin_edge_mean),
            bin_count=sum(bin_count),
            .by=c(directory, group, celltype, bin_index)) %>% 
    group_by(directory, group, celltype) %>%
    # Ensure bins are in order
    arrange(bin_index, .by_group=T) %>%  
    # weighted CDF for binned data
    mutate(cdf=cumsum(bin_count) / sum(bin_count) ) %>%
    ungroup() %>% 
    arrange(celltype) %>% 
    split(., .[['directory']])
  p2 <- hist_cdf %>% 
    imap(\(x, y) x %>% 
           ggplot(aes(bin_index, cdf, color=group)) +
           # geom_step(direction="hv") +
           geom_line() +
           labs(x="Intensity Bins", y="Cumulative Fraction of Pixels", title=y) +
           theme_bw(8) +
           theme(panel.grid=element_blank())
    ) %>% 
    patchwork::wrap_plots(ncol=1)
  # print(p2)
  ggsave(str_glue("{Sys.Date()}_cdf_{use_db}.pdf"), 
         p2, width=4, heigh=6)
  
  # statistics
  stats <- hist_filter %>% 
    split(., .[['directory']]) %>% 
    map(function(x) {
      if("DMSO" %in% x[["group"]] & length(unique(x[["group"]])) > 1) {
        dmso <- filter(x, group == "DMSO") %>% pull(bin_count)
        setdiff(unique(x[["group"]]), "DMSO") %>% 
          set_names(., nm=.) %>% 
          map(\(y) filter(x, group == y) %>% pull(bin_count) %>% 
                ks.test(., dmso) %>% 
                broom::tidy()) %>% 
          list_rbind(names_to="group")
      } else return(NULL)
    }) %>% 
    list_rbind(names_to="directory") %>% 
    dplyr::rename(ks_distance=statistic)
  writexl::write_xlsx(stats, str_glue("{Sys.Date()}_cdf_{use_db}.xlsx"))
  
  
}



# write_rds(hist, str_glue("{Sys.Date()}_histogram_200_40000.rds"))




# ggsave(str_glue("{Sys.Date()}_histogram_selected.pdf"), 
#        patchwork::wrap_plots(p, ncol=1),
#        width=4, heigh=3.5)
# write_rds(hist, str_glue("{Sys.Date()}_histogram_selected.rds"))


# n_count <- cell %>%
#   count(class, group) %>% filter(class=="C1") %>% 
#   summarise(text=paste0(paste(group, "=", n, collapse=" | "), collapse=""), .groups="drop") %>%
#   pull(text) %>% print()
# 
# cell %>%
#   ggplot(aes(x=class, y=count_prop, fill=group)) +
#   geom_boxplot(outlier.shape=NA, width=0.6, staplewidth =0.75,
#                position=position_dodge(width=0.75)) +
#   stat_compare_means(aes(group=group), method="wilcox.test", 
#                      label="p.signif", size=3, label.y.npc=0.92,hide.ns=TRUE) +
#   scale_fill_manual(values=c("DMSO"="#C7D3D4", "NMN"="steelblue")) +
#   theme_bw() +
#   theme(
#     legend.position="top",
#     panel.grid.minor.x=element_blank(),
#     panel.grid.major.x=element_blank(),
#     panel.grid.minor.y=element_blank()) +
#   labs(x="Chromatin Density Classes", y="Pixel Proportion of Total", caption=n_count) +
#   ylim(NA, max(cell$count_prop, na.rm=TRUE) * 0.5)
# 
# ggsave("MEF_CDCs_boxplot.pdf", width=7, heigh=4)
# writexl::write_xlsx(list(cell=cell, stat=stats_table), "MEF_CDCs_data.xlsx")


