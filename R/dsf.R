
library(tidyverse)
library(rio)
options(rio.import.class='tbl')
library(patchwork)
Sys.setenv(JAVA_HOME="")
library(rcdk)

source('C:/Users/haohe/GitHub/myScript/R/helper.R')

# D:\\LGLab\\Project_p53\\2026-03-04_hits_dsf_10_cmpds_w_conc'
root_dir <- norm_path(r'(D:\LGLab\Project_p53\2026-03-24_dsf_more_protein)')
setwd(root_dir); print(getwd())

# set column target/protein -> protein
# set column cmpd/molecule/name -> ligand

# set control
ref_protein <- 'NoProtein'
ref_ligand <- 'DMSO'


# # additional metadata like smiles, chemical properties, etc. can be added here
# metadata <- 
#   import('../2026-02-01_hits_dsf_2nd/2026-02-05_merged.xlsx') %>%
#   # import('../2026-01-25_hits_dsf/2026-01-26_merged.xlsx') %>% 
#   # rename(ligand=name) %>% 
#   select(c(ligand, smiles)) %>%
#   print()

info <- #c('source_cmpd_p1.xlsx', 'source_cmpd_p2.xlsx') %>% 
  c('plate_info.xlsx') %>% 
  # set_names(c('p1','p2')) %>% 
  map(read_metadata) %>% 
  list_rbind(names_to='plate') %>% 
  janitor::remove_constant() %>% 
  # rename(ligand=cmpd) %>% 
  # left_join(previous) %>% 
  print()

# info <- info %>% 
#   mutate(smiles=ifelse(ligand=='PC14586', 'N(C1=C2C(N(CC(F)(F)F)C(C#CCNC3=C(OC)C=C(C(NC)=O)C=C3)=C2)=CC=C1)[C@H]4[C@@H](F)CN(C)CC4', smiles)) %>% 
#   mutate(smiles=ifelse(ligand=='PhiKan083', 'CCN1C3=C(C=CC=C3)C2=C1C=CC(CNC)=C2.Cl', smiles))


## analyzed data ----------------------------

# merge with metadata
raw <- #'DSF_standard_analysis_20260304_180954_AnalysisResults.txt' %>% 
  list.files('.', '_AnalysisResults.txt') %>% 
  import() %>% 
  janitor::clean_names() %>% 
  glimpse()

merged <- raw %>% 
  janitor::remove_empty('cols') %>%
  mutate(well=transform_well_style(well)) %>% #view()
  select(plate=experiment_file_name, 
         # any_of(c('target','ligand')),
         well, 
         tm=tm_d#, poor_fit
  ) %>% 
  mutate(plate=map_chr(plate, \(x) str_replace_all(x, '.*_Admin_', '') %>% 
                         str_replace_all('.eds', '')),
         # target=plate
  ) %>% 
  right_join(info, .) %>% 
  filter(!is.na(ligand), !is.na(target)) %>%
  relocate(plate, well, target, .before=1) %>% 
  relocate(tm, .after=last_col()) %>% 
  glimpse()


# calculate delta TM
ref <- merged %>% 
  filter(ligand %in% ref_ligand) %>%
  reframe(tm_ref=median(tm), .by=c(target)) %>%
  print()
df <- merged %>% 
  left_join(ref) %>% 
  mutate(d_tm=tm - tm_ref) %>% 
  select(!tm_ref) %>% 
  relocate(d_tm, .after=tm) %>% 
  { if(any(str_detect(colnames(.), 'conc'))) {
    mutate(., conc=as.numeric(conc))
  } else . } %>% 
  glimpse()
# filter(df, ligand==ref_ligand) %>% select(tm, target, well)

# merge to previous merged dataset
df_wide <- df %>% 
  # mutate(target=str_c(target, '_2nd')) %>% 
  pivot_wider(id_cols=any_of(setdiff(colnames(info), 'well')), 
              names_from=target, values_from=d_tm,
              names_prefix='dTM_',
              values_fn=median) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>% 
  glimpse()
export(df_wide, str_glue('{Sys.Date()}_merged.xlsx'))



# plot
p <- df %>% 
  # filter(!target %in% c('NoProtein')) %>%
  mutate(d_tm=ifelse(d_tm > 10, 10, ifelse(d_tm < -10, -10, d_tm))) %>% 
  ggplot(aes(ligand, tm, color=target)) +
  geom_point(show.legend=F) +
  facet_wrap(vars(target), ncol=1, scales='fixed', axes='all') +
  labs(x='', y='TM of Derivative') +
  theme_bw() +
  theme(#axis.text.x=element_text(angle=90, vjust=1, size=4),
        axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.ticks.x=element_blank())
ggsave(str_glue('{Sys.Date()}_tm.pdf'), p, width=7, height=9)

p <- df %>% 
  # filter(!target %in% c('NoProtein')) %>%
  # filter(!(ligand %in% c('Nutlin'))) %>% 
  mutate(d_tm=ifelse(d_tm > 10, 10, ifelse(d_tm < -10, -10, d_tm))) %>% 
  ggplot(aes(ligand, d_tm, color=target)) +
  geom_point(show.legend=T) +
  facet_wrap(vars(target), ncol=1, scales='fixed', axes='all') +
  labs(x='', y='Delta TM of Derivative') +
  theme_bw() +
  theme(#axis.text.x=element_text(angle=90, vjust=1, size=4),
        axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.ticks.x=element_blank())
ggsave(str_glue('{Sys.Date()}_d_tm.pdf'), p, width=7, height=5)





## raw data ----------------------------

# calculate derivatives
caculate_deri <- function(data, eps=1e-07, unconditional=F) {
  mod = mgcv::gam(fluorescence ~ s(temperature), data=data)
  res <- gratia::derivatives(mod, n=100, order=1, eps=eps, unconditional=unconditional)
  tibble(temperature=res$temperature, 
         derivative=res$.derivative, 
         derivative_se=res$.se) }

mc1 <- list.files('.', '_RawData_.*.txt') %>% 
  set_names(., nm=map_chr(., \(x) str_replace_all(x, '.*_Admin_', '') %>% 
                            str_replace_all('.eds.txt', ''))) %>% 
  # c('DSF_standard_analysis_20260125_173825_RawData_standard_DSF_SYPRO_Orange_384_20260125_165218_Admin.eds.txt',
  #        'DSF_standard_analysis_20260125_173825_RawData_standard_DSF_SYPRO_Orange_384_20260125_172455_Admin.eds.txt') %>% 
  # set_names(nm=c(
  #   'standard_DSF_SYPRO_Orange_384_20260125_165218_Admin.eds',
  #   'standard_DSF_SYPRO_Orange_384_20260125_172455_Admin.eds'
  # )) %>% 
  map(\(x) import(x)) %>% 
  list_rbind(names_to='plate') %>% 
  janitor::clean_names() %>% 
  select(!c(well)) %>% 
  rename(well=well_position) %>% 
  mutate(well=transform_well_style(well)) %>%
  glimpse()
count(mc1, plate, well)

# df$plate %>% unique()
mc_tidy <- mc1 %>% 
  # rename(target=plate) %>% 
  left_join(info) %>% 
  mutate(temperature=as.numeric(temperature),
         fluorescence=as.numeric(fluorescence)) %>% 
  reframe(fluorescence=mean(fluorescence), .by=c(plate, well, target, ligand, temperature)) %>% 
  filter(!is.na(ligand), !is.na(target)) %>% 
  glimpse()

rois <- mc_tidy$ligand %>% unique() %>% print()

# plot setting
# if color_by_column not existed in df, or NA, NULL, then set color to 'darkred'
# no matter how color_by set, ref_ligand will always be grey color
# if only have more than one plate, then add plate into facetting variable
color_by_column <- 'conc' # plate, target, ligand, NA

is_valid_color_by <- function(df, col) {
  !is.null(col) && !is.na(col) && nzchar(col) && col %in% names(df)
}

plist <- list()
for (roi in rois) {
  print(roi)

  df_sub <- mc_tidy %>%
    filter(ligand %in% c(ref_ligand, roi))

  multi_plate <- n_distinct(df_sub$plate) > 1
  facet_vars <- if (multi_plate) vars(target, plate) else vars(target)

  color_by_used <- if (is_valid_color_by(df_sub, color_by_column) && color_by_column != 'ligand') {
    color_by_column
  } else {
    NA_character_
  }

  is_color_numeric <- FALSE
  if (!is.na(color_by_used)) {
    is_color_numeric <- is.numeric(df_sub[[color_by_used]])
  }

  # raw curve
  p1 <- ggplot() +
    geom_path(
      data = filter(df_sub, ligand == ref_ligand),
      aes(temperature, fluorescence, group = well),
      color = 'grey50',
      show.legend = TRUE
    )

  if (!is.na(color_by_used)) {
    p1 <- p1 +
      geom_path(
        data = filter(df_sub, ligand != ref_ligand),
        aes(temperature, fluorescence, color = .data[[color_by_used]], group = well),
        show.legend = TRUE
      ) +
      (if (is_color_numeric) scale_color_viridis_c(option = 'C', na.value = 'darkred')
       else scale_color_viridis_d(option = 'C', na.value = 'darkred'))
  } else {
    p1 <- p1 +
      geom_path(
        data = filter(df_sub, ligand != ref_ligand),
        aes(temperature, fluorescence, group = well),
        color = 'darkred',
        show.legend = TRUE
      )
  }

  p1 <- p1 +
    guides(color = guide_legend('')) +
    facet_wrap(facet_vars, ncol = 6) +
    labs(title = roi) +
    theme_bw() +
    theme(
      legend.position = 'top',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  p1

  # tm values
  tm <- df %>%
    select(any_of(c('target', 'ligand', 'conc', 'tm'))) %>%
    reframe(tm = mean(tm), .by = any_of(c('target', 'ligand', 'conc'))) %>%
    filter(ligand %in% c(ref_ligand, roi)) %>%
    mutate(tm = ifelse(target == ref_protein, NA, tm))

  # derivative curve
  deri <- df_sub %>%
    nest(.by = c(plate, target, ligand, well)) %>%
    mutate(pred = map(data, caculate_deri)) %>%
    unnest(pred)

  p2 <- ggplot() +
    geom_line(
      data = filter(deri, ligand == ref_ligand),
      aes(temperature, derivative, group = well),
      color = 'grey50',
      show.legend = TRUE)
  if (!is.na(color_by_used)) {
    p2 <- p2 +
      geom_line(
        data = filter(deri, ligand != ref_ligand), show.legend = TRUE,
        aes(temperature, derivative, color = .data[[color_by_used]], group = well)) +
      (if (is_color_numeric) scale_color_viridis_c(option = 'C', na.value = 'darkred')
       else scale_color_viridis_d(option = 'C', na.value = 'darkred'))
  } else {
    p2 <- p2 +
      geom_line(
        data = filter(deri, ligand != ref_ligand), show.legend = TRUE,
        aes(temperature, derivative, group = well), color = 'darkred')
  }

  p2 <- p2 +
    facet_wrap(facet_vars, ncol = 6) +
    coord_cartesian(xlim = c(25, 55)) +
    labs(title = roi) +
    theme_bw() +
    theme(
      legend.position = 'none',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  p2

  # 2D chemical structure
  if (roi %in% info$ligand & 'smiles' %in% colnames(info)) {
    try_res <- try({
      mol <- info %>%
        filter(ligand == roi) %>%
        pull(smiles) %>%
        unique() %>%
        rcdk::parse.smiles() %>%
        .[[1]]
      grob <- grid::rasterGrob(view.image.2d(mol))
      p3 <- ggplot() +
        annotation_custom(grob) +
        labs(title = roi) +
        theme_void()
    })
    if (class(try_res)[1] == "try-error") {
      plist[[roi]] <- p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2))
    } else {
      plist[[roi]] <- p1 + p2 + p3 + plot_layout(widths = c(1, 1, 0.2))
    }
  } else {
    plist[[roi]] <- p1 + p2 + plot_spacer() + plot_layout(widths = c(1, 1, 0.2))
  }
}

pdf(str_glue('{Sys.Date()}_ligand.pdf'), width = 15, height = 4)
for (p in plist) print(p)
while (dev.cur() != 1) dev.off()



  