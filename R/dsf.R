
library(tidyverse)
library(rio)
options(rio.import.class='tbl')
library(patchwork)
Sys.setenv(JAVA_HOME="")
library(rcdk)

source('C:\\Users\\haohe\\GitHub\\myScript\\R\\read_plate_info.R')

setwd('D:\\LGLab\\Project_p53\\2026-02-01_hits_dsf_2nd')

# set column target/protein -> protein
# set column cmpd/molecule/name -> ligand

# set control
ref_protein <- 'NoProtein'
ref_ligand <- 'DMSO'


## functions----------------------------
# transform A01, B13 to A1, B13 style
transform_well_style <- function(wells) sub("^([A-Za-z]+)0+", "\\1", wells)

# transform_well_style <- function(well) {
#   tmp <- str_extract(well, "(?<row>[A-Za-z])(?<column>[0-9]+)", group=1:2)
#   # print(tmp)
#   str_c(tmp[1], as.integer(tmp[2]))
# }

# calculate derivatives
caculate_deri <- function(data, eps=1e-07, unconditional=F) {
  mod = mgcv::gam(fluorescence ~ s(temperature), data=data)
  res <- gratia::derivatives(mod, n=100, order=1, eps=eps, unconditional=unconditional)
  tibble(temperature=res$temperature, 
         derivative=res$.derivative, 
         derivative_se=res$.se) }


# metadata
previous <- 
  import('../Stage1/ASMS_confirm/merged.xlsx') %>%
  # import('../2026-01-25_hits_dsf/2026-01-26_merged.xlsx') %>% 
  rename(ligand=name) %>% 
  # select(!well) %>% 
  print()

info <- read_plate_info('source_cmpd.xlsx') %>% #print()
  rename(ligand=cmpd) %>% 
  left_join(previous) %>% 
  print()


## analyzed result ----------------------------

# merge with metadata
raw <- 'DSF_standard_analysis_20260204_180724_AnalysisResults.txt' %>% 
  import() %>% 
  janitor::clean_names() %>% 
  glimpse()

merged <- raw %>% 
  janitor::remove_empty('cols') %>%
  mutate(well=transform_well_style(well)) %>% 
  select(plate=experiment_file_name, well, tm=tm_d, #poor_fit, 
         any_of(c('protein','ligand'))) %>% 
  right_join(info, .) %>% 
  filter(!is.na(ligand), !is.na(protein)) %>% 
  # select(!well) %>% 
  glimpse()


# calculate delta TM
ref <- merged %>% 
  filter(ligand %in% ref_ligand) %>%
  reframe(tm_ref=median(tm), .by=c(protein)) %>%
  print()
df <- merged %>% 
  left_join(ref) %>% 
  mutate(d_tm=tm - tm_ref) %>% 
  select(!tm_ref) %>% 
  relocate(d_tm, .after=tm) %>% 
  glimpse()
# filter(df, ligand==ref_ligand) %>% select(tm, protein, well)

# merge to previous merged dataset
df_wide <- df %>% 
  # mutate(protein=str_c(protein, '_2nd')) %>% 
  pivot_wider(id_cols=any_of(setdiff(colnames(info), 'well')), 
              names_from=protein, values_from=d_tm,
              names_prefix='dTM_',
              values_fn=median) %>%
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>% 
  glimpse()
export(df_wide, str_glue('{Sys.Date()}_merged.xlsx'))



# plot
p <- df %>% 
  # filter(!protein %in% c('NoProtein')) %>%
  mutate(d_tm=ifelse(d_tm > 10, 10, ifelse(d_tm < -10, -10, d_tm))) %>% 
  ggplot(aes(ligand, d_tm, color=protein)) +
  geom_point(show.legend=F) +
  facet_wrap(vars(protein), ncol=1, scales='fixed', axes='all') +
  labs(x='', y='TM of Derivative') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=1, size=4))
ggsave(str_glue('{Sys.Date()}_d_tm.pdf'), p, width=15, height=12)





## raw data analysis ----------------------------
mc1 <- list.files('.', '_RawData_') %>% 
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
  rename(protein=plate) %>% 
  left_join(df) %>% 
  mutate(temperature=as.numeric(temperature),
         fluorescence=as.numeric(fluorescence)) %>% 
  reframe(fluorescence=mean(fluorescence), .by=c(plate, protein, ligand, temperature)) %>% 
  glimpse()

rois <- df %>% 
  filter(d_tm > 1, protein != ref_protein) %>% 
  slice_max(d_tm, n=20, by=protein) %>% 
  arrange(desc(d_tm)) %>% 
  pull(ligand) %>% 
  unique() %>% 
  print()

plist <- list()
for(roi in rois) {
  # roi <- 'WXHTS0557836'
  print(roi)
  
  # raw curve
  p1 <- mc_tidy %>%
    filter(ligand %in% c(ref_ligand, roi)) %>%
    ggplot(aes(temperature, fluorescence, color=ligand, group=ligand)) +
    geom_path(show.legend=F) +
    scale_color_manual(values=c('grey50', 'darkred') %>% set_names(c(ref_ligand, roi))) +
    facet_wrap(vars(protein), ncol=6) +
    labs(title=roi) +
    theme_bw() +
    theme(legend.position='top',
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  
  # derivatives
  tm <- df %>% 
    select(protein, ligand, tm) %>%
    reframe(tm=mean(tm), .by=c(protein, ligand)) %>% 
    filter(ligand %in% c(ref_ligand, roi)) %>% 
    mutate(tm=ifelse(protein==ref_protein, NA, tm))
  p2 <- mc_tidy %>% 
    filter(ligand %in% c(ref_ligand, roi)) %>% 
    nest(.by=c(protein, ligand)) %>%
    mutate(pred=map(data, caculate_deri)) %>% 
    unnest(pred) %>% 
    ggplot(aes(temperature, derivative, color=ligand, group=ligand)) +
    # geom_ribbon(aes(ymin=derivative-derivative_se, ymax=derivative+derivative_se), 
    #             alpha=0.25, color='grey50', show.legend=F) +
    geom_line(show.legend=F) +
    geom_vline(aes(xintercept=tm, color=ligand), linetype=2, linewidth=0.5, data=tm, show.legend=F) +
    scale_color_manual(values=c('grey50', 'darkred') %>% set_names(c(ref_ligand, roi))) +
    facet_wrap(vars(protein), ncol=6) +
    coord_cartesian(xlim=c(20, 60)) +
    labs(title=roi) +
    theme_bw() +
    theme(legend.position='top',
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  
  # 2D chemical structure
  if(roi %in% info$ligand) {
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
        labs(title=roi) +
        theme_void()
    })
    
    if(class(try_res)[1] == "try-error") {
      plist[[roi]] <- p1 + p2 + plot_spacer() + plot_layout(widths=c(1,1,0.2))
    } else {
      plist[[roi]] <- p1 + p2 + p3 + plot_layout(widths=c(1,1,0.2))
    }
    
  } else {
    plist[[roi]] <- p1 + p2 + plot_spacer() + plot_layout(widths=c(1,1,0.2))
  }

}
# ggsave('test.pdf', p1 + p2 + p3 + plot_layout(widths=c(1,1,0.2)), width=20, height=2.5)
pdf(str_glue('{Sys.Date()}_ligand.pdf'), width=24, height=2.5)
for(p in plist) print(p)
while(dev.cur() != 1) dev.off()



  