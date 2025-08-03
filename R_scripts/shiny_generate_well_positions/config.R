
# Microscopy setup

# Define common config here

# for shulab
common_config_shulab <- c(
  # field of view -> actual length, unit: um
  # use DDDx format, to easy extract mag info
  "100x" = 133, "60x" = 221, "20x" = 600, 
  # coord shift direction
  x_direct=-1, y_direct=1
  )

# for CVRI RM384 nikon spining disk
common_config_rm384 <- c(
  # field of view -> actual length, unit: um
  # use DDDx format, to easy extract mag info
  "100x" = 97, "60x" = 162, "40x" = 220, "20x" = 440, "10x" = 900,
  # coord shift direction
  x_direct=1, y_direct=-1
  )


# Define specific plate here
plate_type <- list(
  "CellVis_P96_RM384" = list(n_row=8, 
                             n_col=12,
                             well_interval_x=9000, 
                             well_interval_y=9000,
                             well_diameter=6000,
                             a1_x=-35062,
                             a1_y=35202,
                             pos_z=5334,
                             pos_pfs=7200) %>%
    c(., common_config_rm384), 
  
  "CellVis_P24_RM384" = list(n_row=4, 
                             n_col=6,
                             well_interval_x=19300, 
                             well_interval_y=19300,
                             well_diameter=15540,
                             a1_x=-32748,
                             a1_y=32633,
                             pos_z=4284,
                             pos_pfs=7200) %>% 
    c(., common_config_rm384),
  
  "Phenoplate_P384" = list(n_row=16, 
                           n_col=24,
                           well_interval_x=4500, # uM
                           well_interval_y=4500, # uM
                           well_diameter=3200, # uM
                           a1_x=98554, # A1 Position X uM
                           a1_y=3966, # A1 Position Y uM
                           pos_z=3721, # default z position
                           pos_pfs=9020 # default PFS position
                           ) %>%
    c(., common_config_shulab), 
  
  "Phenoplate_P96" = list(n_row=8, 
                          n_col=12,
                          well_interval_x=9000, 
                          well_interval_y=9000,
                          well_diameter=6000,
                          a1_x=96412,
                          a1_y=6284,
                          pos_z=3788,
                          pos_pfs=9020) %>%
    c(., common_config_shulab), 
  
  "CellVis_P384" = list(n_row=16, 
                        n_col=24,
                        well_interval_x=4500, 
                        well_interval_y=4500,
                        well_diameter=3200,
                        a1_x=93275,
                        a1_y=6170,
                        pos_z=5512,
                        pos_pfs=11765) %>%
    c(., common_config_shulab), 
  
  "CellVis_P96" = list(n_row=8, 
                       n_col=12,
                       well_interval_x=9000, 
                       well_interval_y=9000,
                       well_diameter=6000,
                       a1_x=96290,
                       a1_y=5089,
                       pos_z=5334,
                       pos_pfs=8720) %>%
    c(., common_config_shulab), 
  
  "CellVis_P24" = list(n_row=4, 
                       n_col=6,
                       well_interval_x=19300, 
                       well_interval_y=19300,
                       well_diameter=15540,
                       a1_x=93586,
                       a1_y=7410,
                       pos_z=4284,
                       pos_pfs=8650) %>%
    c(., common_config_shulab),

  "CellVis_C18" = list(n_row=3, 
                       n_col=6,
                       well_interval_x=7560, 
                       well_interval_y=7140,
                       well_diameter=5650,
                       a1_x=98195,
                       a1_y=28913,
                       pos_z=3674,
                       pos_pfs=8738) %>%
    c(., common_config_shulab),
  
  "CellVis_C8" = list(n_row=2, 
                      n_col=4,
                      well_interval_x=11570, 
                      well_interval_y=11200,
                      well_diameter=8700,
                      a1_x=95795,
                      a1_y=31775,
                      pos_z=3683,
                      pos_pfs=8675) %>%
    c(., common_config_shulab)
)

