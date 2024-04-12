

# Microscopy setup
# field of view -> actual length, unit: um
objective_view <- c("100x" = 133, "60x" = 221, "20x" = 600)


# Define plate here
plate_type <- list(
  "Phenoplate_P384" = list(n_row=16, 
                           n_col=24,
                           well_interval_x=4500, # uM
                           well_interval_y=4500, # uM
                           well_diameter=3200, # uM
                           a1_x=99126, # A1 Position X uM
                           a1_y=3641, # A1 Position Y uM
                           pos_z=5512, # defualt z position
                           pos_pfs=11765), # defualt PFS position
  
  "Phenoplate_P96" = list(n_row=8, 
                          n_col=12,
                          well_interval_x=4500, 
                          well_interval_y=4500,
                          well_diameter=3200,
                          a1_x=93275,
                          a1_y=6170,
                          pos_z=5512,
                          pos_pfs=11765), 
  
  "CellVis_P384" = list(n_row=16, 
                        n_col=24,
                        well_interval_x=4500, 
                        well_interval_y=4500,
                        well_diameter=3200,
                        a1_x=93275,
                        a1_y=6170,
                        pos_z=5512,
                        pos_pfs=11765), 
  
  "CellVis_P96" = list(n_row=8, 
                       n_col=12,
                       well_interval_x=9000, 
                       well_interval_y=9000,
                       well_diameter=6000,
                       a1_x=96629,
                       a1_y=4591,
                       pos_z=5260,
                       pos_pfs=8720), 
  
  "CellVis_P24" = list(n_row=4, 
                       n_col=6,
                       well_interval_x=19300, 
                       well_interval_y=19300,
                       well_diameter=15540,
                       a1_x=98505,
                       a1_y=2435,
                       pos_z=5243,
                       pos_pfs=11442),
  
  "CellVis_C18" = list(n_row=3, 
                       n_col=6,
                       well_interval_x=7560, 
                       well_interval_y=7140,
                       well_diameter=5650,
                       a1_x=98505,
                       a1_y=2435,
                       pos_z=5243,
                       pos_pfs=11442),
  
  "CellVis_C8" = list(n_row=2, 
                      n_col=4,
                      well_interval_x=11570, 
                      well_interval_y=11200,
                      well_diameter=8700,
                      a1_x=98505,
                      a1_y=2435,
                      pos_z=5243,
                      pos_pfs=11442)
)

