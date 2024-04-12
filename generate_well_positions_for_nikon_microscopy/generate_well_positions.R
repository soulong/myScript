
# set path
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |>
  setwd()

# load packages
library(tidyverse)
library(xml2)

# load config
source("config.R")
# load functions
source("functions.R")



## plate info --------------
plate_name <- "CellVis_P96"
fname <- plate_name

# read order priority
image_order <- "row" # column

# choose which objective to use
objective <- "20x"

# 0 means field close next to each other
# 1.5 means each field will interval by 1.5x objective viewfield distance
field_interval <- 0.2
# define field by N x N
N <- 4


if(T) {
  # used setting from config.R
  n_row <- plate_type[[plate_name]]$n_row
  n_col <- plate_type[[plate_name]]$n_col
  well_interval_x <- plate_type[[plate_name]]$well_interval_x
  well_interval_y <- plate_type[[plate_name]]$well_interval_y
  well_diameter <- plate_type[[plate_name]]$well_diameter
  a1_x <- plate_type[[plate_name]]$a1_x
  a1_y <- plate_type[[plate_name]]$a1_y
  pos_z <- plate_type[[plate_name]]$pos_z
  pos_pfs <- plate_type[[plate_name]]$pos_pfs
} else {
  # use defined setting here
  n_row <- 16
  n_col <- 24
  well_interval_x <- 4500 # interval of X [uM]
  well_interval_y <- 4500 # interval of Y [uM]
  well_diameter <- 3200 # well bottom diameter
  a1_x <- 98505 # A1 Position X [uM], for Nikon in lab
  a1_y <- 2435 # A1 Position Y [uM], for Nikon in lab
  pos_z <- 5210 # position Z
  pos_pfs <- 8781 # position PFS
}



## NxN matrix field ---------------
ref_well_fields <- generate_well_field_coords(
  well_diameter, a1_x, a1_y, 
  field_width = objective_view[objective], 
  field_interval = max(0, objective_view[objective] * field_interval),
  N) %>% 
  ## optimize field positions 
  optim_fields()

# check fields
plot_well_fields(c(a1_x, a1_y), 
                 well_diameter, ref_well_fields,
                 font_size = 3)


# > filter fields --------
ref_well_fields <- ref_well_fields %>% 
  filter(field %in% 1:100) %>% 
  filter(field %in% 1:5)

# check fields
plot_well_fields(c(a1_x, a1_y), 
                 well_diameter, ref_well_fields,
                 font_size = 3)


# get all well fields
all_well_fields <- shift_well_coords(
  n_row, n_col, 
  well_interval_x, well_interval_y,
  xtype = -1, 
  ytype = 1) %>% 
  # glimpse()
  crossing(ref_well_fields) %>% 
  transmute(well, field, 
            x = x + x_shift, 
            y = y + y_shift)


## > filter wells ----------------------
if(F) {
  # filter wells method 1
  shiny_select_wells(n_row, n_col)
  interactivate_wells <- 
    paste("selected_wells_", Sys.Date(), ".csv", sep="") %>% 
    read.csv() %>% pull()
  well_fields_filtered <- all_well_fields %>% 
    # separate_wider_delim(strName, "#", 
    #                      names = c("well", "field"), 
    #                      cols_remove = F, 
    #                      too_few = "align_start") %>% 
    filter(well %in% interactivate_wells)
  
} else {
  # filter wells method 2
  well_fields_filtered <- all_well_fields %>% 
    filter(!str_detect(well, "^[A-Z](1)$")) %>% # remove col 1
    filter(!str_detect(well, "^[A-Z](12)$")) %>% # remove col 12
    filter(!str_detect(well, "^[A]")) %>% # remove row A
    filter(!str_detect(well, "^[H]")) %>%  # remove row H
    filter(!str_detect(well, "[CDEF][2-11]G[2-6]$"))
  
}


# optimize well coords
well_fields_filtered_optim <- optim_wells(
  well_fields_filtered, by = image_order)


# add z and pfs info
final <- well_fields_filtered_optim %>% 
  mutate(PositionID = row_number() %>% 
           str_pad(5, "left", "0") %>% str_c("Point", .),
         bChecked = "true",
         strName = str_c(well, "#", field),
         dXPosition = round(x,1),
         dYPosition = round(y,1),
         dZPosition = pos_z,
         dPFSOffset = pos_pfs,
         .before = 1
  )
# check wells
unique(final$well)
print(final)
# write_csv(final, "per_plate.csv")



## > save to xml ---------------
# save_xml(new_points_filtered_optim)
save_xml(final, 
         str_glue("{Sys.Date()}_{fname}.xml"))





