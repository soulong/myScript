

# used for separated metadata in multiple sheets
# usually, use more files to record different plates
# one info file related to one dataset
read_plate_info <- function(
    file, add_directory=F, colname_prefix="") {
  
  if(!file.exists(file)) return(NULL)
  print(str_glue("reading: {file}"))
  
  sheets <- readxl::excel_sheets(file)
  print(str_glue("expected sheets: {str_c(sheets, collapse=', ')}"))
  
  plate_info <- map(sheets, \(x) readxl::read_xlsx(
    file, sheet = x, range = "B2:Y17",  # A1:X16
    col_names = F, col_types = "text", 
    .name_repair = "minimal") %>%
      flatten_chr() %>%
      enframe(name = NULL, value = x) ) %>%
    bind_cols() %>%
    bind_cols(well=map2_chr(
      rep(LETTERS[1:16], 24), 
      rep(1:24, each = 16), # %>% str_pad(2, pad = "0"), 
      ~ str_c(.x, .y)), .) %>% 
    filter(!if_all(!well, is.na)) # remove all empty rows
  
  if(add_directory) {
    plate_info <- plate_info %>%
      mutate(directory=file %>% 
               str_split('[/]', simplify=T) %>% 
               str_subset('.*-Measurement [0-9]') %>%
               str_split('(__)', simplify=T) %>%
               .[1], .before=1)
  }
  
  plate_info <- plate_info %>% 
    rename_with(~ str_c(colname_prefix, .x), everything())
  
  return(plate_info)
}


# used for merged metadata for multiple sheets
# different sheets in one info file correspond to different dataset
# in this format, different metadata should be written into one sheet
read_plate_info_2 <- function(file) {
  
  if(!file.exists(file)) return(NULL)
  print(str_glue("processing:\n {file}"))
  
  sheets <- readxl::excel_sheets(file) %>% 
    set_names(., nm=.)
  print(str_glue("expected sheets: {str_c(sheets, collapse=', ')}"))
  
  plate_info <- sheets %>% 
    map(\(x) readxl::read_xlsx(
      file, sheet = x, range = "B2:Y17",  # A1:X16
      col_names = F, col_types = "text", 
      .name_repair = "minimal") %>%
        as.matrix() %>% 
        as.vector() %>% 
        enframe(name = NULL) %>%
        bind_cols(Metadata_well = map2_chr(
          rep(LETTERS[1:16], 24), rep(1:24, each = 16), 
          \(x, y) str_c(x, y) ),
          .) %>% 
        filter(!is.na(value)) ) %>% 
    bind_rows(.id = "Metadata_sheet")
  
  return(plate_info)
}




