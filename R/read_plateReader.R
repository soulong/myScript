
library(readxl)

read_plateReader <- function(
    f, # file location
    well_type = "96", # plate well counts
    sheet = NULL, # which sheet, NULL will read all sheets
    invert = FALSE # if read plate in vertically and horizontally inverted way
) {
  
  # get data range
  n_row <- case_when(well_type=="24" ~ 4,
                     well_type=="48" ~ 6,
                     well_type=="96" ~ 8,
                     well_type=="384" ~ 16)
  n_col <- as.integer(well_type) / n_row
  
  # get all sheets
  if(is.null(sheet)) sheets <- excel_sheets(f) else sheets <- sheet
  
  res <- list()
  for(sheet in sheets) {
    print(str_glue("read sheet: {sheet}"))
    
    # check where "A" row is 
    start_loc <- tryCatch(
      {read_xlsx(f, sheet, 
                 range = "A1:A100", col_names = F,
                 .name_repair = "minimal") %>% 
          deframe() %>% 
          # "A" must be included in the first column of file
          str_equal("A") %>% 
          which()}, 
      error = function(e) return(NULL)
    )

    if(is.null(start_loc)) {
      print("No 'A' found in first column, check file. skip sheet!")
      next
    }
    
    
    data_range <- 
      str_glue("B{start_loc}:{LETTERS[n_col + 1]}{start_loc + n_row - 1}")
    
    # read data
    raw <- read_xlsx(f, sheet, 
                     range = data_range, col_names = F, 
                     .name_repair = "minimal")
    # invert data
    if(isTRUE(invert)) raw <- raw[nrow(raw):1, ncol(raw):1]
    
    data <- raw %>% 
      as.matrix() %>% 
      as.vector() %>% 
      as.numeric() %>% 
      enframe(name = NULL) %>%
      bind_cols(Metadata_well = map2_chr(
        rep(LETTERS[seq_len(n_row)], n_col), 
        rep(seq_len(n_col), each = n_row), 
        # %>% str_pad(2, pad = "0"),
        ~ str_c(.x, .y)),
        .) %>% 
      filter(!is.na(value))

    # output to list
    res[[sheet]] <- data
  }
  
  res <- list_rbind(res, names_to="Metadata_sheet")
  
  return(res)
}