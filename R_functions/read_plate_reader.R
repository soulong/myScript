

read_plate_reader <- function(f, 
                              sheet = 1, 
                              well_type = "96") {
  # check where "A" row is 
  start_loc <- read_xlsx(f, sheet, 
                         range = "A1:A100", col_names = F,
                         .name_repair = "minimal") %>% 
    deframe() %>% 
    # "A" must be included in the first column of file
    str_equal("A") %>% 
    which()
  if(length(start_loc) == 0) 
    stop("no 'A' found in first column, please check file")
  
  # get data range
  n_row <- case_when(well_type=="24" ~ 4,
                     well_type=="48" ~ 6,
                     well_type=="96" ~ 8,
                     well_type=="384" ~ 16)
  n_col <- as.integer(well_type) / n_row
  data_range <- 
    str_glue("B{start_loc}:{LETTERS[n_col + 1]}{start_loc + n_row - 1}")
  
  # read data
  read_xlsx(f, sheet, 
            range = data_range, col_names = F, 
            .name_repair = "minimal") %>% 
    as.matrix() %>% 
    as.vector() %>% 
    as.numeric() %>% 
    enframe(name = NULL) %>%
    bind_cols(well = map2_chr(rep(LETTERS[seq_len(n_row)], n_col), 
                              rep(seq_len(n_col), each = n_row), 
                              # %>% str_pad(2, pad = "0"),
                              ~ str_c(.x, .y)),
              .) %>% 
    filter(!is.na(value))
}