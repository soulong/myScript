
#' Detect Operating System
#'
#' @description Detects the current operating system, including special
#'   handling for Windows Subsystem for Linux (WSL).
#'
#' @return Character string indicating the OS type:
#'   \itemize{
#'     \item "wsl" - Windows Subsystem for Linux
#'     \item "windows" - Native Windows
#'     \item "linux" - Native Linux
#'     \item Other - Other OS name from Sys.info()
#'   }
#'
#' @examples
#' detect_os()
#'
#' @export
detect_os <- function() {
  sysinfo <- Sys.info()
  
  if (grepl("microsoft|Microsoft", sysinfo["release"], ignore.case = TRUE)) {
    return("wsl")
  }
  
  if (.Platform$OS.type == "windows") {
    return("windows")
  }
  
  if (sysinfo["sysname"] == "Linux") {
    return("linux")
  }
  
  return(sysinfo["sysname"])
}



#' @Title Safe Path Normalization
#' @description Normalize file paths and convert Windows paths to WSL format when needed
#' @param path Character; the file path to process
#' @param standardize Logical; whether to convert Windows paths to WSL format (default: TRUE)
#' @param os Character; the operating system (default: auto-detect)
#' @return Character string with normalized path
norm_path <- function(path=r"(C:\Users\haohe\GitHub\shinyLab)", 
                      standardize=TRUE, 
                      os=detect_os()
                      ) {
  path <- path |> 
    gsub("\\\\", "/", x = _) |> 
    gsub("/+", "/", x = _) |> 
    trimws()
  
  if(standardize & grepl("^[A-Za-z]:[/\\\\]", path) & os == "wsl") {
    path <- file.path(
      "/mnt", 
      tolower(sub("^([A-Za-z]):.*", "\\1", path)), 
      gsub("\\\\", "/", sub("^[A-Za-z]:[/\\\\]", "", path)))
  }
  
  return(path)
}



#' Transform Well Position Style
#'
#' @description Removes leading zeros from well position strings.
#'   Converts formats like "A01" to "A1".
#'
#' @param wells Character vector of well positions (e.g., "A01", "B12").
#'
#' @return Character vector with leading zeros removed from numeric portion.
#'
#' @examples
#' transform_well_style(c("A01", "B12", "C03"))
#' # Returns: c("A1", "B12", "C3")
#'
#' @export
transform_well_style <- function(wells) sub("^([A-Za-z]+)0+", "\\1", wells)




#' Read Plate Metadata from Excel File
#'
#' @description Reads 384-well plate metadata from Excel files with specific
#'   format (sheets containing data in B2:Y17 range). Supports both wide and
#'   long format output.
#'
#' @param file Path to the Excel file.
#' @param format Output format: "wide" (one row per well, columns per sheet)
#'   or "long" (one row per well-sheet combination).
#' @param add_directory Logical. If TRUE, adds a directory column extracted
#'   from the file path.
#' @param colname_prefix String prefix to add to all column names.
#'
#' @return A tibble containing plate metadata. Returns NULL if file doesn't exist.
#'   For wide format: columns are well, sheet names, and optionally directory.
#'   For long format: columns are Metadata_sheet, Metadata_well, and value.
#'
#' @examples
#' \dontrun{
#' read_metadata("plate_info.xlsx", format = "wide")
#' read_metadata("plate_info.xlsx", format = "long", add_directory = TRUE)
#' }
#'
#' @export
read_metadata <- function(
    file,
    format = c("wide", "long"),
    add_directory = FALSE,
    colname_prefix = "") {
  
  format <- match.arg(format)
  
  if (!file.exists(file)) return(NULL)
  print(str_glue("reading: {file}"))
  
  sheets <- readxl::excel_sheets(file) %>% set_names(., nm = .)
  print(str_glue("expected sheets: {str_c(sheets, collapse=', ')}"))
  
  # shared: read each sheet into a named flat vector
  read_sheet <- function(sheet) {
    readxl::read_xlsx(
      file, sheet = sheet,
      range = "B2:Y17",
      col_names = FALSE, col_types = "text",
      .name_repair = "minimal"
    ) %>%
      as.matrix() %>%
      as.vector()
  }
  
  wells <- map2_chr(
    rep(LETTERS[1:16], 24),
    rep(1:24, each = 16),
    ~ str_c(.x, .y)
  )
  
  if (format == "wide") {
    
    plate_info <- sheets %>%
      map(~ enframe(read_sheet(.x), name = NULL, value = .x)) %>%
      bind_cols() %>%
      bind_cols(well = wells, .) %>%
      filter(!if_all(!well, is.na))
    
    if (add_directory) {
      plate_info <- plate_info %>%
        mutate(directory = file %>%
                 str_split("[/]", simplify = TRUE) %>%
                 str_subset(".*-Measurement [0-9]") %>%
                 str_split("(__)", simplify = TRUE) %>%
                 .[1], .before = 1)
    }
    
    if (nchar(colname_prefix) > 0) {
      plate_info <- plate_info %>%
        rename_with(~ str_c(colname_prefix, .x), everything())
    }
    
  } else { # long
    
    plate_info <- sheets %>%
      map(~ enframe(read_sheet(.x), name = NULL, value = "value") %>%
            bind_cols(Metadata_well = wells, .) %>%
            filter(!is.na(value))
      ) %>%
      bind_rows(.id = "Metadata_sheet")
    
  }
  
  return(plate_info)
}




#' Find Best Match from Pattern List
#'
#' @description Searches through a vector of choices using a list of patterns,
#'   returning the first match found. Useful for auto-selecting reference values
#'   from user data.
#'
#' @param choices Character vector of available choices.
#' @param patterns Character vector of patterns to search for (in order of priority).
#'
#' @return The first choice matching any pattern, or the first choice if no
#'   matches are found. Returns empty string if choices is empty.
#'
#' @examples
#' find_best_match(c("DMSO", "PBS", "Water"), c("blank", "dmso"))
#' # Returns: "DMSO"
#'
#' find_best_match(c("Sample1", "Sample2"), c("blank", "control"))
#' # Returns: "Sample1" (no match, returns first choice)
#'
#' @export
find_best_match <- function(choices, patterns) {
  if (length(choices) == 0) return("")
  
  choices_lower <- tolower(choices)
  for (pattern in patterns) {
    matches <- choices[grepl(pattern, choices_lower, ignore.case = TRUE)]
    if (length(matches) > 0) {
      return(matches[1])
    }
  }
  return(choices[1])
}




