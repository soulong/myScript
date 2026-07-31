
library(rio)
options(rio.import.class='tbl')
library(scales)
library(ggsci)
library(rlang)
library(transport)
library(energy)
# library(GRmetrics)
require(bestNormalize)
library(tidyverse)
library(furrr)

# **********************************************************************
# general func -----------------------------------------------------------
# **********************************************************************

se <- function(x) {sd(x, na.rm=T)/sqrt(length(x)-1)}


#' Format tidy well style
#'
#' @param wells a character vector of mixed style wells
#'
#' @returns tidy wells
#' @export
#'
#' @examples
#' transform_well_style(c('A1', 'B01', 'C003', 'D5'))
transform_well_style <- function(wells) {
  sub("^([A-Za-z]+)0+", "\\1", wells)
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
#' 
#' @export
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




list_cbind_efficient <- function(x, y, by=c("uid")) {
  library(data.table)
  
  x <- as.data.table(x)
  y <- as.data.table(y)
  
  common <- setdiff(intersect(names(x), names(y)), by)
  
  if(length(common)) {
    setnames(y, common, paste0(common, ".y"))
  }
  
  y[x, on=by]
}


# **********************************************************************
# read sqlite ------------------------------------------------------
# **********************************************************************

read_sqlite <- function(
    f, skip_tables=c('metadata', 'reduction_pca_variance', 'sqlite_sequence')
) {
  print(paste0('read table: ', f))
  
  conn <- DBI::dbConnect(RSQLite::SQLite(), dbname=f)
  tables <- DBI::dbListTables(conn) %>% 
    setdiff(skip_tables) %>% 
    set_names(., nm=.)
  
  res <- tryCatch(
    map(tables, 
        \(x) tbl(conn, x) %>% collect() %>% 
          janitor::remove_constant(quiet=F)
    ),
    error = function(e) print(e),
    finally = DBI::dbDisconnect(conn)
  )
  
  return(res)
}



# **********************************************************************
# read excel metadata ------------------------------------------------------
# **********************************************************************

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
  read_sheet <- \(sheet) readxl::read_xlsx(
    file, sheet = sheet,
    range = "B2:Y17",
    col_names = FALSE, col_types = "text",
    .name_repair = "minimal"
  ) %>%
    as.matrix() %>%
    as.vector()
  
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




# **********************************************************************
# ggplot2 related ------------------------------------------------------
# **********************************************************************

#' Adaptive color palette generator
#'
#' Adaptive color palette generator for ggsci color palettes using `pal_ramp()`.
#'
#' @param name Color palette name in ggsci
#' @param palette Color palette type in ggsci
#' @param alpha Transparency level, a real number in (0, 1].
#'
#' @details See `names(ggsci:::ggsci_db)` for all color palette names in ggsci.
#' See `names(ggsci:::ggsci_db$"pal")` for available palette types under
#' the palette `pal`.
pal_adaptive <- function(name, palette, alpha = 1) {
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")
  
  raw_cols <- ggsci:::ggsci_db[[name]][[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )
  
  pal_ramp(unname(alpha_cols))
}


#' Adaptive color scales
#'
#' @inheritParams pal_adaptive
#' @param ... additional parameters for [ggplot2::discrete_scale()].
scale_color_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("colour", name, pal_adaptive(name, palette, alpha), ...)
}
scale_fill_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("fill", name, pal_adaptive(name, palette, alpha), ...)
}


## add new geom_xx at Z-order
`-.gg` <- function(plot, layer) {
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}





# adjust x/y axis
# https://github.com/jbengler/tidyplots/blob/master/R/adjust.R#L2
rotate_labels <- function(
    plot,
    axis = "x", 
    angle = 45) {
  
  if (angle == TRUE) angle <- 45
  if (is.numeric(angle) && angle != 0) {
    if(angle >= 90) vjust <- 0.5 else vjust <- 1
    if(axis == "x")
      plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = 1, vjust = vjust))
    if(axis == "y")
      plot <- plot + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = angle, hjust = vjust, vjust = 1))
  }
  return(plot)
}


my_theme <- function(angle=0,
                     fontsize = 7, 
                     family = "sans", 
                     face = NULL, 
                     color = "black",
                     legend_size = 4,
                     unit = "mm") {
  
  # # if(angle == TRUE) angle <- 45
  # if(is.numeric(angle) && angle != 0) {
  #   if(angle >= 90) vjust <- 0.5 else vjust <- 1
  #   if(axis == "x")
  #     plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = 1, vjust = vjust))
  #   if(axis == "y")
  #     plot <- plot + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = angle, hjust = vjust, vjust = 1))
  # }
  
  if(is.numeric(angle) && angle != 0) {
    # for x-axis
    if(angle >= 90) vjust <- 0.5 else vjust <- 1
    th <- theme_classic() %+replace%
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = 1, vjust = vjust))
  } else { th <- theme_classic() }
  
  th %+replace%
    ggplot2::theme(
      text = ggplot2::element_text(size = fontsize, family = family, face = face, color = color),
      plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = fontsize/2, hjust = 0.5, vjust = 0.5),
      plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), unit),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      # legend.title = ggplot2::element_text(size = fontsize/2),
      # legend.text = ggplot2::element_text(size = fontsize/2),
      legend.background = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key.size = ggplot2::unit(legend_size, unit),
      # strip.text = ggplot2::element_text(size = fontsize, family = family, face = face, colour = color),
      strip.background = ggplot2::element_rect(fill = NA, colour = NA),
      # axis.text = ggplot2::element_text(size = fontsize, family = family, face = face, colour = color),
      # axis.title = ggplot2::element_text(size = fontsize, family = family, face = face, colour = color),
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.25)
    )
}





# format numeric to string, to keep length consistent
format_num_auto <- function(v, accuracy=0.001) {
  
  if(typeof(v) == "character") v <- as.numeric(v)
  
  map(v, \(x) if(abs(x) < accuracy) {
    scales::label_scientific()(x) 
  } else {
    round(x, log10(1/accuracy))
  }
  ) %>% as.character() 
}



#' Add an automatic log10 scale to a ggplot
#'
#' Automatically detects the variable mapped to the selected axis from
#' ggplot2 aesthetics, calculates positive data limits, and generates
#' logarithmic breaks.
#'
#' @param axis Character. Which axis to transform. Must be either
#'   `"x"` or `"y"`. Default is `"x"`.
#'
#' @param limits Numeric vector of length two. Optional scale limits.
#'   If NULL (default), limits are automatically calculated from the
#'   mapped variable in ggplot aesthetics.
#'
#' @param breaks Numeric vector. Positions of log-scale breaks.
#'   If NULL (default), decade breaks are generated automatically.
#'
#' @param base Numeric. Logarithm base. Default is 10.
#'
#' @param labels Function or character vector. Labels for scale ticks.
#'   If NULL, labels are generated using scientific notation
#'   (e.g. 10^0, 10^1, 10^2).
#'
#' @param ... Additional arguments passed to `scale_x_log10()` or
#'   `scale_y_log10()`.
#'
#' @return A ggplot2 scale object.
#'
#' @examples
#' ggplot(df, aes(conc, response)) +
#'   geom_point() +
#'   scale_log(axis = "x")
#'
#' @export
scale_log <- function(axis = c("x", "y"),
                      limits = NULL,
                      breaks = NULL,
                      base = 10,
                      labels = NULL,
                      ...) {
  
  axis <- match.arg(axis)
  
  structure(
    list(
      axis = axis,
      limits = limits,
      breaks = breaks,
      base = base,
      labels = labels,
      params = list(...)
    ),
    class = "scale_log_auto"
  )
}

ggplot_add.scale_log_auto <- function(object, plot, object_name) {
  
  axis <- object$axis
  
  # extract mapped variable from aes()
  aes_var <- if (axis == "x") {
    plot$mapping$x
  } else {
    plot$mapping$y
  }
  
  if (is.null(aes_var)) {
    stop("Cannot detect ", axis, " aesthetic from aes()")
  }
  
  var <- as_name(aes_var)
  
  
  # obtain plotting data
  data <- plot$data
  
  if (is.null(data)) {
    # fallback to first layer data
    data <- plot$layers[[1]]$data
  }
  
  if (is.null(data)) {
    stop("Cannot find plot data")
  }
  
  
  # auto range
  if (is.null(object$limits)) {
    
    x <- data[[var]]
    
    x <- x[x > 0 & is.finite(x)]
    
    object$limits <- range(x, na.rm = TRUE)
  }
  
  
  # auto breaks
  if (is.null(object$breaks)) {
    
    rng <- object$limits
    
    object$breaks <-
      object$base^(
        floor(log(rng[1], object$base)):
          ceiling(log(rng[2], object$base))
      )
  }
  
  
  if (is.null(object$labels)) {
    object$labels <- label_math(object$base^.x)
  }
  
  
  scale <- if (axis == "x") {
    
    scale_x_log10(
      limits = object$limits,
      breaks = object$breaks,
      labels = object$labels,
      ...
    )
    
  } else {
    
    scale_y_log10(
      limits = object$limits,
      breaks = object$breaks,
      labels = object$labels,
      ...
    )
  }
  
  
  ggplot_add(scale, plot, object_name)
}



# Polygon-based cell selector for ggplot objects
# Returns indices of selected points within the polygon
ggplot_selector <- function(plot_obj) {
  # Extract data from ggplot object
  if (!"ggplot" %in% class(plot_obj)) {
    stop("Input must be a ggplot object")
  }
  
  # Get the built plot data
  built_plot <- ggplot_build(plot_obj)
  plot_data <- built_plot$data[[1]]
  
  # Extract x and y coordinates
  x_coords <- plot_data$x
  y_coords <- plot_data$y
  
  if (length(x_coords) == 0 || length(y_coords) == 0) {
    stop("No data points found in the plot")
  }
  
  # Ensure graphics device is active and plot is displayed
  if (dev.cur() == 1) {  # If no graphics device is open
    dev.new()
  }
  
  # Extract additional aesthetic information from ggplot
  plot_colors <- plot_data$colour %||% plot_data$fill %||% "black"
  plot_sizes <- plot_data$size %||% 1
  plot_shapes <- plot_data$shape %||% 19
  
  # Convert ggplot colors to R colors if needed
  if (is.factor(plot_colors)) {
    # If colors are mapped to factors, create a color palette
    unique_colors <- unique(plot_colors)
    color_palette <- rainbow(length(unique_colors))
    names(color_palette) <- unique_colors
    plot_colors <- color_palette[as.character(plot_colors)]
  }
  
  # Create base R plot with retained aesthetics
  plot(x_coords, y_coords, 
       xlab = plot_obj$labels$x %||% "x", 
       ylab = plot_obj$labels$y %||% "y",
       main = plot_obj$labels$title %||% "",
       col = plot_colors,
       pch = plot_shapes,
       cex = plot_sizes * 0.5)  # Scale down sizes for base R
  
  cat("Instructions:\n")
  cat("1. Click points to define polygon vertices\n")
  cat("2. Right-click or press ESC when done selecting\n")
  cat("3. The polygon will be automatically closed\n\n")
  
  # Interactive polygon selection
  polygon_points <- locator(type = "o", col = "red", lwd = 2)
  
  if (is.null(polygon_points) || length(polygon_points$x) < 3) {
    cat("No valid polygon selected (need at least 3 points)\n")
    return(integer(0))
  }
  
  # Close the polygon by connecting last point to first
  polygon_x <- c(polygon_points$x, polygon_points$x[1])
  polygon_y <- c(polygon_points$y, polygon_points$y[1])
  
  # Draw the final polygon
  lines(polygon_x, polygon_y, col = "red", lwd = 2)
  
  # Point-in-polygon test using ray casting algorithm
  point_in_polygon <- function(px, py, poly_x, poly_y) {
    n <- length(poly_x) - 1  # Exclude the closing point
    inside <- FALSE
    
    j <- n
    for (i in 1:n) {
      if (((poly_y[i] > py) != (poly_y[j] > py)) &&
          (px < (poly_x[j] - poly_x[i]) * (py - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i])) {
        inside <- !inside
      }
      j <- i
    }
    return(inside)
  }
  
  # Test each point
  selected_indices <- c()
  for (i in seq_along(x_coords)) {
    if (point_in_polygon(x_coords[i], y_coords[i], polygon_x, polygon_y)) {
      selected_indices <- c(selected_indices, i)
    }
  }
  
  # Highlight selected points
  if (length(selected_indices) > 0) {
    points(x_coords[selected_indices], y_coords[selected_indices], 
           col = "blue", pch = 19, cex = 0.8)
  }
  
  cat(sprintf("Selected %d points\n", length(selected_indices)))
  
  return(selected_indices)
}




# **********************************************************************
# dose response ------------------------------------------------------
# **********************************************************************

# func
model_func <- function(data, ...) {
  res <- try( 
    dr4pl(response ~ dose, 
          data = data, 
          # method.init = "logistic",
          # method.robust = "squared",
          # lowerl = c(theta_4 = 0),
          ...))
  if(class(res) == "try-error") res <- NULL
  return(res)
}

coef_tidy <- function(model) {
  if(class(model) != "dr4pl") return(NULL)
  
  ci <- summary(model) %>% 
    .[["coefficients"]] %>%
    as.data.frame()
  ci <- mutate(ci, CI = map_dbl(
    seq_along(ci), ~ (ci[[4]][.x] - ci[[3]][.x])/2))
  # convert logIC50 to IC50
  ci_ic50 <- ci[2, ] %>% 10^.
  rownames(ci_ic50) <- "IC50"
  # combine
  out <- bind_rows(ci[-2, ], ci_ic50) %>% 
    .[c(1,4,2,3), c(1,5)] %>%
    mutate(across(everything(), \(x) format(x, digits=3, scientific=T))) %>%
    mutate(across(everything(), str_trim)) %>%
    # mutate(across(everything(), as.numeric)) %>%
    as_tibble(rownames = "type") %>%
    mutate(type = str_replace_all(type, "Limit", "")) %>%
    pivot_wider(names_from = type, 
                names_glue = "{type}_{.value}", 
                values_from = c(Estimate, CI))
  return(out)
}

pred_func <- function(model, se=F, normalize_residual=T, level=0.95, nboot=200) {
  if(class(model) != "dr4pl") return(NULL)
  
  # model = res$model[[2]]
  range <- sort(unique(model$data$Dose)) %>% na.omit()
  from <- ifelse(range[1]<=0, 0.8*range[2], 0.8*range[1])
  to <- 1.2*last(range)
  xseq <- exp(seq(from=log(from), to=log(to), length.out=200))
  
  # add seq from 0 to min value
  if(range[1]<=0) {
    xseq <- c(xseq,
              exp(seq(from=log(0 + range[2]/1000), to=log(range[2]), length.out=50)))
  }
  
  pred <- MeanResponse(model$parameters, xseq)
  
  if (!se) {
    return(base::data.frame(x=xseq, y=pred))
  }
  
  ## bootstrap residuals
  # use raw residuals
  pred0 <- MeanResponse(model$parameters, model$data$Dose)
  res <- pred0 - model$data$Response
  if(normalize_residual) {
    res <- bestNormalize::yeojohnson(res, standardize = F)$x.t %>%
      scale(center = T,  scale = F) %>%
      as.vector()
  }
  
  # parallel worker
  bootres <- suppressMessages(future_map_dfc(seq(nboot), function(x) {
    response_resample <- pred0 + sample(res, size=length(pred0), replace=TRUE)
    mboot <- try(dr4pl(model$data$Dose, response_resample))
    
    if(class(mboot)=="try-error") {
      return(NA)
    }  else { return(MeanResponse(mboot$parameters, xseq)) 
    } }, .progress = F, .options = furrr_options(seed=NULL))) %>%
    as.matrix()
  
  # pb <- txtProgressBar(max=nboot, style=3)
  # for (i in seq(nboot)) {
  #   setTxtProgressBar(pb, i)
  #   # mboot <- dr4pl(model$data$Dose, pred0 + sample(res, size=length(pred0), replace=TRUE))
  #   # bootres[, i] <- MeanResponse(mboot$parameters, xseq)
  # }
  
  fit <- data.frame(x = xseq,
                    y = pred,
                    ymin = apply(bootres, 1, quantile, probs=(1-level)/2, na.rm=T),
                    ymax = apply(bootres, 1, quantile, probs=(1+level)/2, na.rm=T) )
  return(fit)
}


fit_dose_response <- function(
    df, dose, response, group=NULL, n_worker=8, ...) {
  
  library(furrr)
  library(dr4pl)
  library(bestNormalize)
  
  if(n_worker > 1) plan(multisession, workers = n_worker)
  
  fitted <- df %>% 
    select(dose = !!as.name(dose), 
           response = !!as.name(response),
           any_of(group)) %>% 
    mutate(dose=as.numeric(as.character(dose)), 
           response=as.numeric(as.character(response))) %>% 
    group_by(across(any_of(group))) %>%
    nest() %>% 
    # .[1:10, ] %>% 
    mutate(model = future_map(data, model_func, ...))
  # glimpse(fitted)
  
  plan(sequential)
  
  # get more expansion data
  fitted <- fitted %>%
    mutate(
      convergence = map_lgl(model, \(x) if(is.null(x)) F else x[["convergence"]]),
      method_robust = map_chr(model, \(x) if(is.null(x)) "None" else x[["method.robust"]]),
      coef = map(model, coef_tidy),
      pred = map(model, pred_func)
    ) %>% 
    ungroup()
  
  return(fitted)
}


plot_dose_response_dataPrepare <- function(fitted, group=NULL) {
  # calculate mean and se
  fitted_mean <- fitted %>% 
    # filter(convergence) %>% 
    # filter(!method_robust=="None") %>% 
    mutate(stats = map(
      data, 
      ~ group_by(.x, dose) %>% 
        summarise(mean=mean(response), se=sd(response)/sqrt(n()), .groups = "drop"))) %>%
    mutate(top_inhit = map(pred, function(x) round(100*(1-min(x$y)/max(x$y)), 0))) %>%
    unnest(coef, keep_empty=T) %>%
    {if(!is.null(group)) unite(., "uid", !!!syms(group), sep="\n", remove = F) else mutate(., uid=1) }%>%
    mutate(group_label = str_glue("{uid}\n IC50: {str_replace(IC50_Estimate,'[+]','')} ± {str_replace(IC50_CI,'[+]','')}\n Top Inhibition: {top_inhit}%"))
  
  return(fitted_mean)
}






# **********************************************************************
# population distance ------------------------------------------------------
# **********************************************************************

#' Parallel Population Distance Calculation
#'
#' @param mat matrix or data.frame of features/coordinates (cells/samples x features)
#' @param group vector or factor assigning each row to an experimental condition/group
#' @param method choice of "mmd", "edistance", "wasserstein", or "mahalanobis"
#' @param ref optional character vector of reference group names. If provided, calculates 
#'            distance only relative to these references and returns a tidy tibble. 
#'            If NULL, returns all pairwise distances as a stats::dist object.
#' @param sigma optional numeric; bandwidth parameter for Gaussian RBF kernel in MMD.
#' @param max_cells integer; maximum observations per group to subsample. Default is 2000.
#' @param workers integer; number of CPU cores to use. Default is half of available cores.
#' @return A `tibble` if `ref` is provided, otherwise a `stats::dist` object.
calc_distance <- function(mat, group, 
                          method = c("mmd", "edistance", "wasserstein", "mahalanobis"),
                          ref = NULL,
                          sigma = NULL,
                          max_cells = 2000,
                          workers = floor(future::availableCores() / 2)) {
  
  method <- match.arg(method)
  mat <- as.matrix(mat)
  group_vec <- as.character(group)
  all_groups <- unique(group_vec)
  
  if (!is.null(ref)) {
    ref <- as.character(ref)
    if (!all(ref %in% all_groups)) {
      stop("All elements in 'ref' must exist in unique values of 'group'.")
    }
  }
  
  # Split data once by group
  group_splits <- split(as.data.frame(mat), group_vec)
  
  # Pre-compute pooled covariance if Mahalanobis is requested
  if (method == "mahalanobis") {
    means <- t(sapply(group_splits, colMeans))
    pooled_cov <- stats::cov(mat)
    inv_cov <- solve(pooled_cov)
  }
  
  # Helper for MMD computation
  compute_mmd <- function(X, Y, bandwidth) {
    nX <- nrow(X)
    nY <- nrow(Y)
    
    XY_comb <- rbind(X, Y)
    D2 <- as.matrix(stats::dist(XY_comb))^2
    
    if (is.null(bandwidth)) {
      non_zero_dists <- D2[D2 > 0]
      bandwidth <- sqrt(0.5 * median(non_zero_dists))
    }
    
    K <- exp(-D2 / (2 * bandwidth^2))
    
    K_xx <- K[1:nX, 1:nX]
    K_yy <- K[(nX + 1):(nX + nY), (nX + 1):(nX + nY)]
    K_xy <- K[1:nX, (nX + 1):(nX + nY)]
    
    mmd2_xx <- (sum(K_xx) - nX) / (nX * (nX - 1))
    mmd2_yy <- (sum(K_yy) - nY) / (nY * (nY - 1))
    mmd2_xy <- mean(K_xy)
    
    mmd2 <- mmd2_xx + mmd2_yy - 2 * mmd2_xy
    return(sqrt(max(0, mmd2)))
  }
  
  # Worker kernel execution
  calc_single_dist <- function(g1_name, g2_name) {
    if (g1_name == g2_name) return(0)
    
    if (method == "mahalanobis") {
      diff <- means[g1_name, ] - means[g2_name, ]
      return(sqrt(as.numeric(t(diff) %*% inv_cov %*% diff)))
    }
    
    g1 <- as.matrix(group_splits[[g1_name]])
    g2 <- as.matrix(group_splits[[g2_name]])
    
    if (nrow(g1) > max_cells) g1 <- g1[sample(seq_len(nrow(g1)), max_cells), ]
    if (nrow(g2) > max_cells) g2 <- g2[sample(seq_len(nrow(g2)), max_cells), ]
    
    switch(
      method,
      "mmd" = compute_mmd(g1, g2, bandwidth = sigma),
      "edistance" = energy::edist(rbind(g1, g2), c(nrow(g1), nrow(g2)))[1, 2],
      "wasserstein" = transport::wasserstein(transport::pp(g1), transport::pp(g2), p = 2)
    )
  }
  
  # --- OS-Agnostic Multiprocessing Setup ---
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::multisession, workers = workers)
  }
  
  # --- BRANCH 1: Target Reference Distances (Tidy Tibble Output) ---
  if (!is.null(ref)) {
    pairs <- expand.grid(group_1 = ref, group_2 = all_groups, stringsAsFactors = FALSE)
    
    distances <- furrr::future_map2_dbl(
      pairs$group_1, 
      pairs$group_2, 
      ~calc_single_dist(.x, .y),
      .options = furrr::furrr_options(seed = TRUE)
    )
    
    pairs$distance <- distances
    return(tibble::as_tibble(pairs))
  }
  
  # --- BRANCH 2: All Pairwise Distances (dist object Output) ---
  n_groups <- length(all_groups)
  
  idx <- which(lower.tri(matrix(0, n_groups, n_groups)), arr.ind = TRUE)
  g1_vec <- all_groups[idx[, 1]]
  g2_vec <- all_groups[idx[, 2]]
  
  dist_vals <- furrr::future_map2_dbl(
    g1_vec, 
    g2_vec, 
    ~calc_single_dist(.x, .y),
    .options = furrr::furrr_options(seed = TRUE)
  )
  
  dist_mat <- matrix(0, nrow = n_groups, ncol = n_groups, 
                     dimnames = list(all_groups, all_groups))
  
  for (k in seq_along(dist_vals)) {
    i_name <- g1_vec[k]
    j_name <- g2_vec[k]
    dist_mat[i_name, j_name] <- dist_vals[k]
    dist_mat[j_name, i_name] <- dist_vals[k]
  }
  
  return(as.dist(dist_mat))
}