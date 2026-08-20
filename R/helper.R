
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

#' Find the Best Match from a Pattern List
#'
#' Searches a vector of choices using a list of patterns (in order of
#' priority) and returns the first choice matching any pattern. Useful for
#' auto-selecting reference values (e.g. blank/control samples) from user data.
#'
#' @param choices Character vector of available choices.
#' @param patterns Character vector of (regex) patterns, searched in order of priority.
#'
#' @return The first choice matching any pattern; the first choice if nothing
#'   matches; `""` if `choices` is empty.
#'
#' @examples
#' find_best_match(c("DMSO", "PBS", "Water"), c("blank", "dmso"))
#' find_best_match(c("Sample1", "Sample2"), c("blank", "control"))
#'
#' @export
find_best_match <- function(choices, patterns) {
  if (!length(choices)) return("")
  lower <- tolower(choices)
  for (pat in patterns) {
    hit <- choices[grepl(pat, lower, ignore.case = TRUE)]
    if (length(hit)) return(hit[1])
  }
  choices[1]
}


#' Detect Operating System
#'
#' Detects the current operating system, with special handling for Windows
#' Subsystem for Linux (WSL).
#'
#' @return Character string: "wsl", "windows", "linux", or the raw sysname.
#'
#' @examples
#' detect_os()
#'
#' @export
detect_os <- function() {
  sys <- Sys.info()
  if (grepl("microsoft", sys[["release"]], ignore.case = TRUE)) return("wsl")
  if (.Platform$OS.type == "windows") return("windows")
  if (sys[["sysname"]] == "Linux") return("linux")
  sys[["sysname"]]
}


#' Normalize a File Path
#'
#' Normalizes path separators, collapses repeated slashes and, when running
#' under WSL, converts Windows drive paths (e.g. `C:\foo\bar`) to
#' `/mnt/c/foo/bar`.
#'
#' @param path Character; file path to normalize.
#' @param standardize Logical; whether to convert Windows paths to WSL format
#'   (default TRUE).
#' @param os Character; operating system, auto-detected via [detect_os()] by
#'   default.
#'
#' @return Normalized path as a character string.
#'
#' @examples
#' norm_path("C:\\\\Users\\\\haohe\\\\GitHub\\\\myScript")
#'
#' @export
norm_path <- function(path, standardize = TRUE, os = detect_os()) {
  path <- gsub("/+", "/", gsub("\\\\", "/", trimws(path)))
  if (standardize && grepl("^[A-Za-z]:/", path) && os == "wsl") {
    path <- file.path("/mnt", tolower(substr(path, 1, 1)),
                      sub("^[A-Za-z]:/", "", path))
  }
  path
}


#' Standard Error of the Mean
#'
#' @param x Numeric vector.
#' @param na.rm Logical; remove missing values before computing the standard
#'   deviation (default `TRUE`).
#' @param ... Additional arguments passed to [stats::sd()].
#'
#' @return Standard error of the mean (`sd / sqrt(n)`).
#'
#' @examples
#' se(c(1, 2, 3, 4))
#' se(c(1, 2, NA, 4))
se <- function(x, na.rm = TRUE, ...) {
  if (na.rm) x <- x[!is.na(x)]
  sd(x, ...) / sqrt(length(x))
}


#' Tidy Well Style
#'
#' Removes leading zeros from the numeric part of well positions,
#' e.g. `"A01"` -> `"A1"`.
#'
#' @param wells Character vector of well positions (e.g. `"A01"`, `"B12"`).
#'
#' @return Well positions with leading zeros removed.
#'
#' @examples
#' transform_well_style(c("A1", "B01", "C003", "D5"))
#'
#' @export
transform_well_style <- function(wells) {sub("^([A-Za-z]+)0+", "\\1", wells)}


#' Left-Join Two Tables on a Key, Renaming Conflicting Columns
#'
#' Efficiently joins `y` onto `x` by key column(s) using data.table, renaming
#' any non-key columns shared by both tables with a `.y` suffix to avoid
#' name collisions.
#'
#' @param x Left table (data.frame/tibble); all its rows are kept.
#' @param y Right table whose columns are matched onto `x`; never modified
#'   (copied internally before any renaming).
#' @param by Character vector of join key columns (default `"uid"`).
#'
#' @return A data.table with one row per row of `x` (rows of `y` matched on
#'   the key). Non-key columns shared by both tables are present once, from
#'   `x`, with the `y` copy renamed to a `.y` suffix.
list_cbind_efficient <- function(x, y, by = c("uid")) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)
  common <- setdiff(intersect(names(x), names(y)), by)
  if (length(common)) {
    y <- data.table::copy(y)
    data.table::setnames(y, common, paste0(common, ".y"))
  }
  y[x, on = by]
}



#' Remove Outlier Rows by Group
#'
#' Flags each row as an outlier per group using the IQR rule
#' ([rstatix::is_outlier()]) across a set of (key) columns and removes any row
#' flagged by at least one checked column.
#'
#' @param data Data frame to filter.
#' @param grouping_vars Character vector of columns defining the per-group
#'   context in which outliers are assessed.
#' @param keys Character vector of (regex) fragments matched via
#'   [stringr::str_subset()] against the column names of `data`; all matching
#'   columns are used as outlier features.
#' @param coef Numeric multiplier of the IQR used by the outlier rule.
#'
#' @return `data` with rows flagged as outliers by any checked column removed.
#'   Rows whose checked columns are all missing are kept (no outlier flag).
remove_outlier <- function(
    data,
    grouping_vars,
    keys = c("EquivalentDiameter", "Eccentricity",
             "Intensity_MeanIntensity_", "Intensity_IntegratedIntensity_"),
    coef = 1.5
) {

  col_keys <- keys %>%
    map(\(x) str_subset(colnames(data), x)) %>%
    list_c()

  data_outlier <- data %>%
    mutate(across(any_of(col_keys),
                  \(x) rstatix::is_outlier(x, coef)),
           .by = any_of(grouping_vars)) %>%
    dplyr::select(any_of(col_keys))

  kept_rows <- rowSums(data_outlier, na.rm = TRUE) < 1

  data[kept_rows, ]
}



#' Calculate DR Curve Derivatives from a GAM Fit
#'
#' Fits `mgcv::gam(fluorescence ~ s(temperature))` to the data and returns the
#' first derivative over a 100-point temperature grid, as computed by
#' [gratia::derivatives()].
#'
#' @param data Data frame with `fluorescence` and `temperature` columns.
#' @param eps Numeric; finite-difference epsilon passed to
#'   [gratia::derivatives()].
#' @param unconditional Logical; pass to [gratia::derivatives()] to account for
#'   smoothing-parameter uncertainty.
#'
#' @return A tibble with columns `temperature`, `derivative` and
#'   `derivative_se`.
calculate_deri <- function(data, eps = 1e-07, unconditional = FALSE) {
  mod <- mgcv::gam(fluorescence ~ s(temperature), data = data)
  res <- gratia::derivatives(mod, n = 100, order = 1,
                             eps = eps, unconditional = unconditional)
  tibble(temperature = res$temperature,
         derivative = res$.derivative,
         derivative_se = res$.se)
}



# **********************************************************************
# read sqlite ------------------------------------------------------
# **********************************************************************

#' Read All Tables from a SQLite Database
#'
#' Reads every table in a SQLite file into a named list of tibbles, dropping
#' constant columns. Internal/technical tables can be skipped.
#'
#' @param f Path to the SQLite database file.
#' @param skip_tables Character vector of table names to skip.
#'
#' @return Named list of tibbles (one per table); `NULL` on error.
read_sqlite <- function(
    f, skip_tables = c("metadata", "reduction_pca_variance", "sqlite_sequence")
) {
  print(paste0("read table: ", f))
  conn <- DBI::dbConnect(RSQLite::SQLite(), dbname = f)
  on.exit(DBI::dbDisconnect(conn), add = TRUE)

  tables <- setdiff(DBI::dbListTables(conn), skip_tables)
  tryCatch(
    setNames(
      lapply(tables, \(x) dplyr::tbl(conn, x) %>%
        dplyr::collect() %>%
        janitor::remove_constant(quiet = FALSE)),
      tables),
    error = function(e) { warning(e); NULL })
}


# **********************************************************************
# read excel metadata ------------------------------------------------------
# **********************************************************************

#' Read Plate Metadata from Excel File
#'
#' Reads 384-well plate metadata from Excel files with a specific layout
#' (data in the `B2:Y17` range of each sheet). Supports wide and long output.
#'
#' @param file Path to the Excel file.
#' @param format Output format: "wide" (one row per well, one column per
#'   sheet) or "long" (one row per well-sheet combination).
#' @param add_directory Logical; adds a `directory` column parsed from the
#'   file path: the path segment matching `.*-Measurement [0-9]` with any
#'   `__...` suffix removed.
#' @param colname_prefix String prefix prepended to all column names in the
#'   wide format (including the `well` column).
#'
#' @return A tibble of plate metadata; `NULL` if `file` does not exist.
#'
#' @details Metadata is read from the `B2:Y17` range of every sheet (a 384-well
#'   plate: 16 rows x 24 columns). All-`NA` wells are dropped in wide format.
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

  sheets <- readxl::excel_sheets(file)
  print(str_glue("expected sheets: {str_c(sheets, collapse=', ')}"))

  # each sheet -> flat vector of values (B2:Y17)
  read_sheet <- \(s) readxl::read_xlsx(
    file, sheet = s, range = "B2:Y17",
    col_names = FALSE, col_types = "text", .name_repair = "minimal"
  ) %>% as.matrix() %>% as.vector()

  wells <- str_c(rep(LETTERS[1:16], 24), rep(1:24, each = 16))

  if (format == "wide") {
    out <- bind_cols(
      well = wells,
      map_dfc(sheets, \(s) enframe(read_sheet(s), name = NULL, value = s))
    ) %>% filter(!if_all(!well, is.na))

    if (add_directory) {
      out <- out %>% mutate(
        directory = file %>%
          str_split("[/\\\\]", simplify = TRUE) %>%
          str_subset(".*-Measurement [0-9]") %>%
          str_split("(__)", simplify = TRUE) %>% .[1],
        .before = 1)
    }
    if (nzchar(colname_prefix)) out <- out %>% rename_with(~ str_c(colname_prefix, .x))

  } else { # long
    out <- bind_rows(
      map(sheets, \(s) bind_cols(
        Metadata_well = wells,
        enframe(read_sheet(s), name = NULL, value = "value")
      ) %>% filter(!is.na(value))),
      .id = "Metadata_sheet")
  }
  out
}


#' Read PlateReader Excel Files
#'
#' Reads plate-reader (e.g. PerkinElmer) Excel files: locates the "A" row of
#' each sheet, extracts the well grid, and returns long-format data with a
#' `Metadata_well` column.
#'
#' @param f Path to the Excel file.
#' @param well_type Plate format: "24", "48", "96" or "384".
#' @param sheet Sheet name(s) to read; `NULL` reads all sheets.
#' @param invert Logical; flip the plate rows and columns before reading, to
#'   correct for plates that were placed/read in the wrong orientation.
#'
#' @return A tibble with columns `Metadata_sheet`, `Metadata_well`, `value`.
#'
#' @export
read_plateReader <- function(
    f, well_type = "96", sheet = NULL, invert = FALSE
) {
  well_type <- match.arg(as.character(well_type), c("24", "48", "96", "384"))
  n_row <- c("24" = 4, "48" = 6, "96" = 8, "384" = 16)[[well_type]]
  n_col <- as.integer(well_type) / n_row

  sheets <- sheet %||% readxl::excel_sheets(f)
  res <- vector("list", length(sheets)) %>% set_names(sheets)

  for (s in sheets) {
    print(str_glue("read sheet: {s}"))

    # find the row containing "A"
    start_loc <- tryCatch(
      readxl::read_xlsx(f, s, range = "A1:A100", col_names = FALSE,
                        .name_repair = "minimal") %>%
        tibble::deframe() %>% stringr::str_equal("A") %>% which(),
      error = function(e) NULL)
    if (is.null(start_loc)) {
      print("No 'A' found in first column, check file. skip sheet!")
      next
    }

    data_range <- str_glue("B{start_loc}:{LETTERS[n_col + 1]}{start_loc + n_row - 1}")
    raw <- readxl::read_xlsx(f, s, range = data_range, col_names = FALSE,
                             .name_repair = "minimal")
    if (isTRUE(invert)) raw <- raw[nrow(raw):1, ncol(raw):1]

    res[[s]] <- raw %>%
      as.matrix() %>% as.vector() %>% as.numeric() %>%
      enframe(name = NULL) %>%
      bind_cols(Metadata_well = str_c(
        rep(LETTERS[seq_len(n_row)], n_col),
        rep(seq_len(n_col), each = n_row))) %>%
      filter(!is.na(value))
  }
  list_rbind(res, names_to = "Metadata_sheet")
}


# **********************************************************************
# dose response ------------------------------------------------------
# **********************************************************************

#' Fit a 4-Parameter Logistic Dose-Response Model
#'
#' Wrapper around [dr4pl::dr4pl()] that returns `NULL` (instead of throwing)
#' when the fit fails.
#'
#' @param data Data frame with `dose` and `response` columns.
#' @param ... Additional arguments passed to [dr4pl::dr4pl()].
#'
#' @return A `dr4pl` model object, or `NULL` on failure.
model_func <- function(data, ...) {
  res <- try(dr4pl::dr4pl(response ~ dose, data = data, ...), silent = TRUE)
  if (inherits(res, "try-error")) NULL else res
}


#' Tidy dr4pl Coefficients
#'
#' Extracts estimates and half-width confidence intervals from a `dr4pl`
#' model, converts log10(IC50) to IC50, and returns one tidy row.
#'
#' @param model A `dr4pl` model object.
#'
#' @return A one-row tibble with columns such as `Upper_Estimate`,
#'   `IC50_Estimate`, `IC50_CI`, `Slope_Estimate`, `Lower_Estimate`
#'   (formatted as strings); `NULL` if `model` is not a `dr4pl` object.
#'
#' @details Expects the summary coefficient table of `dr4pl` (v2) with rows
#'   `UpperLimit, Log10(EC50), Slope, LowerLimit` and columns
#'   `Estimate, StdErr, 2.5 %, 97.5 %`. `Log10(EC50)` is in log10 units: the
#'   estimate is back-transformed to IC50 and its half-width CI rebuilt from
#'   the (log10-scale) confidence limits, while the StdErr is left untouched.
coef_tidy <- function(model) {
  if (!inherits(model, "dr4pl")) return(NULL)

  ci <- as.data.frame(summary(model)[["coefficients"]])  # theta_1..4 x Estimate/StdErr/Limits
  ci$CI <- (ci[[4]] - ci[[3]]) / 2                       # half-width of CI

  ic50 <- ci[2, , drop = FALSE]
  ic50$Estimate <- 10^ic50$Estimate                     # log10 IC50 -> IC50
  ic50$CI <- (10^ic50[[4]] - 10^ic50[[3]]) / 2          # CI rebuilt in IC50 units
  rownames(ic50) <- "IC50"
  out <- rbind(ci[-2, ], ic50)[c(1, 4, 2, 3), c(1, 5)]   # theta1, IC50, theta3, theta4

  out[] <- lapply(out, \(x) trimws(format(x, digits = 3, scientific = TRUE)))
  out %>%
    as_tibble(rownames = "type") %>%
    mutate(type = str_remove_all(type, "Limit")) %>%
    pivot_wider(names_from = type, names_glue = "{type}_{.value}",
                values_from = c(Estimate, CI))
}


#' Predict a dr4pl Curve (with Optional Bootstrap Confidence Band)
#'
#' Evaluates a fitted `dr4pl` model on a log-spaced dose grid; optionally
#' resamples residuals to obtain a percentile confidence band.
#'
#' @param model A `dr4pl` model object.
#' @param se Logical; whether to compute the bootstrap confidence band.
#' @param normalize_residual Logical; Yeo-Johnson-transform residuals before
#'   resampling (when `se = TRUE`).
#' @param level Confidence level of the band (when `se = TRUE`).
#' @param nboot Number of bootstrap replicates (when `se = TRUE`).
#'
#' @return A data.frame with `x` (dose) and `y` (predicted response); with
#'   `se = TRUE` also `ymin` and `ymax`.
pred_func <- function(model, se = FALSE, normalize_residual = TRUE,
                      level = 0.95, nboot = 200) {
  if (!inherits(model, "dr4pl")) return(NULL)

  dose <- sort(unique(model$data$Dose))
  single_dose <- length(dose) == 1
  from <- if (dose[1] <= 0 && !single_dose) 0.8 * dose[2] else 0.8 * dose[1]
  to   <- 1.2 * last(dose)
  # degenerate fits (single dose, or non-positive doses) -> safe positive grid
  if (!is.finite(from) || from <= 0) from <- 0.01
  if (!is.finite(to) || to <= 0) to <- 1
  xseq <- exp(seq(log(from), log(to), length.out = 200))
  if (dose[1] <= 0 && !single_dose) {
    xseq <- c(xseq, exp(seq(log(dose[2] / 1000), log(dose[2]), length.out = 50)))
  }

  pred <- dr4pl::MeanResponse(model$parameters, xseq)
  if (!se) return(data.frame(x = xseq, y = pred))

  # bootstrap residuals
  pred0 <- dr4pl::MeanResponse(model$parameters, model$data$Dose)
  res <- pred0 - model$data$Response
  if (normalize_residual) {
    res <- bestNormalize::yeojohnson(res, standardize = FALSE)$x.t %>%
      scale(center = TRUE, scale = FALSE) %>% as.vector()
  }

  bootres <- suppressMessages(future_map_dfc(seq_len(nboot), \(i) {
    mboot <- try(dr4pl::dr4pl(model$data$Dose,
                              pred0 + sample(res, length(pred0), replace = TRUE)),
                 silent = TRUE)
    if (inherits(mboot, "try-error")) rep(NA_real_, length(xseq))
    else dr4pl::MeanResponse(mboot$parameters, xseq)
  }, .options = furrr_options(seed = NULL))) %>% as.matrix()

  data.frame(
    x = xseq, y = pred,
    ymin = apply(bootres, 1, quantile, (1 - level) / 2, na.rm = TRUE),
    ymax = apply(bootres, 1, quantile, (1 + level) / 2, na.rm = TRUE))
}


#' Fit Dose-Response Curves per Group
#'
#' Fits a 4-parameter logistic model to each group of `df` in parallel,
#' returning models plus tidy coefficients and predictions.
#'
#' @param df Data frame with `dose`, `response`, and optional grouping columns.
#' @param dose Character; name of the dose column.
#' @param response Character; name of the response column.
#' @param group Character vector of grouping columns (one model per combination).
#' @param n_worker Number of parallel workers (`> 1` enables multisession).
#' @param ... Additional arguments passed to [dr4pl::dr4pl()] (via [model_func()]).
#'
#' @return A tibble with columns `data`, `model`, `convergence`,
#'   `method_robust`, `coef`, `pred` (plus grouping columns).
fit_dose_response <- function(df, dose, response, group = NULL,
                              n_worker = 8, ...) {
  old_plan <- future::plan()
  if (n_worker > 1) future::plan(future::multisession, workers = n_worker)
  on.exit(future::plan(old_plan), add = TRUE)

  df %>%
    select(any_of(c(dose, response, group))) %>%
    rename(dose = all_of(dose), response = all_of(response)) %>%
    mutate(across(c(dose, response), \(x) as.numeric(as.character(x)))) %>%
    group_by(across(any_of(group))) %>%
    nest() %>%
    mutate(model = future_map(data, model_func, ...)) %>%
    mutate(
      convergence   = map_lgl(model, \(m) !is.null(m) && m[["convergence"]]),
      method_robust = map_chr(model, \(m) if (is.null(m)) "None" else m[["method.robust"]]),
      coef          = map(model, coef_tidy),
      pred          = map(model, pred_func)
    ) %>%
    ungroup()
}


#' Prepare Dose-Response Data for Plotting
#'
#' Computes per-dose mean/SE and top inhibition, and builds plot labels with
#' IC50 and its confidence interval.
#'
#' @param fitted Output of [fit_dose_response()].
#' @param group Character vector of grouping columns used for the labels.
#'
#' @return A tibble with `stats`, `top_inhit`, `uid` and `group_label` columns.
plot_dose_response_dataPrepare <- function(fitted, group = NULL) {
  fitted %>%
    mutate(
      stats = map(data, \(d) d %>%
                    group_by(dose) %>%
                    summarise(mean = mean(response),
                              se = sd(response) / sqrt(n()), .groups = "drop")),
      top_inhit = map_dbl(pred, \(p) round(100 * (1 - min(p$y) / max(p$y)), 0))
    ) %>%
    unnest(coef, keep_empty = TRUE) %>%
    { if (is.null(group)) mutate(., uid = 1)
      else unite(., "uid", all_of(group), sep = "\n", remove = FALSE) } %>%
    mutate(group_label = str_glue(
      "{uid}\n IC50: {str_replace(IC50_Estimate, '[+]', '')} ± ",
      "{str_replace(IC50_CI, '[+]', '')}\n Top Inhibition: {top_inhit}%"))
}


# **********************************************************************
# calc  distance ------------------------------------------------------
# **********************************************************************

#' Parallel Population Distance Calculation
#'
#' Computes pairwise (or reference-vs-all) distances between groups of rows,
#' with optional permutation-based p-values. Supported distances: MMD
#' (Gaussian RBF), energy distance, Wasserstein (p = 2) and Mahalanobis.
#'
#' @param mat Matrix or data.frame of features (rows = samples/cells).
#' @param group Vector or factor assigning each row to a group.
#' @param method One of "mmd", "edistance", "wasserstein", "mahalanobis".
#' @param ref Optional reference group name(s); distance is then computed only
#'   from each reference to every group and returned as a tidy tibble. If NULL,
#'   all pairwise distances are returned (as a `stats::dist` object).
#' @param sigma Optional bandwidth for the Gaussian RBF kernel in MMD; if NULL,
#'   the median pairwise distance is used.
#' @param max_cells Maximum observations per group (groups are subsampled).
#' @param workers Number of CPU cores to use.
#' @param n_perm Number of permutations for a permutation-based p-value
#'   (0 = skip, default).
#' @param seed Seed for the permutation null (deterministic per group pair).
#' @param alternative Alternative hypothesis of the permutation test:
#'   "two.sided" or "greater". (Mahalanobis uses a chi-square p-value instead.)
#'
#' @return A tibble (when `ref` is given) or a `stats::dist` object (a list of
#'   two `dist` objects — `distance` and `p_value` — when `n_perm > 0`).
calc_distance <- function(
    mat, group,
    method = c("mmd", "edistance", "wasserstein", "mahalanobis"),
    ref = NULL, sigma = NULL, max_cells = 2000,
    workers = floor(future::availableCores() / 2),
    n_perm = 0, seed = 123,
    alternative = c("two.sided", "greater")) {

  method <- match.arg(method)
  alternative <- match.arg(alternative)
  n_perm <- max(0L, as.integer(n_perm))
  mat <- as.matrix(mat)
  group_vec <- as.character(group)
  all_groups <- unique(group_vec)
  if (!is.null(ref)) {
    ref <- as.character(ref)
    if (!all(ref %in% all_groups)) {
      stop("All elements in 'ref' must exist in unique values of 'group'.")
    }
  }

  splits <- split(as.data.frame(mat), group_vec)

  # ---- helpers
  mmd2_from_K <- function(K, nX, nY) {
    sqrt(max(0,
      (sum(K[1:nX, 1:nX]) - nX) / (nX * (nX - 1)) +
      (sum(K[(nX + 1):(nX + nY), (nX + 1):(nX + nY)]) - nY) / (nY * (nY - 1)) -
      2 * mean(K[1:nX, (nX + 1):(nX + nY)])))
  }
  compute_mmd <- function(X, Y, bw) {
    D2 <- as.matrix(dist(rbind(X, Y)))^2
    if (is.null(bw)) bw <- sqrt(0.5 * median(D2[D2 > 0]))
    mmd2_from_K(exp(-D2 / (2 * bw^2)), nrow(X), nrow(Y))
  }
  perm_p <- function(null, obs) {
    if (obs == 0 && all(null == 0)) return(1)
    pg <- (sum(null >= obs) + 1) / (length(null) + 1)
    if (alternative == "greater") pg
    else 2 * min(pg, (sum(null <= obs) + 1) / (length(null) + 1))
  }
  subsample <- function(g, m) {
    if (nrow(g) > m) g[sample.int(nrow(g), m), , drop = FALSE] else g
  }
  pinv <- function(M) {                  # Moore-Penrose (robust to singular cov)
    s <- svd(M)
    w <- ifelse(s$d > max(dim(M)) * max(s$d) * .Machine$double.eps, 1 / s$d, 0)
    s$v %*% (t(s$u) * w)
  }

  if (method == "mahalanobis") {
    means <- t(sapply(splits, colMeans))
    inv_cov <- pinv(stats::cov(mat))
  }

  # ---- one pair
  calc_single_dist <- function(g1n, g2n) {
    if (g1n == g2n) {
      return(list(distance = 0, statistic = 0,
                  p_value = if (n_perm > 0) NA_real_ else NULL))
    }
    g1 <- subsample(as.matrix(splits[[g1n]]), max_cells)
    g2 <- subsample(as.matrix(splits[[g2n]]), max_cells)

    if (method == "mahalanobis") {
      diff <- means[g1n, ] - means[g2n, ]
      d <- sqrt(as.numeric(t(diff) %*% inv_cov %*% diff))
      # squared Mahalanobis distance ~ chi-square(df = #features)
      p <- if (n_perm > 0) pchisq(d^2, ncol(mat), lower.tail = FALSE) else NULL
      return(list(distance = d, statistic = d, p_value = p))
    }
    if (n_perm <= 0) {
      d <- switch(method,
        "mmd"         = compute_mmd(g1, g2, sigma),
        "edistance"   = as.matrix(energy::edist(rbind(g1, g2), c(nrow(g1), nrow(g2))))[1, 2],
        "wasserstein" = transport::wasserstein(transport::pp(g1), transport::pp(g2), p = 2))
      return(list(distance = d, statistic = d, p_value = NULL))
    }

    # ---- permutation
    set.seed(seed + as.integer(floor(abs(sum(utf8ToInt(paste0(g1n, g2n))) %% 1e6))))
    nX <- nrow(g1); nY <- nrow(g2)

    if (method == "mmd") {
      D2 <- as.matrix(dist(rbind(g1, g2)))^2
      bw <- if (is.null(sigma)) sqrt(0.5 * median(D2[D2 > 0])) else sigma
      K <- exp(-D2 / (2 * bw^2))
      obs <- mmd2_from_K(K, nX, nY)
      null <- replicate(n_perm, {
        idx <- sample.int(nX + nY)
        mmd2_from_K(K[idx, idx], nX, nY)
      })
    } else {
      m <- if (method == "edistance") 400 else 1000
      ga <- subsample(g1, m); gb <- subsample(g2, m)
      m <- min(nrow(ga), nrow(gb), m)         # same pool size for obs and null
      ga <- ga[seq_len(m), , drop = FALSE]
      gb <- gb[seq_len(m), , drop = FALSE]
      obs <- if (method == "edistance") {
        as.matrix(energy::edist(rbind(ga, gb), c(m, m)))[1, 2]
      } else {  # wasserstein needs equal sample sizes
        transport::wasserstein(transport::pp(ga), transport::pp(gb), p = 2)
      }
      null <- replicate(n_perm, {
        idx <- sample.int(2 * m)
        a <- rbind(ga, gb)[idx[1:m], , drop = FALSE]
        b <- rbind(ga, gb)[idx[(m + 1):(2 * m)], , drop = FALSE]
        if (method == "edistance") as.matrix(energy::edist(rbind(a, b), c(m, m)))[1, 2]
        else transport::wasserstein(transport::pp(a), transport::pp(b), p = 2)
      })
    }
    list(distance = obs, statistic = obs, p_value = perm_p(null, obs))
  }

  # ---- os-agnostic parallel plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  strategy <- if (.Platform$OS.type == "unix") future::multicore else future::multisession
  future::plan(strategy, workers = workers)
  opts <- furrr::furrr_options(seed = TRUE)

  # ---- branch 1: reference vs all (tidy tibble)
  if (!is.null(ref)) {
    pairs <- expand.grid(group_1 = ref, group_2 = all_groups, stringsAsFactors = FALSE)
    res <- furrr::future_map2(pairs$group_1, pairs$group_2,
                              ~ calc_single_dist(.x, .y), .options = opts)
    pairs$distance <- vapply(res, \(r) r$distance, numeric(1))
    if (n_perm > 0) {
      pairs$statistic <- vapply(res, \(r) r$statistic, numeric(1))
      pairs$p_value   <- vapply(res, \(r) r$p_value, numeric(1))
    }
    return(tibble::as_tibble(pairs))
  }

  # ---- branch 2: all pairwise (dist objects)
  n_groups <- length(all_groups)
  idx <- which(lower.tri(matrix(0, n_groups, n_groups)), arr.ind = TRUE)
  g1v <- all_groups[idx[, 1]]; g2v <- all_groups[idx[, 2]]
  res <- furrr::future_map2(g1v, g2v, ~ calc_single_dist(.x, .y), .options = opts)

  d <- vapply(res, \(r) r$distance, numeric(1))
  dist_mat <- matrix(0, n_groups, n_groups, dimnames = list(all_groups, all_groups))
  dist_mat[cbind(idx[, 1], idx[, 2])] <- d
  dist_mat[cbind(idx[, 2], idx[, 1])] <- d

  if (n_perm > 0) {
    p <- vapply(res, \(r) r$p_value, numeric(1))
    p_mat <- matrix(NA_real_, n_groups, n_groups, dimnames = list(all_groups, all_groups))
    p_mat[cbind(idx[, 1], idx[, 2])] <- p
    p_mat[cbind(idx[, 2], idx[, 1])] <- p
    return(list(distance = as.dist(dist_mat), p_value = as.dist(p_mat)))
  }
  as.dist(dist_mat)
}


# **********************************************************************
# df to heatmap ------------------------------------------------------
# **********************************************************************

#' Heatmap from a Data Frame (Self-Contained)
#'
#' Turns pre-aggregated data into a named list of ggplot heatmaps (one per
#' metric), with optional z-score scaling, hierarchical clustering, and split
#' panels.
#'
#' @param df Pre-aggregated data frame (at most one row per
#'   `(row_by, col_by, split_by)` cell).
#' @param metrics Character vector of value columns to plot (one heatmap each).
#' @param row_by Character; column plotted on the y-axis.
#' @param col_by Character; column plotted on the x-axis.
#' @param color Named numeric vector of 2 or 3 colors (names = colors, values =
#'   color positions); NULL for an automatic scale.
#' @param border Color of the tile borders.
#' @param split_by Optional character vector of columns used to split the
#'   heatmap into panels.
#' @param scale Optional "row"/"col" to z-score each row/column.
#' @param scale_scope Whether z-scoring is computed per split panel ("split")
#'   or across all data ("pooled").
#' @param row_cluster, col_cluster Logical; hierarchically cluster the axis
#'   objects (objects without finite values keep their original order).
#' @param cluster_method Method for [stats::hclust()].
#' @param cluster_distance Distance method for [stats::dist()].
#' @param legend_position Legend position: "left", "right", "top" or "bottom".
#' @param x_text_angle Rotation of the x-axis labels.
#' @param na_color Color for missing values.
#'
#' @return A named list of ggplot objects (one per metric).
df2heatmap <- function(df, metrics, row_by, col_by, color = NULL,
                       border = "grey90", split_by = NULL, scale = NULL,
                       scale_scope = c("split", "pooled"),
                       row_cluster = FALSE, col_cluster = FALSE,
                       cluster_method = "complete",
                       cluster_distance = "euclidean",
                       legend_position = "right", x_text_angle = 0,
                       na_color = "grey90") {
  
  scale_scope <- match.arg(scale_scope)
  legend_position <- match.arg(tolower(legend_position),
                               c("left", "right", "top", "bottom"))
  if (!is.null(scale)) scale <- match.arg(tolower(scale), c("row", "col"))
  cluster_method <- match.arg(cluster_method,
                              c("complete", "single", "average", "ward.D", "ward.D2",
                                "mcquitty", "median", "centroid"))
  if (!is.logical(row_cluster) || length(row_cluster) != 1 || is.na(row_cluster))
    stop("`row_cluster` must be TRUE or FALSE.")
  if (!is.logical(col_cluster) || length(col_cluster) != 1 || is.na(col_cluster))
    stop("`col_cluster` must be TRUE or FALSE.")
  if (!is.character(cluster_distance) || length(cluster_distance) != 1 || is.na(cluster_distance))
    stop("`cluster_distance` must be a single character string.")
  
  key_cols <- unique(c(row_by, col_by, split_by))
  missing_cols <- setdiff(unique(c(key_cols, metrics)), names(df))
  if (length(missing_cols))
    stop("`df` is missing columns: ", paste(missing_cols, collapse = ", "))
  if (any(c(row_by, col_by) %in% split_by))
    stop("`row_by` / `col_by` must not also appear in `split_by`")
  
  # df2heatmap expects pre-aggregated data
  n_dup <- df %>%
    dplyr::count(dplyr::across(dplyr::all_of(key_cols))) %>%
    dplyr::filter(n > 1) %>% nrow()
  if (n_dup > 0) {
    stop(stringr::str_glue(
      "Found {n_dup} ({paste(key_cols, collapse = ' / ')}) cell(s) with >1 df row. ",
      "df2heatmap expects pre-aggregated data; aggregate first, e.g. ",
      "df %>% group_by({paste(key_cols, collapse = ', ')}) %>% ",
      "summarise(across(all_of(metrics), mean, na.rm = TRUE), .groups = 'drop')"))
  }
  
  # preserve input order
  df <- df %>% dplyr::mutate(dplyr::across(
    dplyr::all_of(key_cols),
    \(x) if (is.factor(x)) x else factor(x, levels = unique(x))))
  
  split_combos <- if (is.null(split_by)) list(tibble::tibble())
  else df %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(split_by))) %>%
    split(seq_len(nrow(.)))
  
  # ---- z-score along an axis (optionally per split panel)
  zscore_axis <- function(long, axis, by_split) {
    grp <- if (by_split) unique(c(split_by, axis)) else axis
    long %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
      dplyr::mutate(display = as.numeric(base::scale(value))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(display = replace(display, is.nan(display), 0))  # constant axis -> 0
  }
  
  # ---- hierarchical clustering of one axis (shared across split panels)
  cluster_order <- function(long, axis) {
    other <- setdiff(c(row_by, col_by), axis)
    mat <- long %>%
      dplyr::select(dplyr::all_of(c(axis, other)), display) %>%
      tidyr::pivot_wider(names_from = dplyr::all_of(other), values_from = display) %>%
      tibble::column_to_rownames(axis) %>%
      as.matrix()
    if (axis == col_by) mat <- t(mat)  # work on rows uniformly
    
    # cluster objects with any finite value; append the rest in original order
    ok <- rowSums(is.finite(mat)) > 0
    keep <- rownames(mat)[ok]; dropped <- rownames(mat)[!ok]
    if (length(keep) > 1) {
      m <- mat[keep, , drop = FALSE]
      # replace NA with the feature mean (0 if the feature is all NA)
      m[] <- apply(m, 2, \(col) {
        if (all(!is.finite(col))) 0 else replace(col, !is.finite(col), mean(col[is.finite(col)]))
      })
      ord <- rownames(m)[stats::hclust(stats::dist(m, method = cluster_distance),
                                       method = cluster_method)$order]
      c(ord, dropped)
    } else {
      c(keep, dropped)
    }
  }
  
  # ---- one metric -> one (list of) heatmap
  make_one_metric <- function(metric) {
    long <- df %>%
      dplyr::select(dplyr::all_of(c(key_cols, metric))) %>%
      dplyr::rename(value = dplyr::all_of(metric))
    
    long <- if (is.null(scale)) {
      dplyr::mutate(long, display = value)
    } else {
      zscore_axis(long, if (scale == "row") row_by else col_by,
                  by_split = scale_scope == "split")
    }
    
    row_order <- if (row_cluster) cluster_order(long, row_by) else NULL
    col_order <- if (col_cluster) cluster_order(long, col_by) else NULL
    if (!is.null(row_order)) long[[row_by]] <- factor(long[[row_by]], levels = row_order)
    if (!is.null(col_order)) long[[col_by]] <- factor(long[[col_by]], levels = col_order)
    
    # shared color scale (1%/99% quantiles)
    dvals <- long$display
    q02 <- if (all(is.na(dvals))) 0 else unname(stats::quantile(dvals, 0.01, na.rm = TRUE))
    q98 <- if (all(is.na(dvals))) 1 else unname(stats::quantile(dvals, 0.99, na.rm = TRUE))
    
    color_spec <- color
    if (is.null(color_spec)) {
      color_spec <- if (is.null(scale) && q02 >= 0)
        c("white" = q02, "#762a83" = q98)
      else
        c("#2166ac" = q02, "white" = 0, "#b2182b" = q98)
    } else {
      if (!is.numeric(color_spec)) stop("`color` must be a named numeric vector.")
      if (is.null(names(color_spec)) || anyNA(names(color_spec)) || any(!nzchar(names(color_spec))))
        stop("`color` must be a named numeric vector, e.g. c('white' = 0, '#762a83' = 1)")
      if (!length(color_spec) %in% c(2, 3)) stop("`color` must contain exactly 2 or 3 colors.")
      if (!all(is.finite(color_spec))) stop("Values in `color` must be finite.")
      if (is.unsorted(unname(color_spec), strictly = TRUE))
        stop("Values in `color` must be strictly increasing.")
    }
    
    fill_name <- if (is.null(scale)) "value" else paste(scale, "z-score")
    fill_scale <- if (length(color_spec) == 2) {
      ggplot2::scale_fill_gradient(low = names(color_spec)[1], high = names(color_spec)[2],
                                   limits = unname(color_spec), oob = scales::squish,
                                   na.value = na_color, name = fill_name)
    } else {
      ggplot2::scale_fill_gradient2(low = names(color_spec)[1], mid = names(color_spec)[2],
                                    high = names(color_spec)[3], midpoint = color_spec[[2]],
                                    limits = c(color_spec[[1]], color_spec[[3]]),
                                    oob = scales::squish, na.value = na_color, name = fill_name)
    }
    
    # split panels
    plots <- lapply(split_combos, function(combo) {
      if (is.null(split_by)) {
        sub <- long; ttl <- ""
      } else {
        sub <- dplyr::semi_join(long, combo, by = split_by)
        ttl <- paste(split_by, "=", as.character(unlist(combo)), collapse = ", ")
      }
      ggplot2::ggplot(sub, ggplot2::aes(x = .data[[col_by]], y = .data[[row_by]], fill = display)) +
        ggplot2::geom_tile(color = border, linewidth = 0.2) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::scale_y_discrete(drop = FALSE) +
        fill_scale +
        theme_hh(axis_text_x_angle = x_text_angle) +
        ggplot2::labs(x = col_by, y = row_by, title = ttl)
    })
    
    dir <- if (legend_position %in% c("top", "bottom")) "horizontal" else "vertical"
    (patchwork::wrap_plots(plots, nrow = 1, guides = "collect") &
        ggplot2::theme(legend.position = legend_position, legend.direction = dir)) +
      patchwork::plot_annotation(title = metric)
  }
  
  setNames(lapply(metrics, make_one_metric), metrics)
}


# **********************************************************************
# ggplot2 utils ------------------------------------------------------
# **********************************************************************


#' Add an Automatic Log Scale to a ggplot
#'
#' A lazy scale: when added to a plot with `+`, it detects the variable mapped
#' to the selected axis from the plot's aesthetics, computes positive data
#' limits, and generates logarithmic breaks and labels.
#'
#' @param axis Character. Axis to transform, "x" or "y".
#' @param limits Numeric vector of length two; optional scale limits
#'   (auto-computed from the mapped variable if NULL).
#' @param breaks Numeric vector; positions of log-scale breaks
#'   (decade breaks generated automatically if NULL).
#' @param base Numeric. Logarithm base (default 10).
#' @param labels Function or character vector; tick labels
#'   (base^exponent if NULL).
#' @param ... Additional arguments passed to [ggplot2::scale_x_continuous()] /
#'   [ggplot2::scale_y_continuous()].
#'
#' @return An object of class `scale_log_auto`, interpreted by
#'   [ggplot_add.scale_log_auto()] when added to a plot.
#'
#' @examples
#' \dontrun{
#' ggplot(df, aes(conc, response)) + geom_point() + scale_log("x")
#' }
#'
#' @export
scale_log <- function(axis = c("x", "y"), limits = NULL, breaks = NULL,
                      base = 10, labels = NULL, ...) {
  structure(list(axis = match.arg(axis), limits = limits, breaks = breaks,
                 base = base, labels = labels, params = list(...)),
            class = "scale_log_auto")
}


#' @method ggplot_add scale_log_auto
#' @export
ggplot_add.scale_log_auto <- function(object, plot, object_name) {
  axis <- object$axis
  aes_var <- if (axis == "x") plot$mapping$x else plot$mapping$y
  if (is.null(aes_var)) stop("Cannot detect ", axis, " aesthetic from aes()")
  var <- rlang::as_name(aes_var)
  
  data <- plot$data
  if (is.null(data) || nrow(data) == 0) data <- plot$layers[[1]]$data
  if (is.null(data)) stop("Cannot find plot data")
  
  if (is.null(object$limits)) {
    x <- data[[var]]
    x <- x[x > 0 & is.finite(x)]
    if (!length(x)) {
      stop("Cannot determine positive finite range for ", axis, " variable: ", var)
    }
    object$limits <- range(x, na.rm = TRUE)
  }
  if (is.null(object$breaks)) {
    rng <- object$limits
    object$breaks <- object$base^(floor(log(rng[1], object$base)):ceiling(log(rng[2], object$base)))
  }
  if (is.null(object$labels)) {
    object$labels <- \(x) as.expression(lapply(log(x, base = object$base),
                                               \(e) bquote(.(object$base)^.(e))))
  }
  
  scale_fun <- if (axis == "x") scale_x_continuous else scale_y_continuous
  scale <- do.call(scale_fun, c(list(transform = scales::log_trans(object$base),
                                     limits = object$limits, breaks = object$breaks,
                                     labels = object$labels), object$params))
  ggplot_add(scale, plot, object_name)
}

#' Minimal Publication Theme
#'
#' A `theme_classic()`-based theme with transparent backgrounds, compact
#' margins, optional x-axis label rotation and optional panel box/grid.
#'
#' `theme_hh()` also installs its color palettes as the global default ggplot2
#' scales (`scale_*_discrete()` / `scale_*_continuous()`), so every subsequent
#' plot applies them automatically. Discrete data uses `category_color_qua`
#' (interpolating through the base colors when there are more categories than
#' colors); continuous data auto-picks `continuous_color_div` when any value is
#' negative and `continuous_color_seq` otherwise. A `NULL` palette disables
#' that channel (ggplot2's default scale is used), and passing `NULL` for all
#' three disables all custom color setting.
#'
#' @param axis_text_x_angle Numeric; rotation of x-axis labels in degrees
#'   (0 = no rotation).
#' @param fontsize Base font size.
#' @param family Font family.
#' @param face Font face (e.g. "bold"), or NULL.
#' @param legend_size Legend key size.
#' @param unit Unit for `legend_size` and margins.
#' @param legend_position Legend position (e.g. "top", "right", "bottom",
#'   "left", "none").
#' @param panel_border Logical; draw a full panel border box (including the
#'   top and right sides). Default FALSE (only the left/bottom axis lines are
#'   shown).
#' @param panel_grid_major Logical; show major panel grid lines. Default FALSE.
#' @param category_color_qua Character vector of qualitative colors for
#'   categories (applies to both `fill` and `color`); NULL disables custom
#'   discrete scales.
#' @param continuous_color_seq Character vector of sequential colors used when
#'   all continuous values are >= 0; NULL falls back to
#'   `continuous_color_div`.
#' @param continuous_color_div Character vector of diverging colors used when
#'   continuous values include negatives; NULL falls back to
#'   `continuous_color_seq`.
#'
#' @return A ggplot2 theme object (with a `my_palettes` attribute).
theme_hh <- function(axis_text_x_angle = 0, fontsize = 7, family = "sans",
                     face = NULL, legend_size = 4, unit = "mm",
                     legend_position = "top",
                     panel_border = FALSE, panel_grid_major = FALSE,
                     category_color_qua = c("#e41a1c", "#4daf4a", "#377eb8",
                                            "#984ea3", "#ff7f00", "#ffff33",
                                            "#a65628", "#f781bf"),
                     continuous_color_seq = c("#fff7ec", "#fee8c8", "#fdd49e",
                                              "#fdbb84", "#fc8d59", "#ef6548",
                                              "#d7301f", "#990000"),
                     continuous_color_div = c("#67001f", "#b2182b", "#d6604d",
                                              "#f4a582", "#fddbc7", "#f7f7f7",
                                              "#d1e5f0", "#92c5de", "#4393c3",
                                              "#2166ac", "#053061")) {
  th <- ggplot2::theme_classic()
  if (is.numeric(axis_text_x_angle) && axis_text_x_angle != 0) {
    th <- th %+replace%
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = axis_text_x_angle, hjust = 1,
        vjust = if (axis_text_x_angle >= 90) 0.5 else 1))
  }
  out <- th %+replace%
    ggplot2::theme(
      text = ggplot2::element_text(size = fontsize, family = family, face = face, colour = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, vjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = fontsize / 2, hjust = 0.5, vjust = 0.5),
      plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), unit),
      plot.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.background = ggplot2::element_rect(fill = NA, colour = NA),
      panel.border = if (panel_border) {
        ggplot2::element_rect(fill = NA, colour = "black")
      } else {
        ggplot2::element_blank()
      },
      panel.grid.major = if (panel_grid_major) {
        ggplot2::element_line(colour = "grey90", linewidth = 0.25)
      } else {
        ggplot2::element_blank()
      },
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.background = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key.size = ggplot2::unit(legend_size, unit),
      strip.background = ggplot2::element_rect(fill = NA, colour = NA),
      axis.line = if (panel_border) {
        ggplot2::element_blank()
      } else {
        ggplot2::element_line(linewidth = 0.25, colour = "black")
      },
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.25)
    )
  attr(out, "my_palettes") <- list(
    category_color_qua = category_color_qua,
    continuous_color_seq = continuous_color_seq,
    continuous_color_div = continuous_color_div
  )

  # ---- install the palettes as global default scales
  # remove previously installed scales so NULL cleanly disables a channel
  ours <- c("scale_fill_discrete", "scale_colour_discrete", "scale_color_discrete",
            "scale_fill_continuous", "scale_colour_continuous", "scale_color_continuous")
  rm(list = intersect(ours, ls(envir = .GlobalEnv, all.names = TRUE)),
     envir = .GlobalEnv)

  if (!is.null(category_color_qua)) {
    make_discrete <- function(aes) {
      force(aes)
      function(..., na.value = "grey50") {
        ggplot2::discrete_scale(
          aesthetics = aes,
          palette = function(n) adaptive_pal(n, category_color_qua),
          na.value = na.value,
          ...
        )
      }
    }
    assign("scale_fill_discrete",   make_discrete("fill"),   envir = .GlobalEnv)
    assign("scale_colour_discrete", make_discrete("colour"), envir = .GlobalEnv)
    assign("scale_color_discrete",  make_discrete("colour"), envir = .GlobalEnv)
  }

  seq_cols <- continuous_color_seq %||% continuous_color_div
  div_cols <- continuous_color_div %||% continuous_color_seq
  if (!is.null(seq_cols)) {
    # custom ScaleContinuous that remembers the data sign in train() and picks
    # the diverging (midpoint 0) or sequential gradient in map()
    make_continuous <- function(aes) {
      force(aes)
      function(..., na.value = NA, guide = "colourbar") {
        super <- ggplot2::ggproto(
          "ScaleContinuousAdaptive", ggplot2::ScaleContinuous,
          palette_seq = scales::gradient_n_pal(seq_cols),
          palette_div = scales::gradient_n_pal(div_cols),
          has_negative = FALSE,
          train = function(self, x) {
            self$has_negative <- self$has_negative || any(x < 0, na.rm = TRUE)
            ggplot2::ggproto_parent(ggplot2::ScaleContinuous, self)$train(x)
          },
          map = function(self, x, limits = self$get_limits()) {
            x <- self$oob(x, range = limits)
            x <- if (self$has_negative) {
              scales::rescale_mid(x, to = c(0, 1), from = limits, mid = 0)
            } else {
              scales::rescale(x, to = c(0, 1), from = limits)
            }
            uniq <- unique(x)
            pal <- if (self$has_negative) self$palette_div(uniq) else self$palette_seq(uniq)
            scaled <- pal[match(x, uniq)]
            if (anyNA(scaled)) scaled[is.na(scaled)] <- self$na.value
            scaled
          }
        )
        ggplot2::continuous_scale(
          aesthetics = aes,
          palette = function(x) x,
          na.value = na.value,
          guide = guide,
          super = super,
          ...
        )
      }
    }
    assign("scale_fill_continuous",   make_continuous("fill"),   envir = .GlobalEnv)
    assign("scale_colour_continuous", make_continuous("colour"), envir = .GlobalEnv)
    assign("scale_color_continuous",  make_continuous("colour"), envir = .GlobalEnv)
  }
  out
}


#' Use theme_hh() as the Global ggplot2 Theme
#'
#' Makes every subsequent ggplot2 plot use [theme_hh()] by default (equivalent
#' to calling [ggplot2::theme_set()]). Since [theme_hh()] also installs its
#' palettes as the global default scales, one call configures both the theme
#' and the colors. Passing `NULL` for a palette argument disables that color
#' channel (e.g. `category_color_qua = NULL` keeps ggplot2's default discrete
#' scale); passing `NULL` for all three disables all custom color setting.
#'
#' @param ... Arguments passed to [theme_hh()] (e.g. `fontsize = 8`,
#'   `category_color_qua = NULL`).
#'
#' @return Invisibly, the previously active theme.
#'
#' @examples
#' \dontrun{
#' set_theme_hh(fontsize = 8)
#' set_theme_hh(category_color_qua = NULL)                 # no custom discrete scale
#' set_theme_hh(category_color_qua = NULL,                 # disable all color setting
#'              continuous_color_seq = NULL,
#'              continuous_color_div = NULL)
#' }
set_theme_hh <- function(...) {
  invisible(ggplot2::theme_set(theme_hh(...)))
}


#' Adaptive Discrete Color Palette
#'
#' Generates `n` colors for a discrete scale from a base set of qualitative
#' colors. If `n` exceeds the number of base colors, additional colors are
#' created by interpolating through the base colors, so large category counts
#' never run out of colors.
#'
#' @param n Number of colors to generate.
#' @param colors Character vector of base colors; NULL returns `character()`.
#'
#' @return Character vector of `n` colors.
#'
#' @examples
#' adaptive_pal(4, c("#e41a1c", "#4daf4a", "#377eb8", "#984ea3"))
#' adaptive_pal(30, c("#e41a1c", "#4daf4a", "#377eb8"))
adaptive_pal <- function(n, colors) {
  if (n < 1 || is.null(colors)) return(character())
  colors <- as.character(colors)
  if (n <= length(colors)) unname(colors[seq_len(n)])
  else unname(grDevices::colorRampPalette(colors)(n))
}


#' Format Numbers with Automatic Scientific Notation
#'
#' Formats a numeric vector as strings, using scientific notation for values
#' below `accuracy` and fixed rounding otherwise (keeps label length consistent).
#'
#' @param v Numeric vector (or character/factor coercible to numeric).
#' @param accuracy Threshold below which scientific notation is used.
#'
#' @return Character vector.
#'
#' @examples
#' format_num_auto(c(1e-4, 0.5, 1234.567))
format_num_auto <- function(v, accuracy = 0.001) {
  if (is.factor(v)) v <- as.character(v)
  if (is.character(v)) v <- as.numeric(v)
  map_chr(v, \(x) {
    if (is.na(x)) NA_character_
    else if (abs(x) < accuracy) scales::label_scientific()(x)
    else as.character(round(x, log10(1 / accuracy)))
  })
}


#' Select Points from a ggplot via Polygon Lasso
#'
#' Renders the first layer of a ggplot as a base-R scatter plot and lets the
#' user click polygon vertices (right-click or ESC to finish). Returns the row
#' indices of the points inside the polygon.
#'
#' @param plot_obj A ggplot object (first layer must have x/y aesthetics).
#'
#' @return Integer vector of selected point indices (empty if cancelled).
ggplot_selector <- function(plot_obj) {
  if (!inherits(plot_obj, "ggplot")) stop("Input must be a ggplot object")

  pd <- ggplot_build(plot_obj)$data[[1]]
  x <- pd$x; y <- pd$y
  if (!length(x) || !length(y)) stop("No data points found in the plot")
  if (dev.cur() == 1) dev.new()

  plot(x, y,
       xlab = plot_obj$labels$x %||% "x", ylab = plot_obj$labels$y %||% "y",
       main = plot_obj$labels$title %||% "",
       col = pd$colour %||% pd$fill %||% "black",
       pch = pd$shape %||% 19,
       cex = (pd$size %||% 1) * 0.5)

  cat("Instructions:\n1. Click points to define polygon vertices\n",
      "2. Right-click or press ESC when done\n",
      "3. The polygon will be closed automatically\n\n")
  pts <- locator(type = "o", col = "red", lwd = 2)
  if (is.null(pts) || length(pts$x) < 3) {
    cat("No valid polygon selected (need at least 3 points)\n")
    return(integer(0))
  }

  px <- c(pts$x, pts$x[1]); py <- c(pts$y, pts$y[1])
  lines(px, py, col = "red", lwd = 2)

  n <- length(px) - 1
  inside <- vapply(seq_along(x), function(i) {
    in_poly <- FALSE
    j <- n
    for (k in seq_len(n)) {
      if (((py[k] > y[i]) != (py[j] > y[i])) &&
          (x[i] < (px[j] - px[k]) * (y[i] - py[k]) / (py[j] - py[k]) + px[k])) {
        in_poly <- !in_poly
      }
      j <- k
    }
    in_poly
  }, logical(1))

  sel <- which(inside)
  if (length(sel)) points(x[sel], y[sel], col = "darkblue", pch = 19, cex = 0.8)
  cat(sprintf("Selected %d points\n", length(sel)))
  sel
}



