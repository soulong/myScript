

library(ggplot2)
library(ggsci)

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


my_theme <- function(angle=45,
                     fontsize = 7, 
                     family = NULL, 
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
      plot.subtitle = ggplot2::element_text(fontsize/2, hjust = 0.5, vjust = 0.5),
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




se <- function(x) {sd(x, na.rm=T)/sqrt(length(x)-1)}


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





cp_analysis_plot <- function(
    data, # # each row is a summarized data, mean, sd, se
    x_axis, # x variable
    y_axis, # y variable
    error_axis = "se", # se variable
    ..., # other aes mapping
    # aes_mapping = NULL, # other aes mapping
    # color_var = NULL, # fill variable
    fact_var1 = NULL, # facet variables
    fact_var2 = NULL, # facet variables
    facet_scale = "fixed", # facet scale
    facet_independent = "none", # facet independent 
    facet_axes = "all", # display all axes
    fact_drop = FALSE, # discard item if not existed in facet
    plot_type = c("line", "dot", "dotline", "bar"), # plot type
    dodge_preserve = "single",
    point_size = 0.5, # point size
    line_width = 0.5, # line width
    bar_width = 0.8, # bar width
    show_legend = F, # show legend
    x_size = 6,
    x_rotatin = 0,
    x_lab = NULL,
    y_lab = NULL,
    raster = F,
    clip_y_range = "atuo"
) {
  
  # process additional aes
  args <- list(...)
  args <- lapply(args, function(x) if (rlang::is_string(x)) sym(x) else x)
  # print(args)
  
  plot_type = plot_type[1]
  
  ## base plot
  if(is.null(error_axis)) {
    p <- ggplot(data, 
                aes(x = .data[[x_axis]], y = .data[[y_axis]]))
  } else {
    p <- ggplot(data, 
                aes(x = .data[[x_axis]], y = .data[[y_axis]],
                    ymin = .data[[y_axis]] - .data[[error_axis]], 
                    ymax = .data[[y_axis]] + .data[[error_axis]] ))
  }
  
  ## plot type
  # dot
  if(plot_type == "dot") {
    p <- p + 
      geom_pointrange(aes(!!!args),
                      position = position_dodge(width = bar_width+0.05, preserve = dodge_preserve),
                      size = point_size, 
                      linewidth = line_width,
                      show.legend = show_legend) }
  # line
  if(plot_type == "line") {
    p <- p + 
      geom_line(aes(!!!args),
                # position = position_dodge(width = bar_width+0.05, preserve = "single"),
                linewidth = line_width, show.legend = show_legend)
    if(!is.null(error_axis)) {
      p <- p +
        geom_errorbar(aes(!!!args),
                      # position = position_dodge(width = bar_width+0.05, preserve = "single"), 
                      width = round(bar_width*0.6, 1),
                      linewidth = line_width,
                      show.legend = show_legend) }
  }
  
  # dotline
  if(plot_type == "dotline") {
    p <- p + 
      geom_pointrange(aes(!!!args),
                      size = point_size, 
                      linewidth = line_width,
                      show.legend = show_legend) +
      geom_line(aes(!!!args),
                linewidth = line_width, show.legend = show_legend) }
  
  # bar
  if(plot_type == "bar") {
    p <- p + 
      geom_col(aes(!!!args),
               # stat = "identity",
               color = "transparent",
               position = position_dodge(width = bar_width+0.05, preserve = dodge_preserve),
               width = bar_width, 
               show.legend = show_legend)
    if(!is.null(error_axis)) {
      p <- p +
        geom_errorbar(aes(!!!args),
                      position = position_dodge(width = bar_width+0.05, preserve = dodge_preserve), 
                      width = round(bar_width*0.6, 1), 
                      linewidth = line_width,
                      show.legend = show_legend) }
  }
  
  
  ## facet
  if(!is.null(fact_var1)) {
    if(!is.null(fact_var2)) {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          str_c(fact_var1, collapse = "+"), " ~ ", str_c(fact_var2, collapse = "+"))),
        scales = facet_scale,
        independent = facet_independent,
        axes = facet_axes,
        drop = fact_drop) 
    } else {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          str_c(fact_var1, collapse = "+"), " ~ ", ".")),
        scales = facet_scale,
        independent = facet_independent,
        axes = facet_axes,
        drop = fact_drop) 
    }
  } else {
    if(!is.null(fact_var2)) {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          ".", " ~ ", str_c(fact_var2, collapse = "+"))),
        scales = facet_scale,
        independent = facet_independent,
        axes = facet_axes,
        drop = fact_drop)
    }
  }
  
  # set y range
  if(is.character(clip_y_range)[1] & length(clip_y_range) == 1) {
    if(!is.null(error_axis)) {
      y_range <- c(min(data[[y_axis]] - data[[error_axis]]),
                   max(data[[y_axis]] + data[[error_axis]]))
    } else {y_range <- range(data[[y_axis]])}
    p <- p + coord_cartesian(expand=T, ylim=c(y_range[1], y_range[2]))
  } else {
    if(is.numeric(clip_y_range[1]) & length(clip_y_range) == 2)
      p <- p + coord_cartesian(expand=T, ylim=clip_y_range)
  }
  
  ## theme
  if(!is.null(x_lab)) p <- p + labs(x = x_lab)
  if(!is.null(y_lab)) p <- p + labs(y = y_lab)
  p <- p + 
    theme_bw(8) +
    theme(strip.background = element_rect(fill=NA), 
          strip.text = element_text(size = x_size),
          axis.text.x = element_text(
            size = x_size, angle = x_rotatin, 
            vjust=ifelse(x_rotatin>=45, 0.5, 0),
            hjust=ifelse(x_rotatin<=45, 0.5, 1)
          ),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.key.size = unit(0.4, "mm"))
  
  
  # rasterize
  if(is.integer(raster)) {
    p <- ggrastr::rasterize(p, dpi = raster)
  }
  
  return(p)
}





cp_analysis_plot_individual <- function(
    data, # each row is a individual measurement
    x_axis, # x variable
    y_axis, # y variable
    ..., # other aes mapping
    # aes_mapping = NULL, # other aes mapping
    # color_var = NULL, # fill variable
    fact_var1 = NULL, # facet variables
    fact_var2 = NULL, # facet variables
    facet_scale = "fixed", # facet scale
    facet_independent = "none", # facet independent 
    facet_axes = "all", # display all axes
    fact_drop = FALSE,
    plot_type = c("boxplot", "ridge"), # plot type
    dodge_preserve = "single",
    na_rm = TRUE,
    box_width = 0.9,
    box_staple_width = 0.5,
    box_outlier = FALSE,
    box_outlier_color = "grey50",
    box_outlier_size = 0.5,
    linewidth = 0.5,
    ridge_scale = 1, 
    fill_alpha=1,
    show_legend = F, # show legend
    x_size = 6,
    y_size = 6,
    x_rotatin = 0,
    x_lab = NULL,
    y_lab = NULL,
    raster = F
) {
  
  # process additional aes
  args <- list(...)
  args <- lapply(args, function(x) if (rlang::is_string(x)) sym(x) else x)
  # print(args)
  
  
  ## base plot
  p <- ggplot(data, aes(x = .data[[x_axis]], y = .data[[y_axis]]))
  
  ## plot type
  plot_type = plot_type[1]
  
  # boxplot
  if(plot_type == "boxplot") {
    p <- p + geom_boxplot(
      aes(!!!args), linewidth=linewidth, alpha = fill_alpha, outliers = box_outlier, 
      outlier.color = box_outlier_color, outlier.size = box_outlier_size,
      # position = position_dodge(preserve = dodge_preserve),
      width = box_width, na.rm = na_rm, show.legend = show_legend) }
  
  # ridge
  if(plot_type == "ridge") {
    p <- p + ggridges::geom_density_ridges(
      aes(!!!args), scale = ridge_scale, linewidth=linewidth, alpha = fill_alpha, 
      # position = position_dodge(preserve = dodge_preserve),
      na.rm = na_rm, show.legend = show_legend) }
  
  ## facet
  if(!is.null(fact_var1)) {
    if(!is.null(fact_var2)) {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          str_c(fact_var1, collapse = "+"), " ~ ", str_c(fact_var2, collapse = "+"))),
        scales = facet_scale,
        independent = facet_independent,
        axes = facet_axes,
        drop = fact_drop)
    } else {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          str_c(fact_var1, collapse = "+"), " ~ ", ".")),
        scales = facet_scale,
        independent = facet_independent,
        axes = facet_axes,
        drop = fact_drop) 
    } 
  } else {
    if(!is.null(fact_var2)) {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          ".", " ~ ", str_c(fact_var2, collapse = "+"))),
        scales = facet_scale,
        independent = facet_independent,
        axes = facet_axes,
        drop = fact_drop)
    }  }
  
  ## theme
  if(!is.null(x_lab)) p <- p + labs(x = x_lab)
  if(!is.null(y_lab)) p <- p + labs(y = y_lab)
  p <- p + 
    theme_bw(8) +
    theme(strip.background = element_rect(fill=NA), 
          strip.text = element_text(size = x_size),
          axis.text.x = element_text(
            size = x_size, angle = x_rotatin, 
            vjust=ifelse(x_rotatin>=45, 0.5, 0),
            hjust=ifelse(x_rotatin<=45, 0.5, 1)
          ),
          panel.grid = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          legend.key.size = unit(0.4, "mm"))
  
  # rasterize
  if(is.integer(raster)) {
    p <- ggrastr::rasterize(p, dpi = raster) }
  
  return(p)
}

