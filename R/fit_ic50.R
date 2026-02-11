
library(dr4pl)
library(furrr)
# library(GRmetrics)
library(bestNormalize)
library(tidyverse)


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



plot_dose_response <- function(
    fitted_prepared, 
    ...,
    fact_var1 = NULL, # facet variables
    fact_var2 = NULL, # facet variables
    facet_scale = "free",
    facet_independent = "none", # facet independent 
    facet_axes = "all", # display all axes
    point_size = 0.75, # point size
    line_width = 1, # line width
    errorbar_width = 0.01,
    strip_size = 4,
    legend_key_size = 0.2,
    x_expand = c(0.05, 0.1),
    ylim = c(NA, NA),
    xlab="Dose [uM]", ylab="Response") {
  
  # fitted_prepared <- fitted_prepared[1:10,]
  
  # process additional aes
  args <- list(...)
  args <- lapply(args, function(x) if (rlang::is_string(x)) sym(x) else x)
  
  ## base plot
  p <- fitted_prepared %>% 
    unnest(stats) %>% 
    ggplot(aes(dose, mean, !!!args)) +
    geom_point(size=point_size) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                  linewidth=line_width/2, width=errorbar_width) +
    # ggplot(aes(x, y, color=Metadata_name + Metadata_model_class)) +
    geom_line(aes(x, y, group=uid), linewidth=line_width,
              data=unnest(fitted_prepared, pred)) +
    # # geom_ribbon(data=unnest(res2, pred), aes(x, y, ymin=ymin, ymax=ymax, fill=group), alpha=0.2) +
    # ggrepel::geom_text_repel(data=mutate(res2, pred2=map2(pred, group, ~ mutate(.x, label=ifelse(row_number()==nrow(.x), .y, NA)))) %>% unnest(pred2), 
    #                          aes(x,y,label=label), 
    #                          max.overlaps=10, nudge_x=1, nudge_y=1, size=2, 
    #                          segment.curvature = -0.1, segment.square = TRUE, segment.color = 'grey',
    #                          # arrow = arrow(length = unit(0.02, "npc")),
    #                          min.segment.length=0, show.legend=F, na.rm=T) +
    scale_x_log10(expand = expansion(mult = x_expand),
                  labels = function(lab) {
                    do.call(expression, lapply(paste(lab), function(x) bquote(.(x)))) })
  
  ## facet
  if(!is.null(fact_var1)) {
    if(!is.null(fact_var2)) {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          str_c(fact_var1, collapse = "+"), " ~ ", str_c(fact_var2, collapse = "+"))),
        scales = facet_scale,
        independent = facet_independent,
        axes  = facet_axes) 
    } else {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          str_c(fact_var1, collapse = "+"), " ~ ", ".")),
        scales = facet_scale,
        independent = facet_independent,
        axes  = facet_axes) }
  } else {
    if(!is.null(fact_var2)) {
      p <- p + ggh4x::facet_grid2(
        as.formula(str_c(
          ".", " ~ ", str_c(fact_var2, collapse = "+"))),
        scales = facet_scale,
        independent = facet_independent,
        axes  = facet_axes) }
  }
  
  ## theme
  p <- p +
    coord_cartesian(ylim = ylim) +
    labs(x=xlab, y=ylab)
  
  return(p)
}
