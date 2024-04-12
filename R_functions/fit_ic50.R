
library(dr4pl)
library(lemon)
library(furrr)
# library(GRmetrics)


# func
model_func <- function(data, ...) {
  dr4pl(response ~ dose, 
        data = data, 
        # method.init = "logistic",
        # method.robust = "squared",
        # lowerl = c(theta_4 = 0),
        ...)
}

coef_tidy <- function(model) {
  ci <- summary(model) %>% 
    .[["coefficients"]] %>%
    as.data.frame()
  ci <- mutate(ci, CI = map_dbl(seq_along(ci), 
                                ~ (ci[[4]][.x] - ci[[3]][.x])/2))
  # convert logIC50to IC50
  ci_ic50 <- ci[2, ] %>% 10^.
  rownames(ci_ic50) <- "IC50"
  # combine
  out <- bind_rows(ci[-2, ], ci_ic50) %>% 
    .[c(1,4,2,3), c(1,5)] %>%
    mutate(across(everything(), format, digits=3, scientific=T)) %>%
    mutate(across(everything(), str_trim)) %>%
    # mutate(across(everything(), as.numeric)) %>%
    as_tibble(rownames = "type") %>%
    mutate(type = str_replace_all(type, "Limit", "")) %>%
    pivot_wider(names_from = type, names_glue = "{type}_{.value}", values_from = c(Estimate, CI))
  return(out)
}

pred_func <- function(model, se=F, normalize_residual=T, level=0.95, nboot=200) {
  
  # model = res$model[[2]]
  range <- sort(unique(model$data$Dose)) %>% na.omit()
  from <- ifelse(range[1]<=0, 0.8*range[2], 0.8*range[1])
  to <- 1.2*last(range)
  xseq <- exp(seq(from=log(from), to=log(to), length.out=200))
  
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
  
  # parellal worker
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

fit_dose_response <- function(df, dose, response, group=NULL, n_worker=12, ...) {
  
  plan(multisession, workers = n_worker)
  
  data <- ungroup(df) %>% 
    transmute(dose = {{dose}}, #Metadata_Concentration, 
              response = {{response}}, #spots_count,
              group = {{group}}) # Metadata_Compound
  
  # glimpse(data)
  
  # model data
  data <- group_by(data, group) %>%
    nest() %>%
    mutate(model = map(data, model_func, ...),
           convergence = map_lgl(model, ~ .x[["convergence"]]),
           method_robust = map_chr(model, ~ .x[["method.robust"]]),
           coef = map(model, coef_tidy)
    )
  # glimpse(data)
  # return(data)
  
  # subset converged model for prediction
  data_converged <- data %>%
    # filter(isTRUE(convergence)) %>%
    mutate(pred = map(model, pred_func))
  
  return(data_converged)
}

plot_dose_response <- function(res, ncol=4, save=F, name="",
                               type="single", scales="free_x",
                               xlab="Dose [uM]", ylab="Viablity by Cell Titer Glo (%)") {
  
  # calculate mean and se
  res2 <- group_by(res, group) %>%
    mutate(stats = map(
      data, 
      ~ group_by(.x, dose) %>% 
        summarise(mean=mean(response), se=sd(response)/sqrt(n()), .groups = "drop"))) %>%
    mutate(top_inhit = map(pred, function(x) round(100*(1-min(x$y)/max(x$y)), 0))) %>%
    unnest(coef) %>%
    mutate(group_label = str_glue("{group}\n IC50: {str_replace(IC50_Estimate,'[+]','')} ± {str_replace(IC50_CI,'[+]','')}\n Top Inhibition: {top_inhit}%"))
  
  # align plot by single group
  if(type=="single") {
    p <- 
      unnest(res2, stats) %>%
      ggplot(aes(dose, mean)) +
      geom_point(size=0.75, color="grey25") +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), color="grey25", width=0.1) +
      geom_line(data=unnest(res2, pred), aes(x, y), color="#b2182b", size=1) +
      # geom_ribbon(data=unnest(res2, pred), aes(x, y, ymin=ymin, ymax=ymax), fill="#b2182b", alpha=0.2) +
      scale_x_log10(expand = expansion(mult = c(0.05, 0.05)),
                    labels = function(lab) {
                      do.call(expression, lapply(paste(lab), function(x) bquote(.(x)))) }) +
      facet_rep_wrap(~ group_label, repeat.tick.labels=T, scales=scales, ncol=ncol) +
      # coord_trans(x="log10") + # why failed?
      ggpubr::theme_pubr(8, border = T, x.text.angle = 0) +
      coord_cartesian(ylim = c(0, NA)) +
      labs(x=xlab, y=ylab)
    
    if(save) {
      ggsave2(str_c(Sys.Date(), "_", name, ".pdf"), p,
              width = ncol*2, height = ceiling(length(res$group)/ncol) * 2) }
    
  } else {
    p <-
      unnest(res2, pred) %>%
      # ggplot(aes(dose, mean, color=group)) +
      # geom_point(size=1) +
      # geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.075) +
      ggplot(aes(x, y, color=group)) +
      geom_line(size=1, show.legend = F) +
      # geom_ribbon(data=unnest(res2, pred), aes(x, y, ymin=ymin, ymax=ymax, fill=group), alpha=0.2) +
      ggrepel::geom_text_repel(data=mutate(res2, pred2=map2(pred, group, ~ mutate(.x, label=ifelse(row_number()==nrow(.x), .y, NA)))) %>% unnest(pred2), 
                               aes(x,y,label=label), 
                               max.overlaps=10, nudge_x=1, nudge_y=1, size=2, 
                               segment.curvature = -0.1, segment.square = TRUE, segment.color = 'grey',
                               # arrow = arrow(length = unit(0.02, "npc")),
                               min.segment.length=0, show.legend=F, na.rm=T) +
      scale_x_log10(expand = expansion(mult = c(0.05, 0.25)),
                    labels = function(lab) {
                      do.call(expression, lapply(paste(lab), function(x) bquote(.(x)))) }) +
      # facet_rep_wrap(~ group_label, repeat.tick.labels = T, scales = scales) +
      # coord_cartesian(ylim = c(0, 5000)) +
      ggpubr::theme_pubr(8, border = T, x.text.angle = 0) +
      coord_cartesian(ylim = c(0, NA)) +
      labs(x=xlab, y=ylab)
    
    if(save) {
      ggsave2(str_c(Sys.Date(), "_", name, "_merge.pdf"), p,
              width = 4.5, height = 3.5) }
  }
  
  return(list(data=res2, plot=p))
}

