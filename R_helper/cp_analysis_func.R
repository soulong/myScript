
# for .nd2 file extracted metadata channel info
get_channel_names_from_txt <- function(
    metadata_file,
    map_names=c("BFP"="405","GFP"="488","RFP"="561","IFP"="640")
) {
  
  # auto retrieve channel names from metadata_file
  metadata <- read_lines(metadata_file) %>% str_trim() %>% .[!(.=="")]
  metadata_iloc <- str_detect(metadata, ".nd2") %>% which()
  
  channel_names_list <- list()
  for(i in seq_along(metadata_iloc)) {
    i_loc <- metadata_iloc[i]
    if(i==length(metadata_iloc)) j_loc <- length(metadata) else j_loc <- metadata_iloc[i+1] - 1
    meta <- metadata[i_loc:j_loc]
    nd2_name <- metadata[i_loc] %>% basename() %>% str_replace(".nd2:","")
    ch_iloc <- meta %>% str_detect("Channel:") %>% which()
    mapped_names <- c()
    for(ch in ch_iloc) {
      ch_name <- str_replace(meta[ch+1], "Channel Name: ", "") %>% 
        str_replace_all(" Confocal", "") %>% 
        str_replace_all(" - .*", "") %>% 
        match(map_names) %>% names(map_names)[.] 
      ch_id <- str_replace(meta[ch], "-> Channel: ", "") %>% str_c("ch", .)%>% 
        set_names(ch_name)
      # to vector
      mapped_names <- c(mapped_names, ch_id)
      # print(mapped_names)
    }
    channel_names_list <- list(mapped_names) %>% 
      set_names(nd2_name) %>% 
      c(channel_names_list, .)
  }
  
  return(channel_names_list)
}



se <- function(x, ...) { sd(x, ...)/sqrt(length(x)) }



summarize_to_field <- function(
    cell,
    grouping_vars = c("Metadata_prefix", "Metadata_reporter", 
                      "Metadata_substrate", "Metadata_treat"),
    mean_keys = c("EquivalentDiameter", "Eccentricity",
                  "Intensity_MeanIntensity_", 
                  "Correlation_", 
                  "Granularity_"),
    sum_keys = c("Intensity_IntegratedIntensity_"),
    mean_func = function(x) mean(x, trim=0.05, na.rm=T), 
    sum_func = function(x) sum(x, na.rm=T)
) {
  
  # define summarize items
  mean_key_cols <- mean_keys %>% 
    map(\(x) str_subset(colnames(cell), x)) %>% 
    list_c()
  sum_key_cols <- sum_keys %>% 
    map(\(x) str_subset(colnames(cell), x)) %>% 
    list_c()
  
  # common grouping cols
  common_group_cols <- c(
    "Metadata_well", "Metadata_field", 
    "Metadata_position","Metadata_timepoint")
  
  # summarize
  per_field <- cell %>% 
    reframe(
      across(any_of(mean_key_cols), mean_func),
      across(any_of(sum_key_cols), sum_func),
      .by=any_of(c(grouping_vars, common_group_cols)))
  
  # count per filed at field level
  count_per_field <- cell %>% 
    reframe(
      count_per_field=n(),
      .by=any_of(c(grouping_vars, common_group_cols)))
  # add
  per_field <- right_join(count_per_field, per_field)

  # glimpse(per_field)

  return(per_field)
}



summarize_mean_sd <- function(
    per_image, 
    grouping_vars, 
    mean_trim = 0) {
  
  per_group <- per_image %>% 
    reframe(across(!contains(grouping_vars), 
                   list(mean=~ mean(.x, na.rm=T, trim=mean_trim),
                        sd=~ sd(.x, na.rm=T),
                        se=~ se(.x, na.rm=T))),
            .by=any_of(grouping_vars)) %>% 
    mutate(across(!contains(grouping_vars), ~ round(.x, 3)))
  
  return(per_group)
}



# normalize_data <- function(
#     data,
#     ctrl_colname="Metadata_name",
#     ctrl_value="DMSO",
#     grouping_vars=c("Metadata_prefix"),
#     use_median=FALSE,
#     mean_trim=0, # when calculate ctrl mean
#     const=0 # when divide for ratio, add const to get rid of error
# ) {
#   
#   # items for normalize
#   normCtrl_items <- colnames(data) %>% 
#     str_subset("^bg_intensity_", negate=T) %>% 
#     str_subset(".*Number$", negate=T) %>% 
#     str_subset("Metadata_", negate=T)
#   
#   # get ctrl value
#   if(use_median) {
#     # calculate median
#     print(str_glue("use group median as ctrl"))
#     data_ctrl <- data %>% 
#       reframe(across(any_of(normCtrl_items), ~ median(.x, na.rm=T)),
#               .by=all_of(grouping_vars)) %>% 
#       rename_with(~ str_c(.x, "_ctrl"), any_of(normCtrl_items))
#   } else {
#     # calculate mean from ref
#     print(str_glue("use group mean as ctrl"))
#     data_ctrl <- data %>% 
#       filter(!!as.name(ctrl_colname) == ctrl_value) %>% 
#       reframe(across(any_of(normCtrl_items), 
#                      ~ mean(.x, na.rm=T, trim=mean_trim)),
#               .by=all_of(grouping_vars)) %>% 
#       rename_with(~ str_c(.x, "_ctrl"), any_of(normCtrl_items))
#   }
#   # glimpse(data_ctrl)
#   
#   # get ratio for each item
#   data_norm <- left_join(data, data_ctrl)
#   for(item in normCtrl_items) {
#     data_norm[[str_glue("{item}")]] <- 
#       (data_norm[[str_glue("{item}")]] / 
#          data_norm[[str_glue("{item}_ctrl")]] + const) %>% 
#       round(4)
#   }
#   data_norm <- data_norm %>% select(!ends_with("_ctrl"))
#   # glimpse(data_norm)
#   
#   # calculate 3xMAD for ctrl
#   # usually for check ctrl variance and display control group variance om plot
#   if(F) {
#     data_norm_ctrl_mad <- data_norm %>% 
#       filter(!!as.name(ctrl_colname) == ctrl_value) %>% 
#       reframe(across(any_of(normCtrl_items), ~ mad(.x, na.rm=T)*3),
#               .by=all_of(grouping_vars))
#   }
#   
#   # ## get normalized data per_well
#   # per_well_normCtrl <- per_image_norm %>% 
#   #   reframe(across(!starts_with("Metadata"), 
#   #                  list(mean=~mean(.x, na.rm=T), 
#   #                       sd=~ sd(.x, na.rm=T),
#   #                       se=~ se(.x, na.rm=T))),
#   #           .by=any_of(
#   #             c(grouping_vars, "Metadata_well", "Metadata_timepoint"))) %>% 
#   #   mutate(across(!starts_with("Metadata"), ~ round(.x, 3)))
#   # # merge with p-value
#   # per_well_normCtrl <- per_well_normCtrl %>% 
#   #   {if(exists("stats_all")) left_join(., stats_all) else .}
#   # ## save data
#   # writexl::write_xlsx(
#   #   list(per_image=per_image_norm,
#   #        per_well=per_well_normCtrl),
#   #   str_glue("{Sys.Date()}_{obj}_summary_normCtrl.xlsx"))
#   
#   return(data_norm)
# }


calculate_foldchange <- function(
    data,
    ctrl_colname="Metadata_name",
    ctrl_value="DMSO",
    grouping_vars=c("Metadata_prefix"),
    use_median=FALSE,
    log_trans=F,
    mean_trim=0, # when calculate ctrl mean
    const=0 # when divide for ratio, add const to get rid of error
) {
  # items for normalize
  normCtrl_items <- colnames(data) %>% 
    str_subset("^bg_intensity_", negate=T) %>% 
    str_subset(".*Number$", negate=T) %>% 
    str_subset("Metadata_", negate=T) %>% 
    str_subset("_statistic", negate=T) %>% 
    str_subset("_pvalue", negate=T)
  # print(normCtrl_items)
  
  # get ctrl value
  if(use_median) {
    # calculate median
    print(str_glue("use group median as ctrl"))
    data_ctrl <- data %>% 
      reframe(across(any_of(normCtrl_items), ~ median(.x, na.rm=T)),
              .by=any_of(grouping_vars)) %>% 
      rename_with(~ str_c(.x, "_ctrl"), any_of(normCtrl_items))
  } else {
    # calculate mean from ref
    print(str_glue("use {ctrl_colname}={ctrl_value} as ctrl"))
    data_ctrl <- data %>% 
      filter(!!as.name(ctrl_colname) == ctrl_value) %>% 
      reframe(across(any_of(normCtrl_items), 
                     ~ mean(.x, na.rm=T, trim=mean_trim)),
              .by=all_of(grouping_vars)) %>% 
      rename_with(~ str_c(.x, "_ctrl"), any_of(normCtrl_items))
  }
  # glimpse(data_ctrl)
  
  # get ratio for each item
  data_fc <- left_join(data, data_ctrl)
  # glimpse(data_fc)
  
  # add constant value to minimize cal error
  for(item in normCtrl_items) {
    data_fc[[str_glue("{item}")]] <- (
      data_fc[[str_glue("{item}")]] / data_fc[[str_glue("{item}_ctrl")]] + const) %>% 
      { if(log_trans) round(log2(.), 4) else round(., 4) }
  }
  
  data_fc <- data_fc %>% select(!ends_with("_ctrl")) %>% 
    { if(log_trans) {
      rename_with(., \(x) str_c(x, "_log2FC"), any_of(normCtrl_items))
      } else rename_with(., \(x) str_c(x, "_FC"),  any_of(normCtrl_items)) }

    return(data_fc)
}



calculate_pvalue <- function(
    data,
    ctrl_colname="Metadata_name",
    ctrl_value="DMSO",
    grouping_vars=c("Metadata_prefix"),
    pvalue_items=c("Intensity_MeanIntensity_GFP"),
    use_median=FALSE
) {
  
  require(furrr)
  plan("multisession", workers=16)
  # plan(sequential)
  
  # process each plate one by one
  reformated <- data %>% 
    unite("uni_g", any_of(grouping_vars), sep="_#_")
  uni_g <- reformated %>% pull(uni_g) %>% unique()
  
  # computation should be done at per plate level
  test_func <- function(x, y) {
    tryCatch({
      t.test(x, y) %>% 
        broom::tidy() %>% 
        transmute(statistic=signif(statistic, 4),
                  pvalue=signif(p.value, 4)) },
      error=function(e) tibble(statistic=NA, pvalue=NA) )}
  
  stats_fun <- function(df, df_ctrl, stats_item, grouping_by) {
    groups <- df[[grouping_by]] %>% unique() %>% set_names(., nm=.)
    compare_base <- df_ctrl[[stats_item]]
    res_tidy <- furrr::future_map(
      groups, \(x) filter(df, !!as.name(grouping_by) == x) %>% 
        pull(!!as.name(stats_item)) %>% 
        test_func(., compare_base)) %>% 
      list_rbind(names_to=grouping_by)
    
    return(res_tidy)
  }
  
  stats_all <- list()
  for(item in pvalue_items) {
    stats <- list()
    for(g in uni_g) {
      print(str_glue("process {item} on: {g}"))
      s <- filter(reformated, uni_g == g)
      # calculate median
      if(use_median) {
        s_ctrl <- ifelse(nrow(s) > 1000, slice_sample(s, n=1000), s)
      } else {
        s_ctrl <- filter(s, !!as.name(ctrl_colname) == ctrl_value)
      }
      stats <- c(stats, stats_fun(s, s_ctrl, item, ctrl_colname) %>% list())
      stats <- stats %>% set_names(uni_g) %>% 
        list_rbind(names_to="uni_g") %>% 
        separate(uni_g, grouping_vars, sep="_#_", convert=F) %>% 
        rename("{item}_statistic" := statistic, 
               "{item}_pvalue" := pvalue)
      
      stats_all <- c(stats_all, list(stats))
    } }
  
  # rm(g, uni_g, reformated, s, s_ctrl, item, stats)
  
  stats_all <- reduce(stats_all, left_join)
  # data_pvalue <- left_join(data, stats_all)
  
  return(stats_all)
}








remove_outlier <- function(
    data, 
    grouping_vars, 
    keys=c("EquivalentDiameter", "Eccentricity",
           "Intensity_MeanIntensity_", "Intensity_IntegratedIntensity_"),
    coef=1.5
    ) {
  
  col_keys <- keys %>% 
    map(\(x) str_subset(colnames(data), x)) %>% 
    list_c()
  
  data_outlier <- data %>% 
    mutate(across(any_of(col_keys), 
                  \(x) rstatix::is_outlier(x, coef)),
           .by=any_of(grouping_vars)) %>% 
    dplyr::select(any_of(col_keys))
  
  kept_rows <- rowSums(data_outlier) < 1
  
  data_clean <- data[kept_rows, ]
  
  return(data_clean)
}






# read a single cellprofiler db
# only used in cellprofiler generated sqlite.db
# intensity will be multiply by 65535
read_sqlite <- function(
    db_path="cp_result.db", 
    remove_object_prefix=TRUE,
    # if subset is NULL, means all categories will be returned
    subset_areashape_category=c("AreaShape_Area", "_Eccentricity", "_EquivalentDiameter", "_Compactness", "_FormFactor", "_Solidity"),
    subset_intenstiy_category=c("IntegratedIntensity_", "MeanIntensity_"),
    # merge classify table (within db), generated from cellprofiler-analyst
    classify_table=NULL,
    # model, generated from fitted tidymodel workfolw
    model_workflow=NULL,
    # use tablename data as input data to model_workflow
    model_tablename="cell"
) {
  
  if(!file.exists(db_path)) stop("db path not existed")
  
  # transform A01, B13 to A1, B13 style
  transform_well_style <- function(well) {
    tmp <- str_extract(well, "(?<row>[A-Za-z])(?<column>[0-9]+)", group=1:2)
    # print(tmp)
    str_c(tmp[1], as.integer(tmp[2]))
  }
  
  # only used in cellprofiler generated sqlite.db
  average_glcm <- function(per_object) {
    # subset glcm columns
    glcm <- select(per_object, contains("Texture_"))
    glcm_colname_regrex <- "(?<prefix>[0-9A-Za-z]*)_?(?<item>Texture_\\w*)_(?<channel>ch[0-9]*)_(?<scale>[0-9]+)_[0-9]{2}_(?<level>[0-9]+)"
    # get colname-type
    glcm_compose <- colnames(glcm) %>% 
      set_names(., colnames(glcm)) %>% 
      map(~ str_match(.x, glcm_colname_regrex) %>% 
            .[-1] %>% str_c(collapse="_") %>% str_replace("^[_]", "")) %>%
      unlist() %>%
      enframe()
    # define a matrix
    glcm_mean <- matrix(nrow=nrow(glcm), 
                        ncol=length(unique(glcm_compose$value)))
    colnames(glcm_mean) <- unique(glcm_compose$value)
    # replace value by recurring columns
    for(uni in unique(glcm_compose$value)) {
      cols <- glcm_compose %>%
        filter(value==uni) %>%
        pull(name)
      glcm_mean[,uni] <- select(glcm, all_of(cols)) %>% rowMeans(na.rm=T)
    }
    # add glcm mean to original data
    per_object_rest <- select(per_object, all_of(setdiff(colnames(per_object), colnames(glcm))))
    per_object <- bind_cols(per_object_rest, glcm_mean)
    return(per_object)
  }
  
  # read a single object
  read_per_object <- function(conn, object_name) {
    
    object_prefix <- str_replace(object_name, "Per_", "") %>% 
      str_c("_")
    
    per_object <- tbl(conn, object_name) %>%
      dplyr::select(!c(contains("_Location_"))) %>%
      dplyr::select(!matches("(_X)|(_Y)$")) %>% 
      # filter areashape category
      { if(is.null(subset_areashape_category)) . else {
        dplyr::select(., !setdiff(contains("AreaShape_"), contains(subset_areashape_category)))
        # select(!matches(".*(QuartileIntensity)|(MADIntensity)|(MassDisplacement)|(MaxIntensity)|(MinIntensity)|(MedianIntensity)|(StdIntensity)|(IntensityEdge)_.*")) %>% 
      }  } %>% 
      # filter intensity category
      { if(is.null(subset_intenstiy_category)) . else {
        dplyr::select(., !setdiff(contains("Intensity_"), contains(subset_intenstiy_category)))
        # select(!matches(".*(QuartileIntensity)|(MADIntensity)|(MassDisplacement)|(MaxIntensity)|(MinIntensity)|(MedianIntensity)|(StdIntensity)|(IntensityEdge)_.*")) %>% 
      }  } %>% 
      # intensity to integer
      mutate(across(setdiff(contains("Intensity"), 
                            c(contains("Probabilities"), contains("ratio"))), 
                    \(x) round(x*65535))) %>% 
      # remove object prefix
      { if(remove_object_prefix) {
        rename_with(., \(x) str_replace(x, object_prefix, ""))
      } else . } %>%
      relocate(ends_with("_Object_Number"),
               contains("Parent_"),
               contains("Chirden_"),
               ends_with("_Count"),
               .after = ImageNumber) %>% 
      collect()
    
    # calculate mean of glcm 4-direction data
    if(any(str_detect(colnames(per_object), "Texture_"))) {
      print(">> average GLCM data")
      per_object <- average_glcm(per_object)
    }
    
    return(per_object)
  }
  
  
  # run here 
  conn <- db_path %>%
    # print() %>%
    DBI::dbConnect(RSQLite::SQLite(), dbname=.)
  
  tryCatch(
    {
      # read image metadata
      print('>>>>> read Per_Image')
      image <- tbl(conn, 'Per_Image') %>%
        dplyr::select(ImageNumber, 
                      starts_with("Image_Metadata_"),
                      starts_with("Image_Count_"),
                      starts_with("Image_Intensity_MeanIntensity"),
                      starts_with("Image_Intensity_Percentile"),
                      starts_with("Image_Granularity"),
                      starts_with("Image_Overlap"),
                      # Image_Width = Image_Width_ch1
        ) %>% 
        mutate(across(
          setdiff(contains("Intensity"), 
                  c(contains("Probabilities"), contains("ratio", ignore.case=T))), 
          \(x) round(x*65535))) %>% 
        rename_with(\(x) str_replace_all(x, "Image_Metadata", "Metadata")) %>% 
        collect() %>% 
        # transform A01, B13 to A1, B13 style
        {
          if("Metadata_well" %in% colnames(.)) {
            mutate(., Metadata_well=map_chr(Metadata_well, transform_well_style))
          } else {.}
        }
      # glimpse(image)
      
      
      # get all object tables
      object_names <- DBI::dbListTables(conn) %>% 
        str_subset("^Per_") %>% 
        str_subset("(Per_Experiment)|(Per_Image)", negate = T) %>% 
        set_names(., nm=str_replace_all(., "Per_", ""))
      print(str_glue("detect objects: {str_c(names(object_names), collapse=', ')}"))
      
      # get all object to list
      print('>>>>> read Per_Objects')
      objects <- map(object_names, \(x) read_per_object(conn, x))
      # glimpse(objects)
      
      
      # get classification table
      if(!is.null(classify_table)) {
        for(tab in classify_table) {
          print(str_glue('>>>>> read classify: {classify_table}'))
          objects[[classify_table]] <- 
            tbl(conn, classify_table) %>%
            dplyr::select(!class_number) %>% 
            collect()
        }
      }
      
      # apply classify model
      if(!is.null(model_workflow)) {
        print('>> fit model')
        objects[[model_tablename]] <- model_workflow %>% 
          augment(objects[[model_tablename]], type="prob") %>%
          rename_with(~ str_c("Metadata_model_", str_replace(.x, ".pred_", "prob_")), 
                      starts_with(".pred_")) %>%
          rename(Metadata_model_class=Metadata_model_prob_class) %>%
          relocate(starts_with("Metadata_model"), .before=ImageNumber)
      }
      
      
    }, finally={
      DBI::dbDisconnect(conn)
    }
  )
  
  # output
  return(list(image=image, objects=objects))
}



