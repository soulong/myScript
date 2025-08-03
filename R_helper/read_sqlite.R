

# read a single cellprofiler db ---------------
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


