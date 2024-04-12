

# read a single cellprofiler db ---------------
# only used in cellprofiler generated sqlite.db
# intensity will be multiply by 65535
read_sqlite <- function(db_path="cp_result.db", 
                        classify_table=NULL,
                        remove_object_prefix=FALSE) {
  
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
      dplyr::select(-contains("AreaShape_Orientation"), 
                    -contains("AreaShape_EulerNumber")) %>% 
      # filter intensity category # (IntegratedIntensity)
      select(!matches(".*(QuartileIntensity)|(MADIntensity)|(MassDisplacement)|(MaxIntensity)|(MinIntensity)|(MedianIntensity)|(StdIntensity)|(IntensityEdge)_.*")) %>% 
      mutate(across(setdiff(contains("Intensity"), contains("Probabilities")), 
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
                      Image_Width = Image_Width_ch1) %>% 
        mutate(across(contains("Intensity"), \(x) round(x*65535))) %>% 
        rename_with(\(x) str_replace_all(x, "Image_Metadata", "Metadata")) %>% 
        collect() %>% 
        # transform A01, B13 to A1, B13 style
        mutate(Metadata_well=map_chr(Metadata_well, transform_well_style))
      # glimpse(image)
      
      
      # get all object tables
      object_names <- DBI::dbListTables(conn) %>% 
        str_subset("^Per_") %>% 
        str_subset("(Per_Experiment)|(Per_Image)", negate = T) %>% 
        set_names(., nm=str_replace_all(., "Per_", ""))
      
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
      
      
    }, finally={
      DBI::dbDisconnect(conn)
    }
  )
  
  # output
  return(list(image=image, objects=objects))
}





# # read batch cellprofiler data ---------------
# # merge, clean and aggregate data
# # only used in cellprofiler generated sqlite.db
# # intensity will be multiply by 65535
# read_sqlite_batch <- function(measurement_dir, 
#                               cp_db_filename='cp_result.db',
#                               metadata=NULL,
#                               trainset_filename=NULL, 
#                               train_model_wf=NULL) {
#   
#   # transform A01, B13 to A1, B13 style
#   transform_well_style <- function(well) {
#     tmp <- str_extract(well, "(?<row>[A-Za-z])(?<column>[0-9]+)", group=1:2)
#     # print(tmp)
#     str_c(tmp[1], as.integer(tmp[2]))
#   }
#   
#   # only used in cellprofiler generated sqlite.db
#   average_glcm <- function(per_object) {
#     # subset glcm columns
#     glcm <- select(per_object, contains("Texture_"))
#     glcm_colname_regrex <- "(?<prefix>[0-9A-Za-z]*)_?(?<item>Texture_\\w*)_(?<channel>ch[0-9]*)_(?<scale>[0-9]+)_[0-9]{2}_(?<level>[0-9]+)"
#     # get colname-type
#     glcm_compose <- colnames(glcm) %>% 
#       set_names(., colnames(glcm)) %>% 
#       map(~ str_match(.x, glcm_colname_regrex) %>% 
#             .[-1] %>% str_c(collapse="_") %>% str_replace("^[_]", "")) %>%
#       unlist() %>%
#       enframe()
#     # define a matrix
#     glcm_mean <- matrix(nrow=nrow(glcm), 
#                         ncol=length(unique(glcm_compose$value)))
#     colnames(glcm_mean) <- unique(glcm_compose$value)
#     # replace value by recurring columns
#     for(uni in unique(glcm_compose$value)) {
#       cols <- glcm_compose %>%
#         filter(value==uni) %>%
#         pull(name)
#       glcm_mean[,uni] <- select(glcm, all_of(cols)) %>% rowMeans(na.rm=T)
#     }
#     # add glcm mean to original data
#     per_object_rest <- select(per_object, all_of(setdiff(colnames(per_object), colnames(glcm))))
#     per_object <- bind_cols(per_object_rest, glcm_mean)
#     return(per_object)
#   }
#   
#   print(str_glue(">>> processing:\n {measurement_dir}"))
#   db_path <- file.path(measurement_dir, cp_db_filename)
#   if(!file.exists(db_path)) {
#     print("no db file found, skip")
#     return(NULL)
#   }
#   
#   conn <- file.path(measurement_dir, cp_db_filename) %>%
#     DBI::dbConnect(RSQLite::SQLite(), dbname=.)
#   # DBI::dbListTables(conn)
#   tryCatch(
#     {
#       print('>> read table: Per_Image')
#       # read per_image data
#       per_image <- tbl(conn, 'Per_Image') %>%
#         select(starts_with('Image_Metadata_'), ImageNumber) %>% # glimpse()
#         rename_with(~ str_replace_all(.x, 'Image_Metadata', 'Metadata'), 
#                     starts_with('Image_Metadata')) %>%
#         # keep informative columns only 
#         # select(!c(Metadata_row, Metadata_column, Metadata_stack)) %>%
#         select(Metadata_directory, Metadata_well, Metadata_field, ImageNumber) %>%
#         collect() %>%
#         # clean directory path, must be Operetta style
#         mutate(Metadata_directory=Metadata_directory %>%
#                  str_split('[/]', simplify=T) %>% 
#                  str_subset('.*-Measurement [0-9]') %>%
#                  str_split('(__)', simplify=T) %>%
#                  .[1]) %>%
#         # transform A01, B13 to A1, B13 style
#         mutate(Metadata_well=map_chr(Metadata_well, transform_well_style))
#       
#       # read per_object data
#       print('>> read table: Per_Object')
#       per_object <- tbl(conn, 'Per_Object') %>% # glimpse()
#         # remove non-informative columns
#         select(!c(contains('_Location_'))) %>%
#         select(!matches("(_X)|(_Y)$")) %>%
#         collect()
#       # remove feature object prefix like cell_, cell_resized_
#       object_prefix <- colnames(per_object) %>% 
#         str_subset('_AreaShape_Area$') %>%
#         str_replace('AreaShape_Area', '') %>%
#         str_c("^", .)
#       per_object <- rename_with(per_object, 
#                                 ~ str_replace(.x, object_prefix, ''), everything())
#       
#       # calculate mean of glcm 4-direction data
#       if(any(str_detect(colnames(per_object), "Texture_"))) {
#         print(">> average GLCM data")
#         per_object <- average_glcm(per_object)
#       }
#       
#       # read training set
#       if(!is.null(trainset_filename)) {
#         print('>> read training set')
#         trainset <- str_glue("{measurement_dir}/{trainset_filename}") %>%
#           read_csv() %>%
#           select(ImageNumber, ObjectNumber, Class)
#         # merge with per_object
#         per_object <- right_join(trainset, per_object)
#       }
#       
#       # fit model or not
#       if(!is.null(train_model_wf)) {
#         print('>> fit barcode model')
#         per_object <- augment(train_model_wf, per_object, type="prob") %>%
#           rename_with(~ str_c("Metadata_barcode_", str_replace(.x, ".pred_", "prob_")), 
#                       .pred_class:last_col()) %>%
#           rename(Metadata_barcode=Metadata_barcode_prob_class) %>%
#           relocate(starts_with("Metadata_barcode"), .before=ImageNumber) %>%
#           mutate(Metadata_barcode=as.character(Metadata_barcode))
#       }
#       
#       # merge image metadata and per object
#       print('>> merge oject with image')
#       raw <- right_join(per_image, per_object, 
#                         by=join_by(ImageNumber), multiple = "all") %>%
#         relocate(ImageNumber, .before=ObjectNumber)
#       
#       # clean data
#       print('>> clean object')
#       raw <- raw %>% 
#         select(Metadata_directory:Metadata_field, 
#                Metadata_barcode:Metadata_barcode_prob_mS,
#                AreaShape_Area, AreaShape_Eccentricity, 
#                AreaShape_Extent, AreaShape_FormFactor, AreaShape_Solidity,
#                Granularity_1_ch1, Granularity_2_ch1, Granularity_3_ch1,
#                Intensity_MeanIntensity_ch1, Intensity_MeanIntensity_ch2,
#                RadialDistribution_MeanFrac_ch1_1of5:RadialDistribution_MeanFrac_ch1_5of5
#         ) %>% 
#         mutate(
#           Intensity_MeanIntensity_ch1=as.integer(Intensity_MeanIntensity_ch1*65535),
#           Intensity_MeanIntensity_ch2=as.integer(Intensity_MeanIntensity_ch2*65535),
#           RadialDistribution_MeanFrac_ch1_1_45 = round(
#             RadialDistribution_MeanFrac_ch1_1of5 /
#               (RadialDistribution_MeanFrac_ch1_4of5 + 
#                  RadialDistribution_MeanFrac_ch1_5of5), 1)) %>% 
#         select(!(RadialDistribution_MeanFrac_ch1_1of5:RadialDistribution_MeanFrac_ch1_5of5))
#       
#       # merge metadata
#       if(!is.null(metadata)) {
#         print('>> addd metadata')
#         raw <- right_join(metadata, raw)
#       }
#       
#       
#     }, finally={
#       DBI::dbDisconnect(conn)
#     }
#   )
#   
#   return(raw)
# }


