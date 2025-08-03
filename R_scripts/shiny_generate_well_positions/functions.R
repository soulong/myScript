
library(tidyverse)
library(xml2)
library(progress)

## xml func --------------------------------------------------------------
xml2df <- function(doc) {
  if(!is(doc, "xml_document")) doc <- read_xml(doc)
  
  attrs <- c("bChecked","strName",
             "dXPosition","dYPosition","dZPosition",
             "dPFSOffset") %>% 
    set_names(., nm=.)
  
  df <- map(attrs, 
            \(x) xml_find_all(doc, str_glue("//{x}")) %>% 
              map_chr(\(x) xml_attr(x, "value"))) %>% 
    bind_rows()
  
  df <- df %>% 
    mutate(PositionID = xml2::as_list(doc) %>% 
             .[["variant"]] %>% 
             .[["no_name"]] %>% 
             names() %>% 
             str_subset("^Point"),
           .before = 1)
  
  return(df)
}


remove_position_node <- function(doc, tag_name=NULL) {
  if(!is(doc, "xml_document")) stop("doc must be xml_document")
  
  if(is.null(tag_name)) {
    xml_find_all(doc, ".//*[@runtype='NDSetupMultipointListItem']") %>% 
      xml_remove()
  } else {
    xml_find_all(doc, str_glue("//*/{tag_name}")) %>% 
      xml_remove()
  }
  
  return(NULL)
}


add_position_node <- function(doc, 
                              meta=c("Point00000", # PositionID
                                     "true", # bChecked
                                     "A0", # strName
                                     "1", # dXPosition
                                     "2", # dYPosition
                                     "3", # dZPosition
                                     "0" # dPFSOffset
                              )) {
  if(!is(doc, "xml_document")) stop("doc must be xml_document")
  
  xml_child(doc) %>% 
    xml_add_child(meta[1], runtype="NDSetupMultipointListItem") %>% 
    xml_add_child('bChecked', runtype="bool", value=meta[2]) %>% 
    xml_add_sibling('strName', runtype="CLxStringW", value=meta[3]) %>% 
    xml_add_sibling('dXPosition', runtype="double", value=meta[4]) %>% 
    xml_add_sibling('dYPosition', runtype="double", value=meta[5]) %>% 
    xml_add_sibling('dZPosition', runtype="double", value=meta[6]) %>% 
    xml_add_sibling('dPFSOffset', runtype="double", value=meta[7]) %>% 
    xml_root()
}


create_empty_xml <- function(
    target=NULL, 
    source="template.xml"
) {
  
  if(is.null(target)) target <- str_glue("{Sys.Date()}_new.xml")
  
  # plate meta here
  doc <- read_xml(source)
  write_xml(doc, target)
  
  doc_new <- read_xml(target)
  remove_position_node(doc_new)
  write_xml(doc_new, target)
  
  cat(">> create new empty file:", target, "\n")
  
  return(list(doc=doc_new, fname=target))
}


save_xml <- function(new_points_df, fname=NULL) {
  res <- create_empty_xml(fname)
  fname <- res$fname
  doc <- res$doc
  
  # write to new xml
  pb <- progress_bar$new(
    format = " [:bar] :current/:total (:percent) eta: :eta",
    total = nrow(new_points_df), 
    width = 60,
    clear = FALSE
    )
  for(i in seq_len(nrow(new_points_df))) {
    pb$tick()
    new_points_df[i, ] %>% 
      flatten_chr() %>% 
      add_position_node(doc, meta=.)
  }

  # print(xml2df(doc))
  write_xml(doc, fname)
}




## generate NxN field --------------------------------------------------------------
shift_well_coords <- function(n_row=16, 
                              n_col=24, 
                              well_interval_x=50, 
                              well_interval_y=50, 
                              x_direct=-1, y_direct=1) {
  # Loop over the wells
  shift <- matrix(nrow=n_row*n_col, ncol=2)
  colnames(shift) <- c("x_shift", "y_shift")
  rownames(shift) <- seq_len(nrow(shift))
  
  i <- 1
  for (row in seq_len(n_row)) {
    for (col in seq_len(n_col)) {
      shift[i, ] <- c(x_direct *(col - 1) * well_interval_x, 
                      y_direct *(row - 1) * well_interval_y)
      rownames(shift)[i] <- paste0(LETTERS[row], col)
      i <- i + 1
    }
  }
  
  return(as_tibble(shift, rownames="well"))
}

# generate N x N fixed square field with field interval
generate_well_field_coords <- function(
    diameter, center_x, center_y, 
    field_width, field_interval, 
    N) {
  
  # Create data frame for the squares
  fields <- data.frame(
    x = numeric(0),
    y = numeric(0),
    width = numeric(0)
  )
  
  # Calculate max NxN squares
  max_N <- floor(diameter / (field_width + field_interval))
  if (N > max_N) N <- max_N
  seqs <- seq_len(N) - mean(seq_len(N))
  
  # Calculate square coordinates
  for (i in seqs) {
    for (j in seqs) {
      x <- center_x + i * (field_width + field_interval)
      y <- center_y + j * (field_width + field_interval)
      
      # Check if the square is within the circle
      if (sqrt((x - center_x)^2 + (y - center_y)^2) < diameter/2) {
        fields <- rbind(fields, data.frame(x = x, y = y, width = field_width))
      }
    }
  }
  
  # rownames(fields) <- as.character(seq_len(N*N))
  
  return(as_tibble(fields))
}



optim_fields <- function(fields) {
  temp_1 <- fields %>% 
    distinct(x) %>% 
    mutate(field_row_odd=ifelse(row_number()%%2==0, 2, 1))
  
  fields_optim <- left_join(fields, temp_1,
                            by = join_by(x)) %>% 
    mutate(field_x_reorder =
             ifelse(field_row_odd==1, x, x),
           field_y_reorder=
             ifelse(field_row_odd==1, y, -y)
    ) %>% 
    arrange(field_x_reorder, field_y_reorder) %>% 
    mutate(field = row_number())
  
  return(fields_optim)
}


# well_fields=well_fields_filtered
optim_wells <- function(well_fields, 
                        by=c("lr","ll","tb","tt"),
                        clip_zero=FALSE
                        ) {
  by = by[1]
  ## arrange well positions
  # to solve auto order problem when loading to software
  # row first, then loop left -> right -> left -> right -> left
  temp_1 <- well_fields %>% 
    separate_wider_position(well, 
                            c("row"=1, "column"=2), 
                            cols_remove = F, 
                            too_few = "align_start") %>% 
    mutate(column = as.integer(column)) %>% 
    arrange(row, column)
  
  temp_2 <- temp_1 %>% 
    distinct(row) %>% 
    mutate(well_row_odd = ifelse(row_number()%%2==0, 2, 1),
           well_row_id = row_number())
  
  temp_3 <- temp_1 %>% 
    distinct(column) %>% 
    mutate(well_column_odd = ifelse(row_number()%%2==0, 2, 1))
  
  temp_4 <- left_join(temp_1, temp_2, 
                      by = join_by(row)) %>% 
    left_join(temp_3, by = join_by(column)) %>% 
    mutate(well_row_order =
             ifelse(well_column_odd==1, well_row_id, -well_row_id),
           well_column_order =
             ifelse(well_row_odd==1, column, -column)
    )
  
  if(by == "lr") {
    # left->right, then right->left
    well_fields_optim <- temp_4 %>% 
      arrange(row, well_column_order, field)
  } else if(by == "ll") {
    # left->right, then left->right
    well_fields_optim <- temp_4 %>% 
      arrange(row, column, field)
  } else if(by == "tb") {
    # top->bottom, then bottom->top
    well_fields_optim <- temp_4 %>% 
      arrange(column, well_row_order, field)
  } else if(by == "tt") {
    # top->bottom, then top->bottom
    well_fields_optim <- temp_4 %>% 
      arrange(column, row, field)
  }
  
  # change all coords < 0 positions to 0
  if(clip_zero) {
    well_fields_optim <- well_fields_optim %>% 
      mutate(x=ifelse(x<0, 0, x), y=ifelse(y<0, 0, y))
  }

  
  return(well_fields_optim)
}



## plot func --------------------------------------------------------------
circleFun <- function(center = c(0, 0),
                      diameter = 1, 
                      npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(tibble(x = xx, y = yy))
}


plot_well_fields <- function(center, diameter,
                             well_fields, 
                             font_size=4) {
  ggplot(well_fields) + 
    geom_tile(aes(-x, -y, width = width, height = width),
              fill = "red3", color = NA) +
    geom_text(aes(-x, -y, label=field), size=font_size, color = "white") +
    geom_path(aes(-x, -y), linewidth = 1,  color = "grey50",
              data = circleFun(center, diameter)) +
    # geom_point(aes(x, y), 
    #            color="steelblue", size=4,
    #            data=well_fields
    # ) +
    coord_fixed() +
    theme_minimal(16) +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())
}




## simple shiny app to select cell -------------------------------------------

# Define a function to create the Shiny app
shiny_select_wells <- function(num_rows, num_cols) {
  # Load the necessary libraries
  library(shiny)
  library(DT)
  
  # Define the UI
  ui <- fluidPage(
    titlePanel("Interactive Well Selection"),
    mainPanel(
      DTOutput("wellTable"),
      verbatimTextOutput("selectedWells"),
      downloadButton("downloadData", "Download Selected Wells")
    )
  )
  
  # Define the server logic
  server <- function(input, output, session) {
    # Create a data frame to represent the well plate
    rows <- LETTERS[1:num_rows]
    cols <- as.character(1:num_cols)
    wellIDs <- outer(rows, cols, paste0)
    wellPlate <- as.data.frame(wellIDs)
    
    # Render the well table
    output$wellTable <- renderDT({
      datatable(wellPlate, 
                rownames = FALSE,
                colnames = rep("", num_cols),
                selection = "none", 
                extensions = c("Select", "Buttons"),
                options = list(
                  dom = "Blfrtip",
                  select = list(style = "os", items = "cell"),
                  buttons = c("selectCells", "selectAll", "selectNone"),
                  searching = F,
                  paging = F,
                  # pageLength = 16, 
                  ordering = F,
                  autoWidth = F,
                  info = F,
                  # fillContainer = F,
                  scrollX = T))
    }, server = FALSE)
    
    # Retrieve selected wells 
    selected_wells <- reactive({
      req(length(input$wellTable_cells_selected) > 0)
      mat <- as.matrix(input$wellTable_cells_selected)
      # fix column loc mis-location
      mat[, 2] <- mat[, 2] + 1
      wellPlate[mat]
    })
    
    # Display the selected wells
    output$selectedWells <- renderPrint({
      cat(paste(selected_wells(), collapse = ", "))
    })
    
    # Download the selected wells
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("selected_wells_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(as.matrix(selected_wells()), file, row.names = FALSE)
      }
    )
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}

