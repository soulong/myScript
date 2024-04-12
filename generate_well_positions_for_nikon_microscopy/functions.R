
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
    source="2023-11-30_template_shulab_empty.xml"
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




# [DEPRECATED] random points 

# # gererate random field within radius
# generate_rand_coords <- function(center_x, center_y, 
#                                  radius, N=5, 
#                                  well_radius=NULL) {
#   
#   # r <- radius * sqrt(runif(N, 0, 1))
#   # theta <- seq(0, 2 * pi, length.out = N + 1)[-1]
#   # theta <- seq(0, 2 * pi - 2 * pi / N, length.out = N)
#   # x_coords <- center_x + r * cos(theta)
#   # y_coords <- center_y + r * sin(theta)
#   
#   # Define the interval between the points
#   I <- 2 * pi / N
#   
#   # Generate the coordinates of the points
#   theta <- seq(0, 2 * pi - I, length.out = N)
#   r <- sqrt(runif(N, 0, 1)) * radius
#   x_coords <- center_x + r * cos(theta)
#   y_coords <- center_y + r * sin(theta)
#   
#   
#   # Check if the generated points are within the circle
#   if(!is.null(well_radius)) {
#     distances <- sqrt((x_coords - center_x) ^ 2 + (y_coords - center_y) ^ 2)
#     while (any(distances > well_radius)) {
#       r[which.max(distances)] <- radius * sqrt(runif(1, 0, 1))
#       x_coords <- center_x + r * cos(theta)
#       y_coords <- center_y + r * sin(theta)
#       distances <- sqrt((x_coords - center_x) ^ 2 + (y_coords - center_y) ^ 2)
#     }
#   }
#   
#   # Combine the x and y coordinates into a matrix
#   coords <- cbind(x_coords, y_coords)
#   rownames(coords) <- seq_len(nrow(coords))
#   
#   coords <- as_tibble(coords, rownames="field") %>% 
#     arrange(x_coords, y_coords)
#   
#   return(coords)
# }
# 
# 
# optim_points <- function(new_points, pos_z, pos_pfs) {
#   ## arrange positions
#   # to solve auto order problem when loading to software
#   # row first, then loop left -> right -> left -> right -> left
#   temp_1 <- new_points %>% 
#     separate_wider_position(well, 
#                             c("row"=1, "column"=2), 
#                             cols_remove = F, 
#                             too_few = "align_start") %>% 
#     mutate(column = as.integer(column)) %>% 
#     arrange(row, column)
#   
#   temp_2 <- temp_1 %>% 
#     distinct(row) %>% 
#     mutate(row_odd=ifelse(row_number()%%2==0, 2, 1))
#   
#   temp_3 <- left_join(temp_1, temp_2, 
#                       by = join_by(row)) %>% 
#     mutate(column_odd_reorder =
#              ifelse(row_odd==1, column, n_col - column),
#            point_x_odd_reorder=
#              ifelse(row_odd==1, 100000 - x, x)
#     )
#   
#   new_points_optim <- temp_3 %>% 
#     arrange(row, column_odd_reorder, point_x_odd_reorder, y) %>% 
#     mutate(point_x_odd_reorder = row_number(), .by=well) %>% 
#     mutate(PositionID=row_number() %>% 
#                 str_pad(5, "left", "0") %>% str_c("Point", .),
#               bChecked="true",
#               strName=str_c(well, "#", field),
#               # strName=str_c(row, LETTERS[column_odd_reorder],
#               #               "_", well, "#", field),
#               # strName=str_c(row, LETTERS[column_odd_reorder],
#               #               "_", well, "#", point_x_odd_reorder),
#               dXPosition=round(x,1),
#               dYPosition=round(y,1),
#               dZPosition=pos_z,
#               dPFSOffset=pos_pfs,
#               .before = 1
#     )
#   
#   return(new_points_optim)
# }
# 
# 
# 
# 
# plot_well_points <- function(center, 
#                              diameter,
#                              points,
#                              point_size=4
# ) {
#   circleFun(center, diameter) %>% 
#     ggplot(aes(x, y)) + 
#     geom_path() +
#     geom_point(aes(x_coords, y_coords), 
#                color="steelblue", size=4,
#                data=points
#     ) +
#     coord_fixed() +
#     labs(x="", y="") +
#     theme_bw()
# }


## generate NxN field --------------------------------------------------------------

shift_well_coords <- function(n_row=16, n_col=24, 
                              well_interval_x=50, 
                              well_interval_y=50, 
                              xtype=-1, ytype=1) {
  # Loop over the wells
  shift <- matrix(nrow=n_row*n_col, ncol=2)
  colnames(shift) <- c("x_shift", "y_shift")
  rownames(shift) <- seq_len(nrow(shift))
  
  i <- 1
  for (row in seq_len(n_row)) {
    for (col in seq_len(n_col)) {
      shift[i, ] <- c(xtype *(col - 1) * well_interval_x, 
                      ytype *(row - 1) * well_interval_y)
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
optim_wells <- function(well_fields, by="row") {
  
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
  
  # read row first then column
  if(by == "row") {
    well_fields_optim <- temp_4 %>% 
      arrange(row, well_column_order, field)
  } else {
    # by = "column"
    # read column first then row
    well_fields_optim <- temp_4 %>% 
      arrange(column, well_row_order, field)
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
    geom_tile(aes(x, y, width = width, height = width),
              fill = "red3", color = NA) +
    geom_text(aes(x, y, label=field), size=font_size, color = "white") +
    geom_path(aes(x, y), linewidth = 1,  color = "grey50",
              data = circleFun(center, diameter)) +
    # geom_point(aes(x, y), 
    #            color="steelblue", size=4,
    #            data=well_fields
    # ) +
    coord_fixed() +
    labs(x="", y="") +
    theme_minimal(16)
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
    
    # Display the selected wells
    output$selectedWells <- renderPrint({
      req(length(input$wellTable_cells_selected) > 0)
      mat <- as.matrix(input$wellTable_cells_selected)
      cat("Selected wells:\n")
      mat[, 2] <- mat[, 2] + 1
      cat(paste(wellPlate[mat], collapse = ", "))
    })
    
    # Download the selected wells
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("selected_wells_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        selected <- input$wellTable_cells_selected
        if (length(selected) > 0) {
          selectedWells <- wellPlate[as.matrix(selected)]
          write.csv(selectedWells, file, row.names = FALSE)
        }
      }
    )
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}

