
  
server <- function(input, output, session) {
  
  num_rows <- reactive(plate_type[[input$p_type]][["n_row"]])
  num_cols <- reactive(plate_type[[input$p_type]][["n_col"]])
  
  objectives <- reactive(str_subset(names(plate_type[[input$p_type]]), "\\d+x"))
  
  # ui for plate setting
  output$plate_info_ui <- renderUI({
    tagList(
      splitLayout(
        numericInput("a1_x", "A1 X", plate_type[[input$p_type]][["a1_x"]]),
        numericInput("a1_y", "A1 Y", plate_type[[input$p_type]][["a1_y"]]),
        cellWidths = c("50%", "50%") ),
      splitLayout(
        radioButtons("x_direct", "X direction", c(-1, 1), plate_type[[input$p_type]][["x_direct"]], inline=T),
        radioButtons("y_direct", "Y direction", c(-1, 1), plate_type[[input$p_type]][["y_direct"]], inline=T),
        cellWidths = c("50%", "50%") ),
      radioButtons("objective", "Objective", objectives(), "20x", inline = T),
      splitLayout(
        numericInput("pos_z", "Postion Z", plate_type[[input$p_type]][["pos_z"]]),
        numericInput("pos_pfs", "PFS offset", plate_type[[input$p_type]][["pos_pfs"]]),
        cellWidths = c("50%", "50%") )
    )
  })
  
  # generate well fields
  ref_well_fields <- reactive({
    req(!is.null(input$a1_x))
    generate_well_field_coords(
      plate_type[[input$p_type]][["well_diameter"]], 
      input$a1_x, 
      input$a1_y,
      plate_type[[input$p_type]][[input$objective]],
      plate_type[[input$p_type]][[input$objective]] * input$field_interval,
      input$n) %>%
      optim_fields()
  })
  # observe({ref_well_fields()})
  
  # ui for subset field
  output$subset_field_ui <- renderUI({
    selectizeInput("subset_field", "Subset Field",
                   choices = unique(ref_well_fields()[["field"]]),
                   selected = unique(ref_well_fields()[["field"]]),
                   multiple = T)
  })
  
  # filter fields
  ref_well_fields_subset <- reactive({
    req(!is.null(input$subset_field))
    ref_well_fields() %>% 
      filter(field %in% input$subset_field)
  })

  # check fields
  output$plot_fields <- renderPlot({
    req(nrow(ref_well_fields_subset()) > 0)
    plot_well_fields(
      c(input$a1_x, input$a1_y),
      plate_type[[input$p_type]][["well_diameter"]],
      ref_well_fields_subset())
  })

  # shift position
  all_well_fields <- reactive({
    req(nrow(ref_well_fields_subset()) > 0)
    shift_well_coords(
      num_rows(),
      num_cols(),
      plate_type[[input$p_type]][["well_interval_x"]],
      plate_type[[input$p_type]][["well_interval_y"]],
      x_direct = as.integer(input$x_direct),
      y_direct = as.integer(input$y_direct)) %>%
      # glimpse()
      crossing(ref_well_fields_subset()) %>%
      transmute(well, field,
                x = x + x_shift,
                y = y + y_shift)
  })
  
  # select wells
  # Create a data frame to represent the well plate
  wellPlate <- reactive({
    outer(LETTERS[1:num_rows()],
          as.character(1:num_cols()),
          paste0) %>%
      as.data.frame()
    })
  # Render the well table
  output$selectWells <- renderDT({
    datatable(wellPlate(),
              rownames = FALSE,
              colnames = rep("", ncol(wellPlate())),
              # selection = list(target = "cell"),
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
    req(length(input$selectWells_cells_selected) > 0)
    mat <- as.matrix(input$selectWells_cells_selected)
    # fix column loc mis-location
    mat[, 2] <- mat[, 2] + 1
    wellPlate()[mat]
  })
  
  output$display_selected_wells <- renderPrint({
    selected_wells()
  })

  # filter wells
  well_fields_filtered <- reactive({
    all_well_fields() %>%
      filter(well %in% selected_wells())
  })
  # observe({print(well_fields_filtered())})
  
  
  # optimize well coords
  well_fields_filtered_optim <- reactive({
    optim_wells(well_fields_filtered(), 
                by = input$image_order,
                clip_zero = ifelse(input$clip_zero=="Yes", TRUE, FALSE))
  })

  # add z and pfs info
  final <- reactive({
    well_fields_filtered_optim() %>%
      mutate(PositionID = row_number() %>%
               str_pad(5, "left", "0") %>% str_c("Point", .),
             bChecked = "true",
             strName = str_c(well, "#", field),
             dXPosition = round(x,1),
             dYPosition = round(y,1),
             dZPosition = input$pos_z,
             dPFSOffset = input$pos_pfs,
             .before = 1)
  })

  # Update view file
  output$preview <- renderDT({
    datatable(final(),
              options = list(
                # dom = '<"top" p>',
                searching = T,
                paging = T,
                pageLength = 5,
                # ordering = F,
                autoWidth = T,
                info = T,
                # fillContainer = FALSE,
                scrollX=T
              )) })
  
  # Download handler
  output$download <- downloadHandler(
    filename = function() {
      str_glue("{Sys.Date()}_{input$fname}.xml")
    },
    content = function(file) {
      
      withCallingHandlers({
        shinyjs::html("pb", "")
        save_xml(final(), file)
      },
      message = function(m) {
        shinyjs::html(id = "pb", html = m$message, add = FALSE)
      })
      
      # res <- create_empty_xml(fname)
      # fname <- res$fname
      # doc <- res$doc
      # 
      # 
    }
  )
  

}
