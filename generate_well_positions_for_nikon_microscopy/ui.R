
ui <- fluidPage(
  shinyjs::useShinyjs(),
  
  # titlePanel(div("Optimized Well Fields Generation for ShuLab Nikon Microscopy", 
  #                style = "font-size: 70%; width: 100%")),
  tags$head(tags$title("🎭 ShuLab")),
  h3("An Optimized Well Fields Generation App for General Microscopy"),
  tags$head(tags$style('h3 {color:#386cb0;}')),
  
  sidebarLayout(
    
    sidebarPanel(
      div(
        selectInput("p_type", "Plate Type", names(plate_type), names(plate_type)[1]),
        uiOutput("plate_info_ui"),
        # numericInput("a1_x", "A1 X", 0),
        # numericInput("a1_y", "A1 Y", 0),
        numericInput("pos_z", "Postion Z", 0),
        numericInput("pos_pfs", "PFS offset", 0),
        radioButtons("objective", "Objective", c("100x","60x","20x"), "100x", inline = T),
        radioButtons("image_order", "Image Order", c("row","column"), "row", inline = T),
        numericInput("field_interval", "Interval between field", 1),
        sliderInput("n", "Field [NxN]", 1, 10, 4),
        textInput("fname", "File Name Suffix"),
        downloadButton("download", "Download XML"),
        
        style = "font-size: 85%; width: 100%"
      ),
      hr(),
      em("For bug report, refer to"),
      strong("haohe90@gmail.com"),
      width = 2
    ),
    
    mainPanel(
      column(8,
             div(DTOutput("selectWells"), style = "font-size: 70%; width: 100%"),
             hr(),
             # br(),
             div(DTOutput("preview"), style = "font-size: 80%; width: 100%")
      ),
      
      column(4,
             plotOutput("plot_fields", width="400px"),
             br(),
             verbatimTextOutput("display_selected_wells"),
             br(),
             textOutput("pb")
      ),
      # div(DTOutput("select_wells"), style = "font-size: 75%; width: 80%"),
      # DTOutput("preview")
      # column(6,
      #        plotOutput("plot_fields")
      #        ),
      # column(6,
      #        DTOutput("preview", width="600px")
      #        )
      # Show downloaded file name
      # verbatimTextOutput("selected_wells"),
      width = 10
    )
    
  )
)
