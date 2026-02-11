
# app.R
library(shiny)
library(dplyr)

# Function: Rank-based Inverse Normal Transformation (INT)
# Preserves original order
rank_normal_transform <- function(x) {
  n <- length(x)
  if (n == 0) return(numeric(0))
  
  # Handle ties properly using average ranks
  ranks <- rank(x, ties.method = "average")
  
  # Transform ranks to normal quantiles (van der Waerden method)
  # Using qnorm((ranks - 0.5)/n) is a common choice
  # Alternatives: (ranks)/(n+1), (ranks - 0.375)/(n + 0.25) [Blom], etc.
  transformed <- qnorm((ranks - 0.5) / n)
  
  # Return in original order
  transformed[order(match(seq_along(x), order(ranks)))] <- transformed
  return(transformed)
}

ui <- fluidPage(
  titlePanel("Rank-Based Normal Transformation (Preserves Order)"),
  
  sidebarLayout(
    sidebarPanel(
      tags$h4("Input Data"),
      tags$p("Enter numbers separated by spaces, commas, tabs, or newlines:"),
      tags$textarea(id = "input_numbers", rows = 10, cols = 40,
                    placeholder = "1.2  5.4  3.1  2.8  10.5\n-2  0  1.5"),
      
      br(), br(),
      actionButton("transform_btn", "Transform to Normal", class = "btn-primary"),
      
      br(), br(),
      downloadButton("download_data", "Download Transformed Data (.csv)")
    ),
    
    mainPanel(
      h3("Original vs Transformed Data"),
      tableOutput("comparison_table"),
      
      br(),
      h4("Summary Statistics"),
      verbatimTextOutput("summary_orig"),
      verbatimTextOutput("summary_trans"),
      
      br(),
      plotOutput("density_plot")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive: parse input numbers
  input_data <- eventReactive(input$transform_btn, {
    text <- input$input_numbers
    if (trimws(text) == "") return(numeric(0))
    
    # Parse numbers flexibly
    nums <- as.numeric(unlist(strsplit(text, "[,\\s\\n\\t]+")))
    nums <- nums[!is.na(nums)]
    return(nums)
  })
  
  # Reactive: transformed data
  transformed_data <- reactive({
    x <- input_data()
    if (length(x) == 0) return(numeric(0))
    rank_normal_transform(x)
  })
  
  # Comparison table
  output$comparison_table <- renderTable({
    x <- input_data()
    y <- transformed_data()
    if (length(x) == 0) return(data.frame())
    
    data.frame(
      Index = 1:length(x),
      Original = x,
      Transformed = round(y, 4)
    )
  }, digits = 6)
  
  # Summary stats
  output$summary_orig <- renderPrint({
    summary(input_data())
  })
  
  output$summary_trans <- renderPrint({
    summary(transformed_data())
  })
  
  # Density plot
  output$density_plot <- renderPlot({
    x <- input_data()
    y <- transformed_data()
    if (length(x) < 2) return(NULL)
    
    par(mfrow = c(1, 2))
    hist(x, main = "Original Data", xlab = "Value", col = "lightblue", border = "white")
    hist(y, main = "Transformed (Normal)", xlab = "Value", col = "lightcoral", border = "white")
    
    # Overlay normal curve
    curve(dnorm(x, mean(y), sd(y)), add = TRUE, col = "red", lwd = 2)
  })
  
  # Download transformed data
  output$download_data <- downloadHandler(
    filename = function() {
      paste("normal_transformed_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      df <- data.frame(
        Original = input_data(),
        Transformed = transformed_data()
      )
      write.csv(df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)