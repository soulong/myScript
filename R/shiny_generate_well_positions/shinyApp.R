
# set path
if(interactive()) {
  rstudioapi::getActiveDocumentContext()$path |>
    dirname() |>
    setwd()
}

# load packages
library(tidyverse)
library(shiny)
library(DT)


# load config
source("config.R")

# load funtions
source("functions.R")
source("server.R")
source("ui.R")

# app = shinyApp(ui, server)
app = shinyApp(ui, server, options = list(host='0.0.0.0', port=5556))


# # run from terminal
# cd ~/Documents/GitHub/myScript/generate_well_positions_for_nikon_microscopy

# # can close whole terminal
# # this should be for general use
# nohup R -e "shiny::runApp('shinyApp.R')" &

# # need to keep terminal activate, just free mouse in the termial
# R -e "shiny::runApp('shinyApp.R')"
# R -e "shiny::runApp('shinyApp.R')" > run.log 2>&1 &



