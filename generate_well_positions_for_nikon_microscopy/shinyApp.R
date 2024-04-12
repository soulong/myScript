
# set path
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |>
  setwd()

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


runApp(host="0.0.0.0", port=5556L, launch.browser = T, quiet = T)
# shinyApp(ui, server, 
#          options = list(host="0.0.0.0", port=5556L))





