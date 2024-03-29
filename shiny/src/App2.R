
# input is fixed to uploading one results file exported from Quantstudio 7, .csv format
# ! will have to check naming of the standards and controls!!!
# output files are in a .zip
# pic of the standard curve
# Excel workbook with all processed info (PMR for all targets, Ct+SD, DNA input, analysis reps)
# When WID-qEC targets are on the plate, the final csv file with the results for the commercialized test is included

# load libraries
library(shiny)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggpubr)

# Determine pipeline version
version <- read.table("../../version.csv", header=FALSE) %>% pull(V1)
version <- paste0("Shiny ",version)

# Defined Settings
threshold_COL2A1 <- 35
threshold_targets <- 38
external_curve <- NULL
fix_intercept <- 36.9 #these values come from meta-analysis, but will not be used
fix_slope <- -3.4
calib_fixed <- FALSE
inconclusive <- FALSE #when true, inconclusive reactions will be removed, not yet implemented
inconsistent <- FALSE #when true, only CTs for targets for which only one rep amplified will be set to 0

# source main wrapper function
source("./wrapPMR.R")

# Define UI ----
ui <- fluidPage(
  # App title ----
  titlePanel(h3("Calculate PMR for WID-qGYN")),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      fileInput("file1", "Choose Results Analysis File", # Input: Select a file ----
                multiple = FALSE,
                accept = c(".csv")),
      tableOutput("resultsQuantstudio"), # for display uploaded data
      tags$hr(), # Horizontal line ----
      downloadButton("download", "Download PMR Result File") #Download button
    ),

    # Main panel for displaying part of output ----
    mainPanel(
      imageOutput("stdCurve") # display standard curve
    ))
)

# Define server logic ----
server <- function(input, output) {
  #Calculate PMR and standard curve upon uploading of files exported from QuantStudio
  myResults <- eventReactive({input$file1}, {
    # avoid displaying error message
    req(input$file1)
    
    #read in inputfiles
    data <- read.table(file = input$file1$datapath,
                        sep = ",",
                        header = TRUE)
    
    # ensure sample IDs are in correct format
    data$Sample <- as.character(data$Sample)
   
    # calculate PMR
    results <- calculate_pmr(data, threshold_COL2A1, threshold_targets, external_curve, fix_intercept, fix_slope, calib_fixed, inconclusive, inconsistent, version)
    results
  })
  
  # display uploaded data
  output$resultsQuantstudio <- renderTable({
    req(input$file1)
    data <- read.table(file = input$file1$datapath,
                       sep = ",",
                       header = TRUE)
    head(data)[,1:6]
  })
  
  # display calculated curve
  output$stdCurve <- renderPlot({
    myResults()[[1]]
  })
  
  # download batch results
  output$download <- downloadHandler(
    
    filename = function() {paste0("Final_results_", Sys.Date(), ".zip")},
    content = function(file){
      
      #locate temporary dir created in this R session & create sub-dir with time within
      temp_dir <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_dir)
      
      # save cal curve to temporary dir
      ggsave(myResults()[[1]],
             file = paste0(temp_dir, "/","calibration_curve.png"),
             width = 5,
             height = 4) 
      
      # save CT plot to temporary dir
      ggsave(myResults()[[10]],
             file = paste0(temp_dir, "/","CTCOL2A1_MeanStdev.png"),
             width = 85,
             height = 80,
             unit = "mm") 
      
      # save batched results to temporary dir
      wb <- createWorkbook()
      
      addWorksheet(wb, "PMR Values")
      writeDataTable(wb = wb, sheet = 1, x = myResults()[[2]], rowNames=FALSE) #PMRs
      addWorksheet(wb, "Cq mean")
      writeDataTable(wb = wb, sheet = 2, x = myResults()[[3]], rowNames=FALSE) #mean Cq
      addWorksheet(wb, "Cq SD")
      writeDataTable(wb = wb, sheet = 3, x = myResults()[[4]], rowNames=FALSE) #Cq stdev
      addWorksheet(wb, "Input conc")
      writeDataTable(wb=wb, sheet = 4, x = myResults()[[5]], rowNames=FALSE) #log(copy number/5uL)
      addWorksheet(wb, "low DNA input")
      writeDataTable(wb=wb, sheet = 5, x = myResults()[[6]], rowNames=FALSE) #samples for which COL2A1 failed in both reps, should have negative controls
      addWorksheet(wb, "Reprocessing needed")
      writeDataTable(wb = wb, sheet = 6, x = myResults()[[7]], rowNames=FALSE) # samples for which for only one of two reps COL2A1 failed
      addWorksheet(wb, "Reprocessing recommended")
      writeDataTable(wb = wb, sheet =7, x = myResults()[[8]], rowNames=FALSE) # samples for which for only one of two reps target amplified
      addWorksheet(wb, "Warning COL2A1 SD high")
      writeDataTable(wb=wb, sheet = 8, x= myResults()[[9]], rowNames=FALSE)# samples for which SD CT COL2A1 > 1.5
      addWorksheet(wb,"Log")
      writeDataTable(wb=wb, sheet = 9, x=myResults()[[11]], rowNames=FALSE, colNames=FALSE) # save pipeline version and settings

      saveWorkbook(wb, paste0(temp_dir, "/","batch_results.xlsx"))
      
      # When WID-qEC targets are on the plate, save final result summary separately to temporary dir
      if(length(myResults())>10){
        write.csv(myResults()[[11]], 
                  file=paste0(temp_dir, "/", "final_results_commercial_WIDqEC.csv"), 
                  row.names=FALSE, 
                  fileEncoding = "UTF-8")
      }
      
      # zip up temporary dir to file
      zip::zip(
        zipfile = file,
        files = dir(temp_dir),
        root = temp_dir
      )
    },
    contentType = "application/zip"
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)