library(shiny)
library(readxl)
library(stringr)
library(dplyr)

layout <- read.csv("./layout_template_c_CIN.csv") # template layout
layout <- layout %>% arrange(SORT.A)

# Define UI ----
ui <- fluidPage(
  # App title ----
  titlePanel(h3("Generate plate layout qWID-CIN")),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      fileInput("file1", "Choose Excel File", # Input: Select a file ----
                multiple = FALSE,
                accept = c(".xlsx")),
      tableOutput("list"),
      tags$hr(), # Horizontal line ----
      downloadButton("download", "Download plate layout") #Download button
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      tableOutput("layout") #output generated plate layout
    ))
)

# Define server logic ----
server <- function(input, output) {
  #Generate plate layout upon upload sample list
  GeneratePlate <- eventReactive({input$file1}, {
    # avoid displaying error message
    req(input$file1)
    
    #read in sample list
    sheet <- readxl::read_xlsx(input$file1$datapath,sheet = 1) 
    sheet <- sheet %>% mutate_all(as.character) # set all variables as character
    
    # warnings to add
    ### column "Number" does not exist
    ### column "Sample_name" does not exist
    ### no Samples entered
    ### more than 90 Samples entered
    
   # deal with <90 samples
    if (length(sheet$Number) < 90){
        x <- c()
        y <- c()
        for (i in 1:(90-length(sheet$Number))){
          x <- c(x, paste("Sample", formatC(length(sheet$Number)+i, width=2, flag="0"), sep="_"))
          y <- c(y, "empty")}
        df <- data.frame(x,y) %>% dplyr::rename(Number = x, Sample_name=y)
        sheet <- bind_rows(sheet, df)
      }
      
    # Paste in sample names
    df <- layout %>% slice(1:720) %>% select(Sample.Name)
    df <- left_join(df,sheet, by= c("Sample.Name"="Number"))
    layout$Sample.Name[1:720] <- c(df$Sample_name)
      
    # Some tidying, including emptying layout where samples/targets are missing
    layout <- layout %>% arrange(SORT.B)
    layout <- layout %>% select(-c(SORT.A,SORT.B))
    layout$Biogroup.Color <- ""
    layout$Biogroup.Name <- ""
    layout$Comments <- ""
    layout$Quantity <- ""
    layout <- layout %>% 
        mutate(Sample.Color = ifelse(Sample.Name == "empty","",Sample.Color),
               Target.Name = ifelse(Sample.Name == "empty", "",Target.Name),
               Target.Color = ifelse(Sample.Name == "empty", "",Target.Color),
               Task = ifelse(Sample.Name == "empty", "", Task),
               Reporter = ifelse(Sample.Name == "empty", "", Reporter),
               Quencher = ifelse(Sample.Name == "empty", "", Quencher))
    layout["Sample.Name"][layout["Sample.Name"] == "empty"] <- ""
    layout <- layout %>% distinct() # remove duplicated rows
    layout$Sample.Name <- as.character(layout$Sample.Name) #ensure these are characters
    colnames(layout) <- stringr::str_replace(colnames(layout), "[.]", " ")
    layout
  })
  
  # display uploaded data
  output$list <- renderTable({
    req(input$file1)
    sheet <- readxl::read_xlsx(input$file1$datapath,
                               sheet = 1)
    sheet
  })
  
  # display calculated layout
  output$layout <- renderTable({ GeneratePlate() })
  
  #generate download file
  output$download <- downloadHandler(
    
    filename = function() {"plate_layout.csv"},
    content = function(file){
      #add header, required for correct normalization with passive reference dye, and save file
      writeLines(c("* Block Type = 384-Well Block",
                   paste0("* Date Created = ", Sys.Date()),
                   "* Passive Reference = ROX",
                   "* Barcode = ",
                   ""), 
                 file)
      write.table(GeneratePlate(),  
                  file,
                  append=T, 
                  sep=',',
                  qmethod="double",
                  row.names=F, 
                  col.names=T )
    }
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)