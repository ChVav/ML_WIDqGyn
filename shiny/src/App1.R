library(shiny)
library(readxl)
library(stringr)
library(dplyr)

layout <- read.csv("./layout_template.csv") # template layout
exclude_targets <- c() #init exclude targets
all_assays <- list("assay1"=c("ZSCAN12","GYPC1"), 
                   "assay2"=c("DPP6","GYPC2"), 
                   "assay3"=c("GSX1","RALYL"),
                   "assay4"=c("COL2A1","ImC","EpC"))

# Define UI ----
ui <- fluidPage(
  # App title ----
  titlePanel(h3("Generate plate layout qWID")),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      checkboxGroupInput("checkGroup", #option to exclude targets
                         strong("Choose targets"), 
                         choices = list("ZSCAN12 & GYPC1" = "assay1", 
                                        "DPP6 & GYPC2" = "assay2", 
                                        "GSX1 & RALYL" = "assay3",
                                        "COL2A1 & ImC & EpC" = "assay4"),
                         selected = "assay4"),
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
  GeneratePlate <- eventReactive({input$file1
                                  input$checkGroup}, {
    # avoid displaying error message
    req(input$file1)
    
    #read in sample list
    sheet <- readxl::read_xlsx(input$file1$datapath,sheet = 1) 
    sheet <- sheet %>% mutate_all(as.character) # set all variables as character
    
    # warnings to add
    ### column "Number" does not exist
    ### column "Sample_name" does not exist
    ### no Samples entered
    ### more than 42 Samples entered
    
    #define assays to exclude
    if (length(input$checkGroup)<4){
      include_targets <- c() #init
      for (assay in input$checkGroup){
        include_targets <- c(include_targets,all_assays[[assay]])
      }
      exclude_targets <- setdiff(unname(unlist(all_assays)), include_targets)
      }
    
   # deal with <42 samples
    if (length(sheet$Number) < 42){
        x <- c()
        y <- c()
        for (i in 1:(42-length(sheet$Number))){
          x <- c(x, paste("Sample", formatC(length(sheet$Number)+i, width=2, flag="0"), sep="_"))
          y <- c(y, "empty")}
        df <- data.frame(x,y) %>% dplyr::rename(Number = x, Sample_name=y)
        sheet <- bind_rows(sheet, df)
      }
      
    # Paste in sample names
    df <- layout %>% slice(1:756) %>% select(Sample.Name)
    df <- left_join(df,sheet, by= c("Sample.Name"="Number"))
    layout$Sample.Name[1:756] <- c(df$Sample_name)
      
    # Some tidying, including emptying layout where samples/targets are missing
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
    layout <- layout %>%
        mutate(Sample.Name = ifelse(Target.Name %in% exclude_targets, "",Sample.Name),
               Sample.Color = ifelse(Target.Name %in% exclude_targets,"",Sample.Color),
               Target.Color = ifelse(Target.Name %in% exclude_targets, "",Target.Color),
               Task = ifelse(Target.Name %in% exclude_targets, "", Task),
               Reporter = ifelse(Target.Name %in% exclude_targets, "", Reporter),
               Quencher = ifelse(Target.Name %in% exclude_targets, "", Quencher)) %>%
        mutate(Target.Name = ifelse(Target.Name %in% exclude_targets, "",Target.Name)) 
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