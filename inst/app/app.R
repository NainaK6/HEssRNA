# A shiny app for calculating the sample size based on heritability along with other parameters such as Power, FoldChange and Tissue



# source("sourceFunc.R")
library(shiny)
library(HEssRNA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(DESeq2)
library(ssizeRNA)
library(DT)
library(shinythemes)
library(shinyjs)  # Load shinyjs for overlay functionality

js <- "
function openFullscreen(elem) {
  if (elem.requestFullscreen) {
    elem.requestFullscreen();
  } else if (elem.mozRequestFullScreen) { /* Firefox */
    elem.mozRequestFullScreen();
  } else if (elem.webkitRequestFullscreen) { /* Chrome, Safari and Opera */
    elem.webkitRequestFullscreen();
  } else if (elem.msRequestFullscreen) { /* IE/Edge */
    elem.msRequestFullscreen();
  }
}"





# UI for the application
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  #titlePanel("HEss-Shiny: Heritability-based Estimation of Sample Size for RNA-seq"),
  tags$div(
    style = "font-size: 2em; font-weight: bold; margin-top: 20px; background-image: url('bg_blur1.png'); background-size: cover; background-position: center; display: flex; align-items: center;",
    tags$div(
      style = "width: 15%;",
      tags$img(src = "icar.png", width = "60%", alt = "ICAR Logo", style = "margin-right: 10px; margin-left:10px;"),
      
    ),
    tags$div(
      style = "width: 85%; text-align: center;",
      span("HEssRNA-Shiny", style = "color: #8B008B; font-size: 1.5em;"),
      tags$br(),
      span("Heritability-based Estimation of Sample Size for RNA-Seq", style = "color: #000000;")
    ),
    tags$div(
      style = "width: 20%;",
      tags$img(src = "logo1.jpg", width = "60%", style="margin-top:10px; margin-left:70px;")
    ),
    
  ),
  
  tags$head(
    tags$style(
      type="text/css", "#inline label{ display: table-cell; text-align: left; vertical-align: middle; }#inline .form-group { display: table-row;}",
      HTML(".overlay {
        position: fixed;
        top: 50px;
        right: 10px;
        width: 35%;
        height: auto;
        background-color: #eafaf1;
        border: 1px solid #ccc;
        padding: 20px;
        box-shadow: -5px 0 15px rgba(0,0,0,0.3);
        display: none;  /* Initially hidden */
        z-index: 1000;
        overflow-y: auto;
      }
      .closebtn {
        font-size: 18px;
        font-weight: bold;
        cursor: pointer;
        color: #888;
        margin-bottom: 15px;
        text-align: right;
      }
      .closebtn:hover {
        color: #ff0000;
        background-color: #45a049; /* A darker shade for hover effect */
      }
      #countSmplDat {
        height: auto;           /* Allow the height to adjust automatically */
        margin-bottom: 5px;    /* Optional: Add some space below the file input */
        font-size: 15px;
        height: 20px;  /* Adjust height of the entire file input area */
      }
      #countSmplDat .btn-file {
        font-size: 15px;    /* Font size of the 'Browse' button */
        padding: 6px 8px;   /* Adjust padding */
        height: 40px;       /* Adjust height of 'Browse' button */
        background-color: #4CAF50;
        color: white;
      }
      #countSmplDat .form-control-file {
        font-size: 13px;    /* Font size of the file name text */
        height: 15px;       /* Adjust height of the file name text */
      }
      #countSmplDat .shiny-file-input {
        height: 25px;           /* Height of the file input area */
        line-height: 25px;      /* Align text vertically */
        "
      )
    )
  ),
  navbarPage(
    "",
    # Changing themes could be managed here
    theme = shinytheme("readable"), fluid = T,
    tabPanel(tags$b("Home") , icon = icon("home"),
             uiOutput("outHome")
    ), 
    
    tabPanel(tags$b("Tool"), 
             sidebarLayout(
               sidebarPanel(
                 selectInput("heritabilityClass", 
                             "Heritability Class:", 
                             choices = list("High" = "high", "Mid" = "mid", "Low" = "low")),
                 
                 actionButton("powerCalcBtn", "Power Calculation", class = "btn-primary", style="float:right;"), # Button for power calculation
                 
                 numericInput("inptPwr", 
                              "Power (0 to 1):", 
                              value = 0.8, 
                              min = 0, 
                              max = 1, 
                              step = 0.01),
                 
                 numericInput("fc", 
                              "Fold Change:", 
                              value = 1.5, 
                              min = 0, 
                              step = 0.1),
                 
                 textInput("trait", 
                           "Trait (optional):", 
                           value = ""),
                 
                 textInput("tissue", 
                           "Tissue (optional):", 
                           value = ""),
                 
                 fileInput("df4modelFile", 
                           "Upload DataFrame for model generation (optional):", 
                           accept = c(".csv")),
                 
                 fileInput("hIndexMeanDFFile", 
                           "Upload DataFrame for Heritability Index Mean (optional):", 
                           accept = c(".csv")),
                 
                 actionButton("submit", "Submit")
               ),
               
               mainPanel(
                 # h3("Estimated Sample Size:"),
                 # verbatimTextOutput("sampleSize")
                 uiOutput("outPanel")
                 #uiOutput("outPanelSummary"),
                 #uiOutput("outPanelPowerCalc")  # For power calculation results
               )
             )
    ),
    
    tabPanel(tags$b("Help"),
             uiOutput("outDocmnt")
    ),
    tabPanel(tags$b("Team"),
             uiOutput("outTeam")
    ),
  ),
  #tags$br(),
  # Copyright Information
 # column(12, 
        # "copyright @ICAR_IASRI",
         
         
  #),
 # Footer with copyright, logos, and map widget
 tags$footer(
   style = "width: 100%; background-color: lightblue; padding: 5px; z-index: 1000;",
   fluidRow(
     column(
       width = 4,
       tags$div(
         style = "text-align: left; padding-left: 10px; margin-top: 20px",
         
         tags$img(src = "iari.jpg", height = "90px", alt = "IARI Logo", style = "margin-right: 10px;"),
         tags$img(src = "iasri.jpeg", height = "90px", alt = "IASRI Logo", style = "margin-right: 10px;"),
       ),
     ),
     column(
       width = 5,
            tags$div(
              style = "text-align: center; padding-top:20px; padding-right: 10px",
              tags$span(style = "font-size: 16px;", "© 2024 ICAR-Indian Agricultural Statistics Research Institute, Govt. of India", "Library Avenue, PUSA, New Delhi - 110 012 (INDIA)", "All rights reserved")
            )
       ),
     column(
       width = 3,
       tags$div(
         style = "text-align: right; padding-right: 10px; width: 200px; height: 200px;",
         HTML('<script type="text/javascript" src="//rf.revolvermaps.com/0/0/8.js?i=5mtsplw9n07&amp;m=0&amp;c=ff0000&amp;cr1=ffffff&amp;f=arial&amp;l=33" async="async"></script>')
           )
         )
       )
     ),
  
  # Overlay UI for Power Calculation Tool
  div(class = "overlay", id = "powerCalcOverlay",
      div(style = "display: flex; justify-content: space-between; align-items: center; 
             background-color: #4CAF50; color: white; padding: 10px; border-radius: 5px;",
          h4("Power Calculation Options", style = "margin: 0;"),  # Title
          actionButton("closeOverlay", "×", class = "btn btn-light", style = "font-size: 18px;")
      ),
      hr(),
      
      div(style = "height: 400px; overflow-y: auto;",
      fluidRow(id = "countSmplDat",
        column(6,
               fileInput("countDat", tags$text("Upload Count Data (.csv)"), accept = c(".csv")),
               ),
        column(6,
               fileInput("smplDat", "Upload Sample MetaData (.csv)", accept = c(".csv")),
               ),
        br(),
        column(3,
               numericInput("alpha_power", "Alpha (significance level):", value = 0.05, min = 0, max = 1, step = 0.01),
               ),
        column(3,
               numericInput("thrsholdFC", "Log Fold Change Threshold:", value = 2, min = 0, step = 0.1),
               ),
        column(3,
               numericInput("sims", "No of Simulations (sims):", value = 50, min = 1, step = 1),
               ),
        column(3,
               numericInput("inptNoOfReplicates", "No of Replicates (Sample size):", value = 3, min = 0, step = 1),
        ),
        br(),
        column(8,
               # checkboxInput("saveFltrd", "Save Filtered Data?", value = FALSE), style="font-size:15px;",
               ),
        column(4,
               actionButton("submitPowerCalc", "Run Power Calculation", class = "btn-success", style="float:right;")
               )
      ),
      hr(),
      uiOutput("outPanelPowerCalc"),
      )
  )
  
  
)

# Server logic for the application
server <- function(input, output) {
  # Home page to be included while showing on the CoreDECAP tab
  output$outHome <- renderUI({
    includeHTML("www/home2.html")
    #tags$b("Home")
  })
  # OR use html file itself for Documentation tab
  output$outDocmnt <- renderUI({
    includeHTML("www/help.html")
    #tags$b("Help")
  })
  # use html file for Contact tab
  output$outTeam <- renderUI({
    includeHTML("www/team.html")
    #tags$b("Team")
  })
  
  # Show the power calculation overlay when the button is clicked
  observeEvent(input$powerCalcBtn, {
    shinyjs::show("powerCalcOverlay")  # Show overlay
  })
  
  # Close the overlay
  observeEvent(input$closeOverlay, {
    shinyjs::hide("powerCalcOverlay")  # Hide overlay
  })
  
  # Placeholder for df4model and hIndexMeanDF
  df4modelInpt <- reactive({
    if (is.null(input$df4modelFile)) {
      # Use default dataframe
      read.csv("Examples/df4modelInpt.csv", header = TRUE, stringsAsFactors = FALSE)
    } else {
      read.csv(input$df4modelFile$datapath)
    }
  })
  
  hIndexMeanDF <- reactive({
    if (is.null(input$hIndexMeanDFFile)) {
      # Use default dataframe
      read.csv("Examples/hIndexMeanDF.csv", header = TRUE, stringsAsFactors = FALSE)
    } else {
      read.csv(input$hIndexMeanDFFile$datapath)
    }
  })
  
  # Creating model based on the sample size based on inputs
  pred.model <- reactive({
    input$submit  # Wait for user to click the submit button
    
    isolate({
      # Retrieve input values
      heritabilityClass <- input$heritabilityClass
      inptPwr <- input$inptPwr
      fc <- input$fc
      trait <- if (input$trait != "") input$trait else NULL
      tissue <- if (input$tissue != "") input$tissue else NULL
      
      # Call the smplSizPred function
      smplSizPredModel(
        df4model = df4modelInpt(),
        heritabilityClass = heritabilityClass,
        inptPwr = inptPwr,
        fc = fc,
        trait = trait,
        tissue = tissue
      )
    })
  })
  
  # output$pred.model <- renderText({
  #   req(pred.model())  # Ensure pred.model is there
  #   paste(pred.model())
  # })
  
  # Render the predicted model info
  predModelInfo <- reactive({
    predModelCoeff <- pred.model()$coefficients
    predModelCoeff <- as.data.frame(predModelCoeff)
    predModelInfo <- data.frame(Parameter = rownames(predModelCoeff), Value = predModelCoeff[,1])
    return(predModelInfo)
    
  })
  output$pred.model.coeff <- renderDT(
    #as.data.frame(pred.model()$coefficients),
    predModelInfo(),
    options = list(dom = 't',
                   initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': '#990F4B', 'color': 'white'});",
                     "}")
    )
  )
  
  # Plot Model
  output$modelCRplot <- renderPlot({
    crPlots(pred.model())
  })
  
  output$modelPlot <- renderPlot({
    par(mfrow = c(2, 2))
    plot(pred.model())
  })
  
  # Calculate the sample size based on inputs
  sampleSize <- reactive({
    input$submit  # Wait for user to click the submit button
    
    isolate({
      # Retrieve input values
      heritabilityClass <- input$heritabilityClass
      inptPwr <- input$inptPwr
      fc <- input$fc
      trait <- if (input$trait != "") input$trait else NULL
      tissue <- if (input$tissue != "") input$tissue else NULL
      
      # Call the smplSizPred function
      smplSizPred(
        df4model = df4modelInpt(),
        hIndexMeanDFinput = hIndexMeanDF(),
        heritabilityClass = heritabilityClass,
        inptPwr = inptPwr,
        fc = fc,
        trait = trait,
        tissue = tissue
      )
    })
  })
  
  # Render the sample size result
  # output$sampleSize <- renderText({
  #   req(sampleSize())  # Ensure sampleSize has been calculated
  #   paste(sampleSize())
  # })
  output$sampleSize <- renderUI({
    req(sampleSize())  # Ensure sampleSize has been calculated
    
    # Use HTML() to output formatted HTML content
    HTML(
      paste0(
        # First span with background color and curved corners on the right
        "<span style='font-size: 30px; background-color: #F0F0F0; padding: 5px 10px; border-radius: 10px 0 0 10px;'><b>Estimated Sample Size: </b></span>",
        
        # Second span with a different background color and curved corners on the left
        "<span style='color: #FFFFFF; background-color: #990F4B; font-size: 30px; padding: 5px 10px; border-radius: 0 10px 10px 0;'><b>", sampleSize(), "</b></span><br>"
      )
    )
    
  })
  
  # # Dynamically render the UI based on sample size availability
  # output$outPanel <- renderUI({
  #   if(input$submit > 0) {  # If submit button has been clicked
  #     tagList(
  #       # h3("Estimated Sample Size:"),
  #       # verbatimTextOutput("sampleSize")
  #       tags$h3("Result"),
  #       hr(),
  #       uiOutput("sampleSize"),
  #       uiOutput("outPanelSummary")
  #     )
  #   } else {
  #     h3("Execute the tool to estimate the Sample Size")
  #   }
  # })
  
  # Render the df4model
  output$df4model <- renderDT(
    df4modelInpt(), options = list(pageLength = 5),
    callback = JS(
      "$(function(){
      $('#outPanelSummary thead th').css({
        'background-color': '#990F4B',
        'color': 'white'
      });
    });"
    )
  )
  
  # Render the hIndexMeanDF
  output$hIndexMeanDF <- renderDT(
    hIndexMeanDF(), options = list(pageLength = 5),
    callback = JS(
      "$(function(){
      $('#outPanelSummary thead th').css({
        'background-color': '#990F4B',
        'color': 'white'
      });
    });"
    )
  )
  
  # Represent the summary with hr and text
  output$outPanelSummary <- renderUI({
    if(input$submit > 0) {  # If submit button has been clicked
      tagList(
        hr(),  # Add horizontal rule for separation
        tags$h3("Summary"),
        hr(),
        tabsetPanel(
          tabPanel(
            "Model Plot", 
            plotOutput("modelCRplot"),
            tags$b("Generic Plots", style = "display: block; text-align: center;"),
            plotOutput("modelPlot")
          ),
          tabPanel(
            "Predicted Model Coefficients", 
            DTOutput("pred.model.coeff")
          ),
          tabPanel(
            "Input Data Frames", 
            fluidRow(column(
              width = 8,
              div(style = "width:100%; overflow-x: scroll; font-size: 12px; padding: 5px;",
                  tags$h4("Data frame for model (df4model):"),  # Bold text
                  br(),
                  DTOutput("df4model")  # Render the DataTable for df4model
              )
            ),
            # Second table (hIndexMeanDF)
            column(
              width = 4,
              div(style = "width:100%; overflow-x: scroll; font-size: 12px; padding: 5px;",
                  tags$h4("Heritability Index Mean (hIndexMeanDF):"),  # Bold text
                  br(),
                  DTOutput("hIndexMeanDF")  # Render the DataTable for hIndexMeanDF
              )
            )
            ),
          )
        )
      )
    } 
  })
  
  # Logic for Power Calculation Tool
  observeEvent(input$submitPowerCalc, {
    # Check if both files are uploaded
    if (is.null(input$countDat) || is.null(input$smplDat)) {
      showNotification("Please upload both count and sample data files.", type = "error")
      return()
    }
    
    # Read the uploaded files
    # countDat <- read.csv(input$countDat$datapath, row.names = 1)
    # smplDat <- read.csv(input$smplDat$datapath, row.names = 1)
    
    countDat <- read.csv(input$countDat$datapath)
    smplDat <- read.csv(input$smplDat$datapath)
    
    # Collect the inputs for powerCalc function
    alpha_power <- input$alpha_power
    thrsholdFC <- input$thrsholdFC
    inptNoOfReplicates <- input$inptNoOfReplicates
    sims <- input$sims
    # saveFltrd <- input$saveFltrd
    # Show progress bar while running powerCalc
    withProgress(message = 'Running power calculation...', value = 0, {
      
      # Step 1: Estimating size factors
      Sys.sleep(0.5)  # Simulating a task (adjust the time as per actual task duration)
      incProgress(1/12, detail = "Estimating size factors...")
      
      # Step 2: Estimating dispersions
      Sys.sleep(0.5)  # Simulating a task
      incProgress(2/12, detail = "Estimating dispersions...")
      
      # Step 3: Gene-wise dispersion estimates
      Sys.sleep(0.5)  # Simulating a task
      incProgress(9/12, detail = "Gene-wise dispersion estimates...")
      
      # Step 4: Mean-dispersion relationship
      Sys.sleep(0.5)  # Simulating a task
      incProgress(10/12, detail = "Mean-dispersion relationship...")
      
      # Step 5: Final dispersion estimates
      Sys.sleep(0.5)  # Simulating a task
      incProgress(11/12, detail = "Final dispersion estimates...")
      
      # Step 6: Fitting model and testing
      Sys.sleep(0.5)  # Simulating a task
      incProgress(12/12, detail = "Fitting model and testing...")
      
      # Call the powerCalc function (replace with actual function logic)
      # resPowerCalc <- powerCalc(countDat, smplDat, alpha_power, thrsholdFC, inptNoOfReplicates, saveFltrd, sims)
      resPowerCalc <- powerCalc(countDat, smplDat, alpha_power, thrsholdFC, inptNoOfReplicates, sims)
    })
    
    
    # Display results in a DataTable with a message above the table
    output$outPanelPowerCalc <- renderUI({
      tagList(
        hr(),  # Horizontal line on top
        h4("The power calculation provided below output:"),  # Display the message
        div(style = "display: flex; align-items: center; width: 100%; gap: 20px;",  # Create a single row
            div(style = "width: 70%;",  
                DTOutput("powerCalcTable")  # Table in first column
            ),
            div(style = "width: 30%; text-align: center;",  
                downloadButton("downloadFilteredDEGs", "Download Filtered DEGs")  # Download button in second column
            )
        )
      )
    })
    
    
    # Modify the dataframe to have Na instead of no value or NaN value
    resPowerCalcPowerResults <- resPowerCalc$PowerResults
    # print(is.na(resPowerCalcPowerResults$Value))
    # Convert to numeric to ensure proper replacement of NaN, NA, Inf
    resPowerCalcPowerResults$Value <- as.numeric(resPowerCalcPowerResults$Value)
    resPowerCalcPowerResults$Value[is.na(resPowerCalcPowerResults$Value) | is.nan(resPowerCalcPowerResults$Value) | is.infinite(resPowerCalcPowerResults$Value)] <- NA
    # print(resPowerCalcPowerResults)
    # resPowerCalcPowerResults$Value <- as.character(resPowerCalcPowerResults$Value)
    resPowerCalcPowerResults[is.na(resPowerCalcPowerResults)] <- "NaN"
    # Render the DataTable output with simplified options
    output$powerCalcTable <- renderDT({
      datatable(
        # resPowerCalc$PowerResults,
        resPowerCalcPowerResults,
        options = list(
          autoWidth = TRUE,  # Adjust column widths automatically
          scrollX = TRUE,    # Enable horizontal scrolling if needed
          dom = 't'          # Only show the table without any other controls (like search, pagination, etc.)
        ),
        escape = FALSE  # Prevent escaping of special characters
      )
    })
    # print(head(resPowerCalc$FilteredDEGs))
    output$downloadFilteredDEGs <- downloadHandler(
      filename = function() {
        "FilteredDEGs.csv"  # Set file name
      },
      content = function(file) {
        req(resPowerCalc$FilteredDEGs)  # Ensure data exists
        write.csv(resPowerCalc$FilteredDEGs, file, row.names = FALSE)
      }
    )
    
    
    
    
    # Close the overlay after calculation
    #shinyjs::hide("powerCalcOverlay")
  })
  
  
  # Dynamically render the UI based on sample size availability
  output$outPanel <- renderUI({
    if(input$submit > 0) {  # If submit button has been clicked
      tagList(
        # h3("Estimated Sample Size:"),
        # verbatimTextOutput("sampleSize")
        tags$h3("Result"),
        hr(),
        uiOutput("sampleSize"),
        uiOutput("outPanelSummary")
      )
    } 
    else {
      h3("Execute the tool to estimate the Sample Size")
    }
  })
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)


