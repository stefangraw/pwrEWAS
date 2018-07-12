
#' @title Shiny pwrEWAS
#'
#' @description pwrEWAS.shiny provides a user-friendly point-and-click interface for pwrEWAS
#' 
#' @keywords DNAm microarray power Shiny
#' 
#' @export


pwrEWAS.shiny = function(){
  # setwd("S:\\Biostats\\BIO-STAT\\Koestler Devin\\GRAs\\Stefan Graw\\dissertation paper 1\\rcode\\pwrEWAS")
  # library(shiny)
  # library(shinyBS)
  # library(ggplot2)
  # library(parallel)
  
  # user input / default values
  input2 = NULL
  input2$Nmin = 10
  input2$Nmax = 50
  input2$NCntPer = 0.5
  input2$Nsteps = 10
  input2$J = 100000 # simulated CPGs
  input2$targetDmCpGs = 100
  input2$targetDeltaString = "0.2, 0.5"
  input2$targetDelta = as.numeric(unlist(strsplit(input2$targetDeltaString,",")))
  input2$method = "limma"
  input2$detectionLimit = 0.01
  input2$FDRcritVal = 0.05
  input2$cores = round(parallel::detectCores(all.tests = FALSE, logical = TRUE)/2)
  input2$sim = 50
  input2$tissueType = "Saliva"
  # input = input2
  
  
  
  server = function(input,output){
    
    shiny::observeEvent(input$goButton, {
      # eventReactive(input$goButton, {
      # runProgram <- 
      # eventReactive(input$goButton, {
      
      ### here goes all my code
      # source(file = "pwrEWAS_v1.5.R")
      shiny::withProgress(message = 'Program running.', detail = "Please wait!", value = NULL, {
        out = pwrEWAS(minTotSampleSize = input$Nmin,
                      maxTotSampleSize = input$Nmax,
                      SampleSizeSteps = input$Nsteps,
                      NcntPer = input$NCntPer,
                      # targetDelta = input$targetDelta,
                      targetDelta = as.numeric(unlist(strsplit(input$targetDeltaString,","))),
                      J = input$J,
                      targetDmCpGs = input$targetDmCpGs,
                      tissueType = input$tissueType,
                      detectionLimit = input$detectionLimit,
                      DMmethod = input$method,
                      FDRcritVal = input$FDRcritVal,
                      core = input$cores,
                      sims = input$sim)
        
        # source(file = "plotFunction_v1.3.R")
        output$powerPlot <- shiny::renderPlot({myPlotCI3D(out$powerArray)})
        
        meanPowerTable = cbind(rownames(out$meanPower), round(out$meanPower, 2))
        colnames(meanPowerTable)[1] = shiny::HTML("N</sub> \\ &Delta;<sub>&beta;")
        output$meanPower <- shiny::renderTable({meanPowerTable}, sanitize.text.function = function(x) x)
        output$deltaDensity = shiny::renderPlot({myDensityPlots(out$deltaArray, input$detectionLimit)})
      }) # processbar done
    })
    
  }
  
  
  ui = shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(".shiny-notification {
                            height: 70px;
                            width: 400px;
                            position:fixed;
                            font-size: 200%;
                            top: calc(50% - 35px);;
                            left: calc(50% - 100px);;}"))),
    
    shiny::titlePanel("pwrEWAS"),
    shiny::HTML("pwrEWAS is a computationally efficient tool to estimate power in EWAS as a function of sample and effect size 
       for two-group comparisons of DNAm (e.g., case vs control, exposed vs non-exposed, etc.). Detailed description 
       of in-/outputs, instructions and an example, as well as interpretations of the example results are provided in 
       the following vignette: "),
    shiny::tags$a(href="https://github.com/stefangraw/pwrEWAS/blob/master/vignette.pdf", "pwrEWAS vignette"),
    
    
    shiny::HTML("</br></br>Authors: Stefan Graw, Devin Koestler </br>"),
    shiny::HTML("Department of Biostatistics, University of Kansas School of Medicine"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        ### Inputs
        shiny::selectInput(inputId = "tissueType", label = "Tissue Type", choices = c("Adult (PBMC)",
                                                                               "Saliva", 
                                                                               "Lymphoma",
                                                                               "Placenta",
                                                                               "Liver",
                                                                               "Colon",
                                                                               "Blood adult",
                                                                               "Blood 5 year olds",
                                                                               "Blood newborns",
                                                                               "Cord-blood (whole blood)",
                                                                               "Cord-blood (PBMC)")),
        shinyBS::popify(shiny::numericInput(inputId = "Nmin", label = "Minimum total sample size", value = input2$Nmin, min = 4, step = 1),
               'Lower boundary of the range of total sample sizes.'),
        shinyBS::popify(shiny::numericInput(inputId = "Nmax", label = "Maximum total sample size", value = input2$Nmax, min = 4, step = 1),
               'Upper boundary of the range of total sample sizes.'),     
        shinyBS::popify(shiny::numericInput(inputId = "Nsteps", label = "Sample size increments", value = input2$Nsteps, min = 1, step = 1),
               'Sample size increments from lower to upper boundary of total sample sizes.'),     
        shinyBS::popify(shiny::numericInput(inputId = "NCntPer", label = "Percentage samples in group 1", value = input2$NCntPer, min = 0, max = 1, step = 0.1),
               'Percentage of the total sample split into groups (0.5 corresponds to a balanced study)'),     
        
        shiny::numericInput(inputId = "J", label = "Number of CpGs to be tested", value = input2$J, min = 1, step = 10000),
        shiny::numericInput(inputId = "targetDmCpGs", label = "Target number of DM CpGs", value = input2$targetDmCpGs, min = 1, step = 10),
        
        shiny::textInput(inputId = "targetDeltaString", label = "List of target max DM (comma delimited)", value = input2$targetDeltaString),
        shiny::numericInput(inputId = "FDRcritVal", label = "Target FDR", value = input2$FDRcritVal, min = 0, max = 1, step = 0.01), 
        
        shiny::checkboxInput(inputId = "advancedSettings", label = "Advanced settings"),
        shiny::conditionalPanel( 
          condition = "input.advancedSettings == 1",
          shiny::numericInput(inputId = "detectionLimit", label = "Detection Limit", value = input2$detectionLimit, min = 0, max = 1, step = 0.01), # 0.01 teschendorff paper / 0.05 Conference
          shiny::selectInput(inputId = "method", label = "Method for DM analysis", choices = c("limma", "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc")),
          shiny::numericInput(inputId = "sim", label = "Number of simulated data sets", value = input2$sim, min = 1, step = 10), 
          shiny::numericInput(inputId = "cores", label = "Threads", value = input2$cores, min = 1, max = parallel::detectCores(all.tests = FALSE, logical = TRUE)-1, step = 1)
        ),
        
        # submitButton(text = "Simulate"),
        shiny::actionButton(inputId = "goButton", label = "Go!", width = '100%', style='font-size:150%')
      ),
      
      ### Outputs
      shiny::mainPanel(
        shiny::fluidRow(
          shiny::column(12, align="center",
                        shiny::plotOutput("powerPlot"),
                        shiny::br(),shiny::br(),shiny::br(),
                 shiny::tableOutput(outputId = "meanPower"),
                 shiny::plotOutput("deltaDensity")
          )
        )
      )
      
    )
  )
  
  
  shiny::shinyApp(ui = ui, server = server)
}







