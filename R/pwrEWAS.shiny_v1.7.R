
#' @title Shiny pwrEWAS
#'
#' @description pwrEWAS_shiny provides a user-friendly point-and-click interface for pwrEWAS
#' 
#' @keywords DNAm microarray power Shiny
#' 
#' @return pwrEWAS_shiny initializes pwrEWAS's user-interface
#' 
#' @export
#' 
#' @examples
#' 
#' if(interactive()) {
#'   pwrEWAS_shiny()
#' }


pwrEWAS_shiny <- function(){
    # library(shiny)
    # library(shinyBS)
    # library(ggplot2)
    # library(parallel)
    
    # user input / default values
    input2 <- NULL
    input2$Nmin <- 10
    input2$Nmax <- 50
    input2$NCntPer <- 0.5
    input2$Nsteps <- 10
    input2$J <- 100000 # simulated CPGs
    input2$targetDmCpGs <- 100
    input2$targetDeltaString <- "0.2, 0.5"
    input2$tauString <- "0.01, 0.03"
    input2$targetDelta <- as.numeric(unlist(strsplit(input2$targetDeltaString,",")))
    input2$method <- "limma"
    input2$detectionLimit <- 0.01
    input2$FDRcritVal <- 0.05
    input2$cores <- round(parallel::detectCores(all.tests = FALSE, logical = TRUE)/2)
    input2$sim <- 50
    input2$tissueType <- "Saliva"
    # input <- input2
    
    #############################################################
    
    server <- function(input,output){
        
        shiny::observeEvent(input$goButton, {
            
            # reset plots
            output$powerPlot <- NULL
            output$meanPower <- NULL
            output$probTP <- NULL
            output$deltaDensity <- NULL
            output$log <- NULL
            
            shiny::withProgress(message = 'Program running. Please wait.', detail = "This can take several minutes. Progress will be displayed in R console.", value = NULL, {
                
                runTimeStart <- Sys.time()
                if(input$switchTargetDmSd == 1){
                    out <- pwrEWAS(minTotSampleSize = input$Nmin,
                        maxTotSampleSize = input$Nmax,
                        SampleSizeSteps = input$Nsteps,
                        NcntPer = input$NCntPer,
                        targetDelta = as.numeric(unlist(strsplit(input$targetDeltaString,","))),
                        J = input$J,
                        targetDmCpGs = input$targetDmCpGs,
                        tissueType = input$tissueType,
                        detectionLimit = input$detectionLimit,
                        DMmethod = input$method,
                        FDRcritVal = input$FDRcritVal,
                        core = input$cores,
                        sims = input$sim)
                } else if(input$switchTargetDmSd == 2){
                    out <- pwrEWAS(minTotSampleSize = input$Nmin,
                        maxTotSampleSize = input$Nmax,
                        SampleSizeSteps = input$Nsteps,
                        NcntPer = input$NCntPer,
                        deltaSD = as.numeric(unlist(strsplit(input$tauString,","))),
                        J = input$J,
                        targetDmCpGs = input$targetDmCpGs,
                        tissueType = input$tissueType,
                        detectionLimit = input$detectionLimit,
                        DMmethod = input$method,
                        FDRcritVal = input$FDRcritVal,
                        core = input$cores,
                        sims = input$sim)
                }
                
                output$powerPlot <- shiny::renderPlot({isolate(pwrEWAS_powerPlot(out$powerArray, sd = ifelse(input$switchTargetDmSd == 1, FALSE, TRUE)))})
                
                # mean power table
                meanPowerTable <- cbind(rownames(out$meanPower), round(out$meanPower, 2))
                if(input$switchTargetDmSd == 1){
                    colnames(meanPowerTable)[1] <- shiny::HTML("N</sub> \\ &Delta;<sub>&beta;")
                } else if(input$switchTargetDmSd == 2){
                    colnames(meanPowerTable)[1] <- shiny::HTML("N</sub> \\ SD(&Delta;<sub>&beta;)")
                }
                positionToAddTitle <- ceiling(dim(meanPowerTable)[2]/2)
                colnames(meanPowerTable)[positionToAddTitle] <- paste0(shiny::HTML("Power<br/>"), colnames(meanPowerTable)[positionToAddTitle])
                output$meanPower <- shiny::renderTable({meanPowerTable}, sanitize.text.function = function(x) x)
                
                # delta density plot
                output$deltaDensity <- shiny::renderPlot({isolate(pwrEWAS_deltaDensity(out$deltaArray, input$detectionLimit, sd = ifelse(input$switchTargetDmSd == 1, FALSE, TRUE)))})
                
                # probability of detecting at least one TP
                probTPTable <- cbind(rownames(out$metric$probTP), round(out$metric$probTP, 2))
                if(input$switchTargetDmSd == 1){
                    colnames(probTPTable)[1] <- shiny::HTML("N</sub> \\ &Delta;<sub>&beta;")
                } else if(input$switchTargetDmSd == 2){
                    colnames(probTPTable)[1] <- shiny::HTML("N</sub> \\ SD(&Delta;<sub>&beta;)")
                }
                colnames(probTPTable)[positionToAddTitle] <- paste0(shiny::HTML("P(#TP&ge;1) <br/>"), colnames(probTPTable)[positionToAddTitle])
                output$probTP <- shiny::renderTable({probTPTable}, sanitize.text.function = function(x) x)
                
                # run time
                runTimeStop <- difftime(Sys.time(), runTimeStart, units = "auto") 
                
                # log 
                logString <- paste0(
                    "Tissue type = ", input$tissueType, "\n",
                    "Minimum total sample size = ", input$Nmin, "\n",
                    "Maximum total sample size = ", input$Nmax, "\n",
                    "Sample size increments = ", input$Nsteps, "\n",
                    "Percentage samples in group 1 = ", input$NCntPer, "\n",
                    "Number of CpGs to be tested = ", input$J, "\n",
                    "Target number of DM CpGs = ", input$targetDmCpGs, "\n",
                    if(input$switchTargetDmSd == 1){
                        paste0("'Target max Delta' was selected \n", 
                            "Target maximal difference in DNAm (comma delimited) = ", input$targetDeltaString)
                    } else if(input$switchTargetDmSd == 2){
                        paste0("'SD(&Delta;)Delta)' was selected \n", 
                            "Std. dev. of difference in DNAm (comma delimited) = ", input$tauString)}, "\n",
                    "Target FDR = ", input$FDRcritVal, "\n",
                    "Detection Limit = ", input$detectionLimit, "\n",
                    "Method for DM analysis = ", input$method, "\n",
                    "Number of simulated data sets = ", input$sim, "\n",
                    "Threads = ", input$cores, "\n",
                    "Run time = ", round(runTimeStop,1), " ", attr(runTimeStop, "units"))
                
                
                output$log <- renderText({HTML(logString)})
            }) # processbar done
        })
        
    }
    
    
    ui <- shiny::fluidPage(
        shiny::tags$head(shiny::tags$style(shiny::HTML(".shiny-notification {
            height: 150px;
            width: 400px;
            position:fixed;
            font-size: 200%;
            top: calc(50% - 35px);;
            left: calc(50% - 100px);;}"))),
        shiny::tags$style(type='text/css', '#log {text-align: left;}'),
        
        shiny::titlePanel("pwrEWAS"),
        shiny::HTML("pwrEWAS is a computationally efficient tool to estimate power in EWAS as a function of sample and effect size 
            for two-group comparisons of DNAm (e.g., case vs control, exposed vs non-exposed, etc.). Detailed description 
            of in-/outputs, instructions and an example, as well as interpretations of the example results are provided in 
            the following vignette: "),
        shiny::tags$a(href="https://github.com/stefangraw/pwrEWAS/blob/master/vignettes/vignette.pdf", "pwrEWAS vignette"),
        
        
        shiny::HTML("</br></br>Authors: Stefan Graw, Devin Koestler </br>"),
        shiny::HTML("Department of Biostatistics, University of Kansas School of Medicine"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                ### Inputs
                shinyBS::popify(shiny::selectInput(inputId = "tissueType", label = "Tissue Type", choices = c("Adult (PBMC)",
                        "Saliva", 
                        "Sperm", 
                        "Lymphoma",
                        "Placenta",
                        "Liver",
                        "Colon",
                        "Blood adult",
                        "Blood 5 year olds",
                        "Blood newborns",
                        "Cord-blood (whole blood)",
                        "Cord-blood (PBMC)")),
                    'Heterogeneity of different tissue types can have effects on the results. Please select your tissue type of interest or one you believe is the closest.', placement =  "top"),
                shinyBS::popify(shiny::numericInput(inputId = "Nmin", label = "Minimum total sample size", value = input2$Nmin, min = 4, step = 1),
                    'Lowest total sample sizes to be considered.'),
                shinyBS::popify(shiny::numericInput(inputId = "Nmax", label = "Maximum total sample size", value = input2$Nmax, min = 4, step = 1),
                    'Highest total sample sizes to be considered.'),     
                shinyBS::popify(shiny::numericInput(inputId = "Nsteps", label = "Sample size increments", value = input2$Nsteps, min = 1, step = 1),
                    'Steps with which total sample size increases from "Minimum total sample size" to "Maximum total sample size".'),     
                shinyBS::popify(shiny::numericInput(inputId = "NCntPer", label = "Samples rate for group 1", value = input2$NCntPer, min = 0, max = 1, step = 0.1),
                    'Rate by which the total sample size is split into groups (0.5 corresponds to a balanced study; rate for group 2 is equal to 1 rate of group 1)'),     
                
                shinyBS::popify(shiny::numericInput(inputId = "J", label = "Number of CpGs tested", value = input2$J, min = 1, step = 10000),
                    'Number of CpG site that will simulated and tested (increasing Number of CpGs tested will require increasing RAM (memory)).'),
                shinyBS::popify(shiny::numericInput(inputId = "targetDmCpGs", label = "Target number of DM CpGs", value = input2$targetDmCpGs, min = 1, step = 10),
                    'Target number of CpGs simulated with meaningful differences (differences greater than detection limit)'),
                
                shinyBS::popify(shinyWidgets::radioGroupButtons(inputId = "switchTargetDmSd",choiceValues = c(1,2), justified = TRUE, choiceNames = c(shiny::HTML("Target max &Delta;"), shiny::HTML("SD(&Delta;)"))),
                    shiny::HTML('The expected simulated differences in methylation can be control by "Target max &Delta;" or "SD(&Delta;)". For "Target max &Delta;" standard deviations of the simulated differences is automatically determined such that the 99%til of the simulated differences are within a range around the provided values. If "SD(&Delta;)" is chosen, differences in methylation will be simulated using provided standard deviation.')),
                
                
                shiny::conditionalPanel( 
                    condition = "input.switchTargetDmSd == 1",
                    shinyBS::popify(shiny::textInput(inputId = "targetDeltaString", label = "Target maximal difference in DNAm (comma delimited)", value = input2$targetDeltaString),
                        'Standard deviations of the simulated differences is automatically determined such that the 99%til of the simulated differences are within a range around the provided values.')
                ),
                
                shiny::conditionalPanel( 
                    condition = "input.switchTargetDmSd == 2",
                    shinyBS::popify(shiny::textInput(inputId = "tauString", label = "Std. dev. of difference in DNAm (comma delimited)", value = input2$tauString),
                        'Differnces in methylation will be simulated using provided standard deviation.')
                ),
                
                
                shinyBS::popify(shiny::numericInput(inputId = "FDRcritVal", label = "Target FDR", value = input2$FDRcritVal, min = 0, max = 1, step = 0.01),
                    'Critical value to control the False Discovery Rate (FDR) using the Benjamini and Hochberg method.'), 
                
                shiny::checkboxInput(inputId = "advancedSettings", label = "Advanced settings"),
                shiny::conditionalPanel( 
                    condition = "input.advancedSettings == 1",
                    shinyBS::popify(shiny::numericInput(inputId = "detectionLimit", label = "Detection Limit", value = input2$detectionLimit, min = 0, max = 1, step = 0.01),
                        'Limit to detect changes in methylation. Simulated differences below the detection limit will not be consider as meaningful differentially methylated CpGs.'),
                    shinyBS::popify(shiny::selectInput(inputId = "method", label = "Method for DM analysis", choices = c("limma", "t-test (unequal var)", "t-test (equal var)", "Wilcox rank sum", "CPGassoc")),
                        'Method used to perform differential methylation analysis.', placement =  "top"),
                    shinyBS::popify(shiny::numericInput(inputId = "sim", label = "Number of simulated data sets", value = input2$sim, min = 1, step = 10),
                        'Number of repeated simulation/simulated data sets under the same conditions for consistent results.'),
                    shinyBS::popify(shiny::numericInput(inputId = "cores", label = "Threads", value = input2$cores, min = 1, max = parallel::detectCores(all.tests = FALSE, logical = TRUE)-1, step = 1),
                        'Number of cores used to run multiple threads. Ideally, the number of different total samples sizes multiplied by the number of effect sizes should be a multiple (m) of the number of cores (#sampleSizes * #effectSizes = m * #threads). An increasing number of threads will require an increasing amount of RAM (memory).', placement =  "top")
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
                        shiny::splitLayout(cellWidths = c("50%", "50%"), 
                            shiny::tableOutput(outputId = "meanPower"),
                            shiny::tableOutput(outputId = "probTP")),
                        shiny::plotOutput("deltaDensity"),
                        shiny::verbatimTextOutput(outputId = "log")
                    )
                )
            )
            
        )
    )
    
    
    shiny::shinyApp(ui = ui, server = server)
}
# pwrEWAS_shiny()






