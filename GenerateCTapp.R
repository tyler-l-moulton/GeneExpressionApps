#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyalert)

log.gen <- function(ALPHA, BETA, GAMMA,
                    cycles, SD.pct=.01,
                    CTshift=0,
                    qPCR.trials, n.subj,
                    subj.CT.SD){
    
    cycle.vec <- rep(c(1:cycles), qPCR.trials * n.subj) 
    
    fluor.rand <- numeric(cycles * qPCR.trials * n.subj )  #empty vector for fluor readings, length = cycles 
    for(j in c(1:n.subj)){           #for each subject, randomly vary CT by subj.CT.SD
        fluor.curve <- (ALPHA / 
                            (1 + BETA*exp(-GAMMA *cycle.vec +
                                              rnorm(n=1, mean = CTshift, sd = subj.CT.SD))))
        for(k in c(1:qPCR.trials)){
            for(i in c(1:cycles)){
                
                #adds fluoresecence value for each cycle, trial, and subject sequentially
                fluor.rand[i+
                               cycles*((k-1) +
                                           (j-1)*qPCR.trials)] <- rnorm(n=1,
                                                            mean = fluor.curve[i],
                                                            sd = SD.pct*10/(
                                                                (ALPHA^2)*(10+
                                                                    (log(fluor.curve[i]+1))*100)))
                                                           
            }      
        }
    } 
    
    
    trial.vec <- rep(c(1:qPCR.trials), times = n.subj, each = cycles)
    subj.vec <-  rep(c(1:n.subj), each = qPCR.trials * cycles)
    
    gen.dat <- as.data.frame(cbind(subj.vec, trial.vec, cycle.vec, fluor.rand))
    names(gen.dat)<- c("subject", "qPCR.Trial", "cycle", "log10.fluor")
    return(gen.dat)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    useShinyalert(),
    
    # Application title
    titlePanel("Generate qPCR Data"),
    hr(HTML(paste("This app generates data for simulated qPCR experiments.", 
                  " Colors on  plots represent different subjects.", 
            " The logistic function is of the form Fluor = ", 
            tags$em("a"), " / ","(1+", tags$em("b*e"),
            tags$sup(tags$em("-y"), "*cycle"),")"))),
verticalLayout(
    sidebarLayout(
        sidebarPanel(
           
            titlePanel("Logistic Curve parameters"),
            
             # Sidebar with a slider input for number of bins 
            sliderInput(inputId = "ALPHA",
                        "Alpha",
                        min = .01,
                        max = 5,
                        value = 2,
                        step = .01),
            
            sliderInput(inputId = "BETA",
                        "Beta",
                        min = 100,
                        max = 10000,
                        value = 2500,
                        step = 100),
            
            sliderInput(inputId = "GAMMA",
                        "Gamma",
                        min = 0.01,
                        max = 1,
                        value = 0.45,
                        step = .001),
            

            titlePanel("qPCR Parameters"),
            
            numericInput(inputId = "cycles",
                         "Number of Cycles",
                         min = 5,
                         max = 60,
                         value = 40,
                         step = 1),
            

            numericInput(inputId = "qPCR.trials",
                         "qPCR Replicates",
                         min = 1,
                         max = 20,
                         value = 3,
                         step = 1),
            
            numericInput(inputId = "n.subj",
                         "Number of Subjects",
                         min = 1,
                         max = 100,
                         value = 10,
                         step = 1),
            
            numericInput(inputId = "CTshift",
                         "Shift CT",
                         min = -20,
                         max = 40,
                         value = 0,
                         step = .01),
            
            titlePanel("Adjust Variability"),
            
            sliderInput(inputId = "SD.pct",
                        "Noise",
                        min = 0.00,
                        max = 11,
                        value = 01,
                        step = .01),

            sliderInput(inputId = "subj.CT.SD",
                        "Standard Deviation of CT among Subjects",
                        min = .00,
                        max = 2,
                        value = .5,
                        step = .01),
            
            
            textInput(inputId = "treat.label",
                        "Treatment",
                        width = 200,
                        placeholder = "treat"),
            
            textInput(inputId = "gene.label",
                        "Gene",
                        width = 200,
                        placeholder = "gene"),
            
            
            actionButton(inputId = "UPDATE",
                         label = "Update Data")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("rawplot"),
            plotOutput("log10plot")
            
        )
    ),
    titlePanel("qPCR Data"),
    
    textInput(inputId = "dataname",
              label = "File Name", 
              width = 500,
              placeholder = "your_fname_here"),
    
    downloadButton( outputId = "d.download",
                 label = "Download Data"),
    

    tableOutput(outputId = "The.Data")

    
)
    

)

# Define server logic 
server <- function(input, output) {


v <- reactiveValues(the.data = NULL)
output$rawplot <- renderPlot({
    # generate bins based on input$bins from ui.R
   par(mfrow=c(2,1))
     plot(1,1, type = "n",
         main = "qPCR Fluorescence",
         ylab = "Fluorescence")

     plot(1,1, type = "n",
         main = "log10 qPCR Fluorescence",
         ylab = "log10 Fluorescence")
    
}, width = 1000, height = 1300)

    observeEvent(eventExpr = input$UPDATE,
            { 
                 if(input$treat.label == ""){
                     shinyalert(title = "Input Error", 
                                text =    "Before clicking 'Update Data', specify 'treatment'",
                                type = "error")
                   
                 }else if(input$gene.label == ""){
                     shinyalert(title = "Input Error", 
                                text =    "Before clicking 'Update Data', specify 'gene'",
                                type = "error")
                 }else{
                     v$the.data <-log.gen(input$ALPHA, input$BETA, input$GAMMA,
                                          input$cycles, input$SD.pct,
                                          input$CTshift,
                                          input$qPCR.trials, input$n.subj,
                                          input$subj.CT.SD)
                     v$the.data$fluor <- 10 ^v$the.data$log10.fluor
                     v$the.data$treatment <- rep(input$treat.label, length(v$the.data[,1]))
                     v$the.data$gene <-  rep(input$gene.label, length(v$the.data[,1]))
                     
                     
                     output$rawplot <- renderPlot({
                         par(mfrow=c(2,1))
                         # generate bins based on input$bins from ui.R
                         plot(fluor ~ cycle, data = v$the.data, type = "n",
                              main = "qPCR Fluorescence",
                              ylab = "Fluorescence")
                         
                         for(i in 1:length(unique(v$the.data$subject))){
                             d.sub1 <- v$the.data[v$the.data$subject == unique(v$the.data$subject)[i],]
                             for(j in c(1:length(unique(d.sub1$qPCR.Trial)))){
                                 d.sub2 <- d.sub1[d.sub1$qPCR.Trial == unique(d.sub1$qPCR.Trial)[j],]
                                 points(x=d.sub2$cycle, y = d.sub2$fluor, 
                                        col = d.sub2$subject[1], type = "l")
                             }
                         }
                         
                         
                         plot(log10.fluor ~ cycle, data = v$the.data, type = "n",
                              main = "Log 10 Transformed qPCR Fluorescence",
                              ylab = "log10 Fluorescence")
                         
                         for(i in 1:length(unique(v$the.data$subject))){
                             d.sub1 <- v$the.data[v$the.data$subject == unique(v$the.data$subject)[i],]
                             for(j in c(1:length(unique(d.sub1$qPCR.Trial)))){
                                 d.sub2 <- d.sub1[d.sub1$qPCR.Trial == unique(d.sub1$qPCR.Trial)[j],]
                                 points(x=d.sub2$cycle, y = d.sub2$log10.fluor, 
                                        col = d.sub2$subject[1], type = "l")
                             }
                         }
                     }, width = 1000, height = 1300)} 
                 }
             )
 
        output$The.Data <- renderTable({head(v$the.data,100)})
        
        output$d.download <- downloadHandler(
            filename = function(){
                if(input$dataname == ""){
                    shinyalert(title = "Download Error", 
                               text =    "Specify File Name",
                               type = "error")
                }else{
                    paste(input$dataname, ".csv", sep="")    
                }
                    
                
            },
            content = function(file){
                write.csv(v$the.data, file, row.names=FALSE)
            }
        )

    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
