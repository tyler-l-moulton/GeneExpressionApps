#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plyr)
library(ggplot2)
library(shinyalert)
# Designer Functions #
colvec <- c("grey50", "dodgerblue", "coral", "gold", "aquamarine","hotpink", "green", "purple")
thresh.cross<-function(threshold, dat){
    
    datspline <- as.data.frame(spline(dat$log10.fluor))
    
    xy2 <- datspline[datspline[,2] >= threshold,][1,]           #above cross
    xy1 <- tail(datspline[datspline[,2] <= threshold,],n = 1)   #below cross
    
    xy21 <- (xy2-xy1)
    M <- xy21[2]/xy21[1]
    B <- xy2[2]-M*xy2[1]
    
    ThreshX <- (threshold-B)/M
    
    return(ThreshX)
    
}

findCTs <- function(Threshold, Data){
    
    CTframe <- NULL
    for(g in 1:length(unique(Data$treat))){
        dtreat <- Data[Data$treat == unique(Data$treat)[g],]
        
        for(h in 1:length(unique(dtreat$gene))){
            dgene <- dtreat[dtreat$gene == unique(dtreat$gene)[h],]
            
            for(i in 1:length(unique(dgene$SUBJECT))){
                dsubj <- dgene[dgene$SUBJECT == unique(dgene$SUBJECT)[i],]
                
                for(j in 1:length(unique(dsubj$REPLICATE))){
                    drep <- dsubj[dsubj$REPLICATE == unique(dsubj$REPLICATE)[j],]
                    
                    ct <- thresh.cross(threshold = Threshold, dat = drep)
                    
                    CTframe <- rbind(CTframe, 
                                     cbind(
                                           as.character(drep$treat[1]),
                                           as.character(drep$gene[1]),
                                           drep$SUBJECT[1], drep$REPLICATE[1],
                                           ct, Threshold)) #drep$subj.idx[1]
                    # print(CTframe)
                    
                    
                }
            }
        }
    }
    
    CTframe <- as.data.frame(CTframe)
    
    names(CTframe) <- c("treat", "gene", "SUBJECT", "REPLICATE", "ct", "Threshold") #"subj.idx", 
    
    CTframe$treat <- as.character(CTframe$treat)
    CTframe$gene <- as.character(CTframe$gene)
    
    return(CTframe)
    
    
}

tryObserve <- function(x) {
  x <- substitute(x)
  env <- parent.frame()
  observe({
    tryCatch(
      eval(x, env),
      error = function(e) {
        showNotification(paste("Error: ", e$message), type = "error")
      }
    )
  })
}
## try observe by https://community.rstudio.com/t/prevent-observer-crashes/25169/2 


# `%then%` <- shiny:::`%OR%`
# Define UI for application that draws a histogram
ui <- fluidPage(
useShinyalert(),
    # Application title
    titlePanel("Analyze qPCR Data"),

    # Sidebar with a slider input for number of bins 
    verticalLayout(
        sidebarLayout(
            sidebarPanel(
                verticalLayout(

                  fileInput(inputId = "file.name",
                            label = "Choose CSV file",
                            multiple = FALSE,
                            accept = ".csv"),
                  
                  selectInput(inputId = "SUBJECT", label = "Subject or Individual", ""),
                  
                  selectInput(inputId = "REPLICATE", label = "Replicate or Trial", ""),
                  
                  selectInput(inputId = "log10.fluor", label = "Log10 Fluorescence", ""),
                  
                  selectInput(inputId = "treat", label = "Treatment (or species)", ""),
                  
                  selectInput(inputId = "gene", label = "Genes", ""),
                  
                  selectInput(inputId = "cycle", label = "Cycles", ""),
                  
                  sliderInput(inputId = "thresh", label = "Select Threshold",
                              min = -10, max = 10, step = .01, value = 2),
                  
                  actionButton(inputId = "update", label = "Calculate CTs")
                                    
                ),
            
                hr(),
                verticalLayout(
                  selectInput(inputId = "RefGene", label = "Reference Gene", ""),
                  
                  selectInput(inputId = "ctrl", label = "Baseline Treatment Value", "")
                  
                )
                
                
            ),
            
            
            # Show a plot of the generated distribution
            mainPanel(
                verticalLayout(
                    tableOutput("thetable1"),
                    tableOutput("thetable2"),
                    plotOutput("theplot1")
                )
              )
        ),
        sidebarLayout(
            sidebarPanel(verticalLayout(
                

                checkboxInput(inputId = "subj.grp", label= "Average by Subject"),
                
                actionButton(inputId = "makeplot", label = "Build Plot")
                
            ),
                
                hr(),
                
            verticalLayout(
                
                textInput(inputId = "ct.download.name", label = "CT Data Filename"),
                
                downloadButton( outputId = "ct.download",
                                label = "Download CT Data")
                
            )   
            ),
            mainPanel(
                tableOutput("to1"),
                
            )
        )
    ),
    sidebarLayout(
      sidebarPanel(
        actionButton(inputId = "fc.button", label = "Calculate Relative Expression"),
        hr(),
        textInput(inputId = "fc.download.name", label = "Relative Expression Filename"),
        hr(),
        downloadButton( outputId = "fc.download",
                        label = "Download Fold Change Data")
      ),
      mainPanel(
        plotOutput("heatplot")
      )
      )


)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  ##Establish Reactive Values   
  v <- reactiveValues(the.dat = NULL,
                         ctdata = NULL,
                         graph.data = NULL,
                         fcdata = NULL)


 # CTDATA <- reactive({
     UDATA <- reactive({
         req(input$file.name)
         
       #Upload Data#  
       udata <- read.csv(file = input$file.name$datapath,
                           header = TRUE)
        
        ## Select variables from dropdown menu ##
       
         updateSelectInput(session = session, inputId = "SUBJECT",
                           label = "Subject or Individual",
                           choices = names(udata), selected =names(udata))
         
         updateSelectInput(session = session, inputId = "REPLICATE",
                           label = "Replicate or Trial",
                           choices = names(udata), selected =names(udata))
         
         updateSelectInput(session = session, inputId = "log10.fluor",
                           label = "Log10 Fluorescence",
                           choices = names(udata), selected =names(udata))
         
         updateSelectInput(session = session, inputId = "treat",
                           label = "Treatment (or species)",
                           choices = names(udata), selected =names(udata))
         
         updateSelectInput(session = session, inputId = "gene",
                           label = "Genes",
                           choices = names(udata), selected =names(udata))
         
         updateSelectInput(session = session, inputId = "cycle",
                           label = "Cycles",
                           choices = names(udata), selected =names(udata))
         


         return(udata)
         

     })

     ##At this point , the data.frame UDATA() has been created



     #Output table of uploaded data
     output$thetable1 <- renderTable(head(UDATA()))

     ## The following  is implemented every time "calculate CTs" is clicked
     observeEvent(eventExpr = input$update,{

       
       
       var.set.1 <- c(input$SUBJECT,
                      input$REPLICATE,
                      input$log10.fluor,
                      input$treat,
                      input$gene,
                      input$cycle)
       

       if(is.null(UDATA()) == T){
         
         
       }else if(length(unique(var.set.1)) != 6){
          
          shinyalert(title = "Input Error", 
                     text =    "Before clicking 'calculate CTs', Inputs must be 6 unique variable names",
                     type = "error")
      
        
          }else if(is.character(UDATA()[,input$REPLICATE]) == T){
            
            shinyalert(title = "Input Error", 
                       text =    "'Replicate or Trial' must be of class numeric",
                       type = "error")
            
          }else if(is.numeric(UDATA()[,input$log10.fluor]) == F){
            
            shinyalert(title = "Input Error", 
                       text =    "'Log10 Fluorescence' must be of class numeric",
                       type = "error")
            
          }else if(is.numeric(UDATA()[,input$cycle]) == F ||
                   is.integer(UDATA()[,input$cycle]) == F){
            
            shinyalert(title = "Input Error", 
                       text =    "'Cycles' must be of class numeric or integer",
                       type = "error")
            
          }else{
      

          
          
          len <- length(UDATA()[,1])
          
          the.dat <- data.frame(treat = character(len),
                                gene = character(len),
                                SUBJECT = integer(len),
                                REPLICATE = integer(len),
                                log10.fluor = numeric(len),
                                cycle = numeric(len))
          
          the.dat$treat <- as.character(UDATA()[,input$treat])
          the.dat$gene <- as.character(UDATA()[,input$gene])
          the.dat$SUBJECT <- UDATA()[,input$SUBJECT]
          the.dat$REPLICATE <- UDATA()[,input$REPLICATE]
          the.dat$log10.fluor <- UDATA()[,input$log10.fluor]
          the.dat$cycle <- UDATA()[,input$cycle]
          
          
          ## create an index for the internal looping of findCTs function ##
          create.idx <-  rep(c(1:(length(unique(UDATA()[,input$treat]))*
                                    length(unique(UDATA()[,input$gene]))*
                                    length(unique(UDATA()[,input$SUBJECT]))*
                                    length(unique(UDATA()[,input$REPLICATE])))),
                             each = max(UDATA()[,input$cycle]))
          
          if(length(create.idx) != length(the.dat[,1])){
            
            shinyalert(title = "Input Error", 
                       text = "Check that appropriate variable names are selected.",
                       type = "error")
            
          }else{
            the.dat$idx <- create.idx
          }
          
          
          ## make v$the.dat -- This data.frame is evaluated by function findCTs()
          v$the.dat <- the.dat
          
         ## Run findCTs function. If it crashes, generate an error message rather than crash app ##
          v$ctdata <- try(findCTs(Data = the.dat, Threshold = input$thresh),
                          silent = T)
          if(class(v$ctdata) == "try-error"){
            shinyalert(title = "Input Error", 
                       text = "findCTs() crashed. Check that appropriate variable names are selected.",
                       type = "error")
            
          }else{
            
            ref.cts <- v$ctdata[v$ctdata$gene == input$RefGene & 
                                  v$ctdata$treat == input$ctrl, "ct"]
            
            output$thetable2 <- renderTable(head(v$ctdata))

          }
          
          ## Allow user to select reference gene and baseline treatment ##
          ## This is congingent upon the ability of findCTs() to run
          
          updateSelectInput(session=session, inputId = "RefGene",
                            label = "Reference Gene",
                            choices = levels(UDATA()[,input$gene]))
          
          updateSelectInput(session=session, inputId = "ctrl",
                            label = "Baseline Treatment Value",
                            choices = levels(UDATA()[,input$treat]))
          
        }

  })
    

    
    observeEvent(eventExpr = input$makeplot,{
              if(is.null(v$the.dat) == T){
                
                shinyalert(title = "Error", 
                           text = "Calculate CTs before making plot!",
                           type = "error")
                
              }else if(input$subj.grp == T){
                        
                        v$graph.data <- v$the.dat[,-c(8)]
                        
                        # print(str(v$graph.data))
                        v$graph.data<- ddply(.data = v$graph.data,
                                             .variables = .(treat, gene,SUBJECT, cycle),                                                         
                                             summarize,
                                             log10.fluor = mean(x = log10.fluor, na.rm = TRUE ),
                                             idx = min(idx))
                          

                        
                        v$graph.data$treat <- as.factor(as.character(v$graph.data$treat))
                        v$graph.data$gene <- as.factor(as.character(v$graph.data$gene))
                        # print(str(v$graph.data))
                                                 
                    }else{
                        v$graph.data <- na.omit(object = v$the.dat)
                        v$graph.data$gene <- as.factor(v$graph.data$gene)
                        v$graph.data$treat <- as.factor(v$graph.data$treat)
                        
                        }
                     
      
                     output$theplot1<-renderPlot({
                         #Blank plot
                       
                       
                       plot(0, xlim = c(0,max(v$graph.data$cycle)),
                              ylim = c(min(v$graph.data$log10.fluor),
                                       max(v$graph.data$log10.fluor)),
                             xlab = "Cycle Number",
                             ylab = "Fluorescence (log10 transformed)",
                            type = "n")
                         
                         #Draw lines on plot

                         for(i in 1:length(unique(v$graph.data$idx))){
                             
                             datidx <- v$graph.data[v$graph.data$idx == unique(v$graph.data$idx)[i],]
                             
                             if(datidx$gene[1] == as.character(input$RefGene)){
                             ## make the reference gene  grey
                                     the.color <- colvec[1]
                             
                             }else{
                                         the.color <- colvec[as.numeric(datidx$gene[1])+1]
                                    }
                             
                             linetype <- datidx$treat[1]
                             
                             points(datidx$cycle, datidx$log10.fluor,
                                    type = "l", col = the.color, lty = as.numeric(linetype))
                             
                         }
              ## Prepare to make legend very inefficiently
                       legend.frame <- as.data.frame(unique(cbind(as.character(v$graph.data$gene),
                                                                  as.character(v$graph.data$treat))))
                       names(legend.frame) <- c("gene", "treat")
                       legend.frame$gene <-as.factor(legend.frame$gene)
                       legend.frame$treat <-as.factor(legend.frame$treat)
                       

                       the.colors <- character()
                       the.types <- character()
                       the.labs <- character()
                       
                       for(i in 1:length(legend.frame[,1])){
                         
                         
                         
                         if(legend.frame$gene[i] == as.character(input$RefGene)){
                           ## make the reference gene  grey
                           the.colors[i] <- colvec[1]
                           
                         }else{
                           the.colors[i] <- colvec[as.numeric(legend.frame$gene[i])+1]
                         }
                         
                          the.types[i] <- legend.frame$treat[i]
                         
                         the.labs[i] <- paste(legend.frame$treat[i],
                                              legend.frame$gene[i],
                                              sep = ".")
                         
                       }
                       
                       the.types <- as.factor(the.types)
                  
                       ## add legend to plot
                       legend(x = 1, y =max(v$graph.data$log10.fluor),
                              legend = the.labs,
                              col = the.colors,
                              lty = as.numeric(the.types))

                       
                       
                         abline(input$thresh,0, col = "gray")
                     },
                     width = 800, height = 800)
                     
                # Download Handler
                output$ct.download <- downloadHandler(
                     
                    filename = function(){
                         paste(input$ct.download.name, ".csv", sep="")
                    
                         },
                     content = function(file){
                         write.csv(x = v$ctdata, file = file, row.names=FALSE)
                     }
                 )
                         
    })

    
    observeEvent(eventExpr = input$fc.button,{
      
      ## Get mean CT of all technical replicates
      mn.ctdat<- ddply(.data = v$ctdata,
                       .variables = .(treat, gene,SUBJECT),                                                         
                       summarize,
                       mn.ct = mean(x = ct, na.rm = TRUE ))
      
      ##for each treatment, designate the reference gene
      ref.ct <- numeric()
      for(i in 1:length(unique(mn.ctdat$treat))){
        d.sub1 <- mn.ctdat[mn.ctdat$treat == unique(mn.ctdat$treat)[i],]
        refct <- rep(d.sub1[d.sub1$gene == input$RefGene, "mn.ct"], length(unique(d.sub1$gene)))
        ref.ct <- c(ref.ct, refct)
      }
      
      ## incorporate ref.ct into 'mn.ctdat' data.frame
      mn.ctdat$ref.ct <- ref.ct
      
      ## calculate delta CT (d.ct)
      mn.ctdat$d.ct <- mn.ctdat$mn.ct - mn.ctdat$ref.ct
      
      # Establish mean reference d.CT and incorporate into 'mn.ctdat' data.frame
      mn.ref.d.ct <- numeric()
      un.tg <- unique(cbind(mn.ctdat$treat, mn.ctdat$gene))
      for(i in 1:length(un.tg[,1])){

        subvec <-mn.ctdat[mn.ctdat$treat == input$ctrl &
                            mn.ctdat$gene == un.tg[i,2], "d.ct"]
        mrdc <- rep(mean(subvec), length(subvec))
        
        mn.ref.d.ct <- c(mn.ref.d.ct, mrdc)
      }
      mn.ctdat$mn.ref.d.ct <- mn.ref.d.ct 
      
      
      ## calculate delta delta CT, incorporate into 'mn.ctdat' data.frame
      mn.ctdat$dd.ct <- mn.ctdat$d.ct - mn.ctdat$mn.ref.d.ct
      
      ## calculate relative expression
      mn.ctdat$RelExpr <- 2^-(mn.ctdat$dd.ct)
      
      v$mn.ctdat <- mn.ctdat

      
      #Average by treatment and gene, except for reference gene
      
      sumct <- ddply(.data = v$mn.ctdat[v$mn.ctdat$gene != input$RefGene,],
                     .variables = .(treat, gene),                                                         
                     summarize,
                     FC = mean(x = RelExpr, na.rm = TRUE ))

      #Create heat plot
      output$heatplot <- renderPlot({

        ggplot(sumct , aes(x = gene, y = treat)) +
          geom_tile(aes(fill = FC)) +
          scale_fill_gradient2(low="cyan", mid="black", high="red", 
                               midpoint=1, limits=range(c(0,max(sumct$FC)))) +
          theme_classic()
        
        
      }, width = 100*length(sumct[,1]) +100, height = 100*length(sumct[1,])
      ) 
      
    })
    
    #Download fold change data table
    output$fc.download <- downloadHandler(
      
      filename = function(){
        paste(input$fc.download.name, ".csv", sep="")
        
      },
      content = function(file){
        write.csv(x = v$mn.ctdat, file = file, row.names=FALSE)
      }
    )


}

# Run the application 
shinyApp(ui = ui, server = server)
