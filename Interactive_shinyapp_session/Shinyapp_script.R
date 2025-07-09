if (!require("shiny")) install.packages("shiny")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("BiocManager")) install.packages("BiocManager")
install.packages("htmltools")
install.packages("scatterplot3d",repo="http://cran.ma.imperial.ac.uk")

BiocManager::install("affy")
BiocManager::install("limma") 
BiocManager::install("mouse4302.db")
BiocManager::install("annotate")
install.packages("htmltools")

library(mouse4302.db)
library(affy) 
library(limma) 
library(annotate)
library(shiny)
library(pheatmap)
library(shinyjs)
library(scatterplot3d)


setwd("D:/NRP DTP (PhD)/Data/Functional Genomic Technologies/Assigments/FGT_my_shiny")

#load the data. 
# can't upload to git. Too big file.
load("D:/NRP DTP (PhD)/Data/Functional Genomic Technologies/Assigments/FGT_my_shiny/FGT_assignm.Rdata")


# Define UI for application that draws a histogram
#create 3 tabs for 3 different plots
ui <- navbarPage( title = "Gene expression data", #Tomlinson, S. R. FGT_T7_shiny_example. (2020).
                  tabPanel( "Heatmap", # Shiny from RStudio. Learn Shiny. (2017). Available at: https://shiny.rstudio.com/tutorial/. (Accessed: 8th April 2020)
                            
                            # title
                            titlePanel("Most differentially expressed genes"),
                            
                            # Sidebar with a slider input for number of bins, 
                            # a box to type the title,
                            # parameters: to show row names, different scaling
                            
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("font_row0",
                                            label = "Font size of rows:",
                                            min = 6,
                                            max = 16,
                                            value = 10),
                                sliderInput("font_col0",
                                            label= "Font size of columns:",
                                            min = 6,
                                            max = 16,
                                            value = 10),
                                sliderInput("sample_size0",
                                            label= "Select the number of genes displayed in the heatmap:",
                                            min = 5,
                                            max = 50,
                                            value = 25),
                                textInput("Title0", #Shiny from RStudio. Learn Shiny. (2017). Available at: https://shiny.rstudio.com/tutorial/. (Accessed: 8th April 2020)
                                          label = "Can specify the title",
                                          value= "Top 25 differentially expressed genes"),
                                checkboxInput("srownames0", "Show Row Names", FALSE),
                                radioButtons("norm0", "Scale by", choices=c("none","row","column"))
                                
                              ),
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                plotOutput("distPlot0",  width = "100%", height = "600px") # Narasimhan, R. Scale and size of plot in RStudio shiny. (2013). Available at: https://stackoverflow.com/questions/17838709/scale-and-size-of-plot-in-rstudio-shiny. (Accessed: 8th April 2020)
                              )
                            )
                  ),
                  
                  #2nd tab for other plot
                  tabPanel( "Heatmap Log scale",
                            #Tomlinson, S. R. FGT_T7_shiny_example. (2020).
                            # Application title
                            titlePanel("Most differentially expressed genes"),
                            
                            # Sidebar with a slider input for number of bins 
                            # a box to type the title,
                            # parameters: to show row names, different scaling
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("font_row",
                                            label = "Font size of rows:",
                                            min = 6,
                                            max = 16,
                                            value = 10),
                                sliderInput("font_col",
                                            label= "Font size of columns:",
                                            min = 6,
                                            max = 16,
                                            value = 10),
                                sliderInput("sample_size",
                                            label= "Select the number of genes displayed in the heatmap:",
                                            min = 5,
                                            max = 50,
                                            value = 25),
                                textInput("Title",
                                          label = "Can specify the title",
                                          value= "Top 25 differentially expressed genes"),
                                checkboxInput("srownames", "Show Row Names", FALSE),
                                radioButtons("norm", "Scale by", choices=c("none","row","column"))
                                
                              ),
                              
                              # Show a plot of the generated distribution
                              mainPanel(
                                plotOutput("distPlot1",  width = "80%", height = "500px")
                              )
                            )
                  ),
                  
                  #3rd tab for the last plot
                  tabPanel("Volcano Plot",
                           #Tomlinson, S. R. FGT_T7_shiny_example. (2020).
                           # Application title
                           titlePanel("The empirical Bayes moderated t-statistic"),
                           fluidRow(
                             column(2, 
                                    wellPanel(
                                      
                                      # Sidebar with a slider input for number of bins
                                      # a box to type the title, colour for gene names
                                      
                                      sliderInput("genes",
                                                  label = "Number of genes:",
                                                  min = 0,
                                                  max = 10000,
                                                  value = 0),
                                      sliderInput("thresold_y",
                                                  label = "The height of the threshold line:",
                                                  min = 0,
                                                  max = 7,
                                                  value = 1.301),
                                      sliderInput("thresold_x",
                                                  label = "The width of the threshold line:",
                                                  min = 0,
                                                  max = 4,
                                                  value = 1),
                                      sliderInput("thresold_thickness",
                                                  label = "The thickness of the threshold line:",
                                                  min = 0,
                                                  max = 10,
                                                  value = 4),
                                      textInput("gene_col",
                                                label = "The colour of gene names",
                                                value= "Blue"),
                                      textInput("Title_v",
                                                label = "Can specify the title",
                                                value= "Volcano plot of differential gene expression of mouse primary aortic endothelial cells in control vs transgenic V290M mutants")
                                      
                                      
                                      
                                    )
                             ),
                             
                             # Show a plot of the generated distribution
                             
                             column(10, class="well",
                                    h4("Brush and double-click to zoom"), #RStudio. Plot interaction - zoom. Available at: http://www.baoruidata.com/examples/105-plot-interaction-zoom/. (Accessed: 3rd April 2020)
                                    #mainPanel(
                                    plotOutput("distPlot2",width = "100%", height = "800", 
                                               dblclick = "distPlot2_dbclick",
                                               brush=brushOpts(
                                                 id="distPlot2_brush",
                                                 resetOnNew = TRUE
                                               )
                                    )
                             )
                           )
                  )
)

#Tomlinson, S. R. FGT_T7_shiny_example. (2020).
# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  output$distPlot0 <- renderPlot({
    # loads data
    load("FGT_assignm.Rdata")
    
    # plot the heatmap
    pheatmap(expression_matrix_unique[1:input$sample_size0,],
             fontsize_row = input$font_row0,
             fontsize_col=input$font_col0,
             show_rownames=input$srownames0, 
             main= input$Title0, 
             cluster_col=FALSE,
             scale=input$norm0)
  }, execOnResize = F)
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  
  
  output$distPlot1 <- renderPlot({
    #loads data
    load("FGT_assignm.Rdata")
    # plots the heatmap at Log2 scale
    pheatmap(expression_matrix_unique_log2[1:input$sample_size,],
             fontsize_row = input$font_row,
             fontsize_col=input$font_col,
             show_rownames=input$srownames, 
             main= input$Title, 
             cluster_col=FALSE,
             scale=input$norm)
  }, execOnResize = F)
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  # zoomable plot 
  ranges <- reactiveValues(x=NULL, y=NULL)
  
  output$distPlot2 <- renderPlot({
    #loads data
    load("FGT_assignm.Rdata")
    # plots the volcano plot
    volcanoplot(fit2,
                coef =1,
                style = "p-value",
                highlight=input$genes,
                names=fit2$genes$Symbol,     #input$gene_names,
                hl.col=input$gene_col, 
                xlab = "Log2 Fold Change",
                ylab = "-Log10 p-value",
                pch=20,
                cex=1,
                main= input$Title_v,
                xlim=ranges$x,
                ylim= ranges$y)
    abline(h=input$thresold_y,v=c( -input$thresold_x,input$thresold_x), lty=2, lwd=input$thresold_thickness, col= "red")
  }, execOnResize = F)
  
  #RStudio. Plot interaction - zoom. Available at: http://www.baoruidata.com/examples/105-plot-interaction-zoom/. (Accessed: 3rd April 2020)
  observeEvent(input$refresh, { 
    session$invalidate
  })
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  # setting x and y coordinates, depending on the user's action
  observeEvent(input$distPlot2_dbclick, {
    brush <- input$distPlot2_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
