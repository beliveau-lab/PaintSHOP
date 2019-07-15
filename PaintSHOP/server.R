# PaintSHOP server logic
# Elliot Hershberg
# July 15, 2019

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    # a reactive element that is a vector of the RefSeq IDs
    # if the user inputs them manually on the UI
    manual_accessions <- reactive({
        # split the manual comma separated input
        strsplit(input$refseq_manual, ", ")
    })
    
    #
    
    
    
    
    
    output$plot <- renderPlot({
        plot(cars, type=input$plotType)
    })
    
    output$summary <- renderPrint({
        summary(cars)
    })
    
    output$table <- DT::renderDataTable({
        DT::datatable(cars)
    })

})
